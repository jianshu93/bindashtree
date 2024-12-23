use clap::{Arg, ArgAction, Command};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use needletail::{parse_fastx_file, Sequence};
use std::collections::HashMap;
use kmerutils::sketcharg::{SeqSketcherParams, SketchAlgo, DataType};
use kmerutils::base::{
    CompressedKmerT, KmerBuilder,
    kmergenerator::{KmerGenerator, KmerGenerationPattern},
    alphabet::Alphabet2b,
    sequence::Sequence as SequenceStruct,
    kmer::{Kmer32bit, Kmer16b32bit, Kmer64bit}
};
use kmerutils::sketching::setsketchert::*; // Contains SeqSketcherT, OptDensHashSketch, RevOptDensHashSketch
use anndists::dist::{Distance, DistHamming};
use num;
use std::path::Path;
use speedytree::DistanceMatrix;
use speedytree::{Canonical, Hybrid, NeighborJoiningSolver, RapidBtrees};
use std::str::FromStr;
use std::fmt::Debug;
use serde::Serialize;
use rand_distr::uniform::SampleUniform;

// Introduce SeqSketcherFactory trait to provide `new` method.
trait SeqSketcherFactory<Kmer>: SeqSketcherT<Kmer>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer>,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
{
    fn new(params: &SeqSketcherParams) -> Self;
}

// Implement SeqSketcherFactory for OptDensHashSketch
impl<Kmer, S> SeqSketcherFactory<Kmer> for OptDensHashSketch<Kmer, S>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    Kmer::Val: num::PrimInt + Send + Sync + Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    S: num::Float + SampleUniform + Send + Sync + Debug + Serialize,
{
    fn new(params: &SeqSketcherParams) -> Self {
        // Call the existing public new method from OptDensHashSketch
        OptDensHashSketch::<Kmer, S>::new(params)
    }
}

// Implement SeqSketcherFactory for RevOptDensHashSketch
impl<Kmer, S> SeqSketcherFactory<Kmer> for RevOptDensHashSketch<Kmer, S>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    Kmer::Val: num::PrimInt + Send + Sync + Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    S: num::Float + SampleUniform + Send + Sync + Debug + Serialize,
{
    fn new(params: &SeqSketcherParams) -> Self {
        // Call the existing public new method from RevOptDensHashSketch
        RevOptDensHashSketch::<Kmer, S>::new(params)
    }
}

#[derive(Debug, Clone)]
enum TreeAlgo {
    Naive,
    RapidNJ,
    Hybrid,
}

impl FromStr for TreeAlgo {
    type Err = String;
    fn from_str(s: &str) -> Result<TreeAlgo, String> {
        match s.to_lowercase().as_str() {
            "naive" => Ok(TreeAlgo::Naive),
            "rapidnj" => Ok(TreeAlgo::RapidNJ),
            "hybrid" => Ok(TreeAlgo::Hybrid),
            _ => Err(format!("Unknown tree method: {}", s)),
        }
    }
}

fn ascii_to_seq(bases: &[u8]) -> Result<SequenceStruct, ()> {
    let alphabet = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alphabet);
    Ok(seq)
}

fn read_sequences(path: &str) -> Vec<SequenceStruct> {
    let mut sequences = Vec::new();
    let mut reader = parse_fastx_file(path).expect("Invalid FASTA/Q file");
    while let Some(record) = reader.next() {
        let seq_record = record.expect("Error reading sequence record");
        let seq_seq = seq_record.normalize(false).into_owned();
        let seq = ascii_to_seq(&seq_seq).unwrap();
        sequences.push(seq);
    }
    sequences
}

fn sketch_with<Kmer, Sketcher>(
    sketch_args: &SeqSketcherParams,
    genomes: &Vec<String>,
    kmer_size: usize,
    dens: usize
) -> HashMap<String, Vec<f32>>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    Sketcher: SeqSketcherFactory<Kmer, Sig = f32> + Send + Sync + 'static,
{
    let nb_alphabet_bits = 2;
    let sketcher = Sketcher::new(sketch_args);
    let hash_fn = move |kmer: &Kmer| -> <Kmer as CompressedKmerT>::Val {
        let mask: <Kmer as CompressedKmerT>::Val =
            num::NumCast::from::<u64>((1u64 << (nb_alphabet_bits * kmer.get_nb_base())) - 1).unwrap();
        let canonical = kmer.reverse_complement().min(*kmer);
        canonical.get_compressed_value() & mask
    };

    genomes
        .par_iter()
        .map(|path| {
            let sequences = read_sequences(path);
            let sequences_ref: Vec<&SequenceStruct> = sequences.iter().collect();
            let signature = sketcher.sketch_compressedkmer_seqs(&sequences_ref, &hash_fn);
            (path.clone(), signature[0].clone())
        })
        .collect()
}

fn sketch_genomes(
    kmer_size: usize,
    dens: usize,
    sketch_args: &SeqSketcherParams,
    genomes: &Vec<String>
) -> HashMap<String, Vec<f32>> {
    if kmer_size <= 14 {
        if dens == 0 {
            sketch_with::<Kmer32bit, OptDensHashSketch<Kmer32bit, f32>>(sketch_args, genomes, kmer_size, dens)
        } else {
            sketch_with::<Kmer32bit, RevOptDensHashSketch<Kmer32bit, f32>>(sketch_args, genomes, kmer_size, dens)
        }
    } else if kmer_size == 16 {
        if dens == 0 {
            sketch_with::<Kmer16b32bit, OptDensHashSketch<Kmer16b32bit, f32>>(sketch_args, genomes, kmer_size, dens)
        } else {
            sketch_with::<Kmer16b32bit, RevOptDensHashSketch<Kmer16b32bit, f32>>(sketch_args, genomes, kmer_size, dens)
        }
    } else if kmer_size <= 32 {
        if dens == 0 {
            sketch_with::<Kmer64bit, OptDensHashSketch<Kmer64bit, f32>>(sketch_args, genomes, kmer_size, dens)
        } else {
            sketch_with::<Kmer64bit, RevOptDensHashSketch<Kmer64bit, f32>>(sketch_args, genomes, kmer_size, dens)
        }
    } else {
        panic!("kmers cannot be 15 or greater than 32");
    }
}

fn build_distance_matrix(
    sketches: &HashMap<String, Vec<f32>>,
    kmer_size: usize,
    genomes: &Vec<String>,
) -> Vec<u8> {
    let dist_hamming = DistHamming;
    let n = genomes.len();
    let distances: Vec<(usize, usize, f64)> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            (i + 1..n)
                .into_par_iter()
                .map(move |j| {
                    let query_signature = &sketches[&genomes[i]];
                    let reference_signature = &sketches[&genomes[j]];
                    let hamming_distance = dist_hamming.eval(query_signature, reference_signature);
                    let hamming_distance = if hamming_distance == 0.0 {
                        std::f32::EPSILON // Use a small value close to zero
                    } else {
                        hamming_distance
                    };
                    let jaccard = 1.0 - hamming_distance;
                    let numerator = 2.0 * jaccard;
                    let denominator = 1.0 + jaccard;
                    let fraction = (numerator as f64) / (denominator as f64);
                    let distance = -fraction.ln() / (kmer_size as f64);
                    (i, j, distance)
                })
        })
        .collect();

    let mut matrix = vec![vec![0.0_f64; n]; n];
    for &(i, j, dist) in &distances {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    let mut phylip_data = Vec::new();
    writeln!(phylip_data, "{}", n).unwrap();
    for i in 0..n {
        let name = Path::new(&genomes[i])
            .file_name()
            .and_then(|os_str| os_str.to_str())
            .unwrap_or(&genomes[i])
            .to_string();
        write!(phylip_data, "{:10}", name).unwrap();
        for j in 0..n {
            write!(phylip_data, " {:8.6}", matrix[i][j]).unwrap();
        }
        writeln!(phylip_data).unwrap();
    }

    phylip_data
}

fn build_tree(
    tree_algo: &TreeAlgo,
    chunk_size: usize,
    naive_percentage: usize,
    phylip_data: &[u8]
) -> String {
    let distance_matrix =
        DistanceMatrix::read_from_phylip(&phylip_data[..]).expect("Error reading phylip matrix");

    let graph = match tree_algo {
        TreeAlgo::Naive => {
            NeighborJoiningSolver::<Canonical>::default(distance_matrix).solve()
        }
        TreeAlgo::RapidNJ => {
            NeighborJoiningSolver::<RapidBtrees>::build(distance_matrix, chunk_size).solve()
        }
        TreeAlgo::Hybrid => {
            let naive_steps = distance_matrix.size() * naive_percentage / 100;
            NeighborJoiningSolver::<Hybrid>::build(distance_matrix, chunk_size, naive_steps).solve()
        }
    }
    .expect("Error constructing tree");

    speedytree::to_newick(&graph)
}

fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();
    let matches = Command::new("BinDashtree")
        .version("0.1.0")
        .about("Binwise Densified MinHash and Rapid Neighbor-joining Tree Construction")
        .arg(
            Arg::new("input_list")
                .short('i')
                .long("input")
                .value_name("INPUT_LIST_FILE")
                .help("Genome list file (one FASTA/FNA file per line), .gz supported")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("kmer_size")
                .short('k')
                .long("kmer_size")
                .value_name("KMER_SIZE")
                .help("K-mer size")
                .default_value("16")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("sketch_size")
                .short('s')
                .long("sketch_size")
                .value_name("SKETCH_SIZE")
                .help("MinHash sketch size")
                .default_value("10240")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("dens_opt")
                .short('d')
                .long("densification")
                .value_name("DENS_OPT")
                .help("Densification strategy: 0=Optimal Densification, 1=Reverse Optimal Densification/faster Densification")
                .default_value("0")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads to use in parallel")
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("tree_method")
                .long("tree")
                .value_name("TREE_METHOD")
                .help("Tree construction method: naive, rapidnj, hybrid")
                .default_value("rapidnj")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("chunk_size")
                .long("chunk_size")
                .help("Chunk size for RapidNJ/Hybrid methods")
                .default_value("30")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("naive_percentage")
                .long("naive_percentage")
                .help("Percentage of steps naive for hybrid method")
                .default_value("90")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output_matrix")
                .long("output_matrix")
                .value_name("OUTPUT_MATRIX_FILE")
                .help("Output the phylip distance matrix to a file")
                .required(false)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output_tree")
                .long("output_tree")
                .value_name("OUTPUT_TREE_FILE")
                .help("Output the resulting tree in Newick format to a file")
                .required(true)
                .action(ArgAction::Set),
        )
        .get_matches();

    let input_list = matches.get_one::<String>("input_list").unwrap().to_string();
    let kmer_size = *matches.get_one::<usize>("kmer_size").unwrap();
    let sketch_size = *matches.get_one::<usize>("sketch_size").unwrap();
    let dens = *matches.get_one::<usize>("dens_opt").unwrap();
    let threads = *matches.get_one::<usize>("threads").unwrap();
    let tree_method = matches.get_one::<String>("tree_method").unwrap();
    let chunk_size = *matches.get_one::<usize>("chunk_size").unwrap();
    let naive_percentage = *matches.get_one::<usize>("naive_percentage").unwrap();
    let output_matrix = matches.get_one::<String>("output_matrix").cloned();
    let output_tree = matches.get_one::<String>("output_tree").cloned();

    let tree_algo: TreeAlgo = tree_method.parse().expect("Invalid tree method");

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let file = File::open(&input_list).expect("Cannot open input genome list file");
    let reader = BufReader::new(file);
    let genomes: Vec<String> = reader
        .lines()
        .map(|line| line.expect("Error reading genome list"))
        .collect();

    let sketch_args = SeqSketcherParams::new(kmer_size, sketch_size, SketchAlgo::OPTDENS, DataType::DNA);

    println!("Sketching all genomes...");
    let sketches = sketch_genomes(kmer_size, dens, &sketch_args, &genomes);

    println!("Building PHYLIP distance matrix...");
    let phylip_data = build_distance_matrix(&sketches, kmer_size, &genomes);

    if let Some(filename) = output_matrix.as_ref() {
        let mut f = BufWriter::new(File::create(filename).expect("Cannot create matrix file"));
        f.write_all(&phylip_data).expect("Error writing matrix");
    }

    println!("Constructing the tree...");
    let newick = build_tree(&tree_algo, chunk_size, naive_percentage, &phylip_data);

    if let Some(filename) = output_tree {
        let mut f = BufWriter::new(File::create(filename).expect("Cannot create tree file"));
        writeln!(f, "{}", newick).expect("Error writing tree");
    } else {
        println!("{}", newick);
    }
}
