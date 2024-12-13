[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bindashtree/README.html)
![](https://anaconda.org/bioconda/bindashtree/badges/license.svg)
![](https://anaconda.org/bioconda/bindashtree/badges/version.svg)
![](https://anaconda.org/bioconda/bindashtree/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/bindashtree/badges/platforms.svg)



## Build genome phylogenetic tree via Binwise Densified MinHash and Rapid Neighbor-joining
One Permutation Hashing with Optimal Densification (1,2) can be use for genomic distance estimation (1-ANI) and then we can perform rapid neighbor-joining (4) based on the genomic distance. We also provided a new densification strategy called faster densification (or reverse optimal densification) (3), which is more accurate and faster for large sketch size (5).

## Install

### Install via bioconda

```bash
conda install -c bioconda -c conda-forge bindashtree

```
### Pre-built binary
```bash
### Pre-built library (Linux only)
wget https://github.com/jianshu93/bindashtree/releases/download/v0.1.0/bindashtree_Linux_86-64_v0.1.0.zip
unzip bindashtree_Linux_86-64_v0.1.0.zip
chmod a+x ./bindashtree
./bindashtree -h
```

### Install from source
```bash
### compling from source
git clone https://github.com/jianshu93/bindashtree.git
cd bindashtree
cargo build --release
./target/release/bindashtree -h
```

## Usage
```bash

 ************** initializing logger *****************

Binwise Densified MinHash and Rapid Neighbor-joining Tree Construction

Usage: bindashtree [OPTIONS] --input <INPUT_LIST_FILE> --output_tree <OUTPUT_TREE_FILE>

Options:
  -i, --input <INPUT_LIST_FILE>
          Genome list file (one FASTA/FNA file per line), gz supported
  -k, --kmer_size <KMER_SIZE>
          K-mer size [default: 16]
  -s, --sketch_size <SKETCH_SIZE>
          MinHash sketch size [default: 10240]
  -d, --densification <DENS_OPT>
          Densification strategy: 0=Optimal Densification, 1=Reverse Optimal Densification/faster Densification [default: 0]
  -t, --threads <THREADS>
          Number of threads to use in parallel [default: 1]
      --tree <TREE_METHOD>
          Tree construction method: naive, rapidnj, hybrid [default: rapidnj]
      --chunk_size <chunk_size>
          Chunk size for RapidNJ/Hybrid methods [default: 30]
      --naive_percentage <naive_percentage>
          Percentage of steps naive for hybrid method [default: 90]
      --output_matrix <OUTPUT_MATRIX_FILE>
          Output the phylip distance matrix to a file
      --output_tree <OUTPUT_TREE_FILE>
          Output the resulting tree in Newick format to a file
  -h, --help
          Print help
  -V, --version
          Print version

```
## Output explanation

A newick format tree and phylip format distance matrix will be the output depending on your options. Tree can be visualized via [Figtree](http://tree.bio.ed.ac.uk/software/figtree/), [iTOL](https://itol.embl.de) or [ggtree](https://guangchuangyu.github.io/software/ggtree/)


## Testing dataset

```bash
ls ./data/*.fna.gz > name.txt
./target/release/bindashtree -i name.txt -k 16 -s 10240 -d 1 -t 8 --output_tree try.nwk
```

## References

1.Li, P., Owen, A. and Zhang, C.H., 2012. One permutation hashing. Advances in Neural Information Processing Systems, 25.

2.Shrivastava, A., 2017, July. Optimal densification for fast and accurate minwise hashing. In International Conference on Machine Learning (pp. 3154-3163). PMLR.

3.Mai, T., Rao, A., Kapilevich, M., Rossi, R., Abbasi-Yadkori, Y. and Sinha, R., 2020, August. On densification for minwise hashing. In Uncertainty in Artificial Intelligence (pp. 831-840). PMLR.

4.Simonsen, M., Mailund, T. and Pedersen, C.N., 2008. Rapid neighbour-joining. In Algorithms in Bioinformatics: 8th International Workshop, WABI 2008, Karlsruhe, Germany, September 15-19, 2008. Proceedings 8 (pp. 113-122). Springer Berlin Heidelberg.

5.Zhao, J., Zhao, X., Pierre-Both, J. and Konstantinidis, K.T., 2024. BinDash 2.0: new MinHash scheme allows ultra-fast and accurate genome search and comparisons. bioRxiv, pp.2024-03.