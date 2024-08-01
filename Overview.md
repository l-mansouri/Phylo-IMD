# Overview of the repository

The different processes are implemented as [modules](https://github.com/l-mansouri/Phylo-IMD/tree/main/modules) and the pipeline is composed by [subworkflows](https://github.com/l-mansouri/Phylo-IMD/tree/main/subworkflows) and [workflows](https://github.com/l-mansouri/Phylo-IMD/tree/main/workflows). All the workflows and subworkflows are orchestrated by the main script [main.nf](https://github.com/l-mansouri/Phylo-IMD/blob/main/main.nf).

For the sake of reproducibility and seamless deployment, each of the processes has an associated container. The associated dockerfiles available at [Dockerfiles](https://github.com/l-mansouri/Phylo-IMD/tree/main/Dockerfiles).

An example of input files can be found in [data](https://github.com/l-mansouri/Phylo-IMD/tree/main/data).

Available in the pipeline: 
- Computation of multiple sequence alignment.
  Available MSAs: [mTM-align](https://yanglab.qd.sdu.edu.cn/mTM-align/), [3D-Coffee](https://tcoffee.org/Projects/expresso/index.html) and [T-Coffee](https://github.com/cbcrg/tcoffee). 
- (optional) Trimming of the multiple sequence alignment: 
  Available: trimming with [trimAl](https://vicfero.github.io/trimal/). 
- Computation of phylogenetic treees. 
  Sequence based: ME trees (using [FastME](http://www.atgc-montpellier.fr/fastme/)), and ML trees (using [iQtree](http://www.iqtree.org/))
  Structure based: IMD trees (using the IMD distance + [FastME](http://www.atgc-montpellier.fr/fastme/) as described in the manuscript).
- Computation of bootstrap support values using: 
  - Sequence data 
  - Structure data 
  - The combination of sequence and structure data
