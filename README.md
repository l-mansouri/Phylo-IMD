# MULTISTRAP
## Boosting phylogenetic boostrap with structural information
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7437267)

Multistrap is a toolkit designed to calculate and combine phylogenetic bootstrap support values. It generates these support values using both sequence and structural data, and then combines them. 

For more details, see the associated manuscript: ["Boosting phylogenetic bootstrap with structural information"](https://zenodo.org/records/11187505).


## Installation 

### Requirements
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) version >= 23.10.
```
curl -s https://get.nextflow.io | bash
chmod +x nextflow
# /usr/local/bin or any other executable path
sudo mv nextflow /usr/local/bin
```
- Either [Docker](https://docs.docker.com/engine/install/) version >= 20.10 or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#) version >= 3.7.

Multistrap was tested on Scientific Linux release 7.2.

### Get Multistrap
Multistrap is distributed as a Nextflow pipeline. 
To obtain the source code: 
```
wget https://github.com/l-mansouri/Phylo-IMD/archive/refs/heads/main.zip
unzip main.zip
cd Phylo-IMD-main
```
On a normal Desktop computer this step should take seconds.
Now you are ready to run Multistrap!

## Run Multistrap

### DEMO: deploy multistrap on a test dataset
To deploy multistrap on the provided test dataset using docker: 

`nextflow run main.nf -profile multistrap,test,docker --seq_tree ME`

To deploy multistrap on the provided test dataset using singularity: 
`nextflow run main.nf -profile multistrap,test,singularity --seq_tree ME`


This will use the test [data](https://github.com/l-mansouri/Phylo-IMD/tree/main/data) to run multistrap. 
We use --seq_tree ME as ML takes longer and this is meant to be just a basic test. 
In a normal Desktop computer this should take ~10 minutes to complete. 

#### Output
This will produce a folder results: 
 - the **MSAs** computed by mTM-align both in FASTA format (msas folder)
 - the **Trees**  computed using your preferred sequence method (ME or ML) (trees/<ME|ML> folder) and the IMD trees (trees/IMD folder). **Tree replicates** are found in the replicates folder within the ME|ML|IMD folders respectively.
 - the **Bootstrap support values** are stored as node labels in the trees found in multistrap_<ME|ML>_and_IMD folder. Here you will find one file with the tree with the <ME|ML> support values and one with the <IMD> bootstrap support values separatly and one with the multistrap support values.


### Run MULTISTRAP on your dataset

To obtain the combined bootstrap support values in your own dataset please use the multistrap profile as shown in the following lines. 
To see how to properly prepare the input files, look into the example dataset in the [data](https://github.com/l-mansouri/Phylo-IMD/tree/main/data). 

The command line: 

`nextflow run main.nf -profile multistrap -fasta <id.fasta> -templates <id.template> -pdbs mypdbs/* -seq_tree <ML|ME>`


## Pipelines parameters 
- Input parameters
    - `fasta` is a fasta file with the sequences you want to build the tree on. 
    - `pdbs` is all the pdbs associated to the sequences present in your fasta file. 
    - `templates` is a file with the explicit mapping of each sequence in your fasta file and each pdb you are providing.
      The template files should follow the corresponding syntax (mTM-align or 3D-Coffee correspondingly). You can find examples for both in the data folder.
- Parameters for tree computation:
    - `seq_tree` determines the type of sequence based tree to be computed: either ME or ML. Default: ML. 
    - `gammaRate` that determines the gamma rate for FastME tree reconstruction
    - `seedValue` that is the random seed for FastME tree reconstruction
    - `replicatesNum` that determines the number of bootstrap replicates
    - `tree_mode` that determines the distance mode to run the IMD distance matrix computation, default is 10.
- Output parameter:
    - `output` that determines where to store the outputs that the pipeline publishes


Mulistrap per default will: 
- compute the mTMalign MSA
- compute the sequence based tree and corresponding bootstrap replicates (ME or ML tree)
- compute the IMD tree and corresponding bootstrap replicates
- return the tree with the combined (multistrap) bootstrap support values



## Overview of the repository

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

## Analysis
### How to reproduce the results reported in the paper

In the paper we perform an extensive benchmark and produce accessory analyses to assess the robustness and validity of Multistrap. 

You can get the input data for the full dataset at: [https://zenodo.org/records/7447443](https://zenodo.org/records/7447443)

We hereby provide the command lines necessary to reproduce the reported results:

To run the main analysis present in the paper (uses untrimmed mTM-align MSAs to produce ME, ML and IMD trees):
`nextflow run main.nf -profile untrimmed_mTM`

To modify the aligner or trimming to use you can do it with the following parameters: 
  - `align` that determines the alignment method to use. It can be `mTMalign`, `sap_tmalign`, `tcoffee`
  - `trimmer` that determined the type of trimming of the alignment. It can be `untrimmed` or `trimmal`

To run the titration, you should run the pipeline like:
`nextflow run main.nf -profile titration`

To run bootstrap on titration trees:
`nextflow run main.nf -profile titration_bootstrap`

To compute extra statistics needed for the complete analysis shown in the paper, such as the calculation of NiRMSD, you can use the profile analysis.

`nextflow run main.nf -profile analysis -msas <id.aln> -templates <id.template> -pdbs my_pdbs/*`

The scripts to reproduce the figures and tables in the manuscript are available in the [analysis](https://github.com/l-mansouri/Phylo-IMD/tree/main/analysis) folder. 

The pipeline includes additional profiles that set default parameters values: 
- `untrimmed_mTM` that runs the standard analysis with the **mTMalign untrimmed** alignment
- `trimmed_mTM` that runs the standard analysis with the **mTMalign** alignment trimmed by using **trimal**
- `untrimmed_3dc` that runs the standard analysis with the **3d-coffee untrimmed** alignment
- `trimmed_3dc` that runs the standard analysis with the **3d-coffee** alignment trimmed by using **trimal**
- `untrimmed_tc` that runs the standard analysis with the **t-coffee untrimmed** alignment
- `trimmed_tc` that runs the standard analysis with the **t-coffee** alignment trimmed by using **trimal**
- `titration` that runs the titration analysis present in the paper; by default, it computes and uses **untrimmed mTMalign** alignment; if you wish to change the alignment method, you would have to change the `align` parameter and the input paths `fasta` and `template`
- `titration_bootstrap` that runs the bootstrap for the titration-based trees; by default, it computes and uses **untrimmed mTMalign** alignment; if you wish to change the alignment method, you would have to change the `align` and the input paths `fasta` and `template`

