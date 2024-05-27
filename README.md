# MULTISTRAP
## Boosting phylogenetic boostrap with structural information
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7437267-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7437267)

Multistrap is a toolkit designed to calculate and combine phylogenetic bootstrap support values. It generates these support values using both sequence and structural data, and then combines them. 

For more details, see the associated manuscript: ["Boosting phylogenetic bootstrap with structural information"](https://zenodo.org/records/11187505).


## Installation 

### Requirements
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) version >= 23.10.
  
  Make sure you have the [right Java version installed](https://stackoverflow.com/questions/74103638/how-do-i-install-the-correct-version-of-java-for-nextflow) 
```
curl -s https://get.nextflow.io | bash
chmod +x nextflow
# /usr/local/bin or any other executable path
sudo mv nextflow /usr/local/bin
```
- Either [Docker](https://docs.docker.com/engine/install/) version >= 20.10 or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#) version >= 3.7.

**! Remember to start Docker** before starting the pipeline. 

Multistrap was tested on Scientific Linux release 7.2.

### Get Multistrap
Multistrap is distributed as a Nextflow pipeline. 
To obtain the source code: 
```
curl -L -o main.zip https://github.com/l-mansouri/Phylo-IMD/archive/refs/heads/main.zip 
unzip main.zip
cd Phylo-IMD-main
```
On a normal Desktop computer this step should take seconds.
Now you are ready to run Multistrap!

## Run Multistrap

Multistrap per default will: 
- compute the mTMalign MSA
- compute the sequence based tree and corresponding bootstrap replicates (ME or ML tree)
- compute the IMD tree and corresponding bootstrap replicates
- return the tree with the combined (multistrap) bootstrap support values


### On a test dataset

```
nextflow run main.nf -profile multistrap,test,docker
```

If you want to use singularity: 

```
nextflow run main.nf -profile multistrap,test,singularity
```
<details markdown="1">
<summary>More</summary>

This will use the test [data](https://github.com/l-mansouri/Phylo-IMD/tree/main/data) to run multistrap. 
We use `--seq_tree ME` as ML takes longer and this is meant to be just a basic test.
`replicatesNum` is also set to 10, to speed up the run.
In a normal Desktop computer this should take few minutes to complete. 
</details>


### On your dataset

To obtain the combined bootstrap support values in your own dataset please use the multistrap profile as shown in the following lines. 
To see how to properly prepare the input files, look into the example dataset in the [data](https://github.com/l-mansouri/Phylo-IMD/tree/main/data). 

The command line: 

```
nextflow run main.nf -profile multistrap -fasta <id.fasta> -templates <id.template> -pdbs mypdbs/* -seq_tree <ML|ME>
```

- `fasta` is a fasta file with the sequences you want to build the tree on. 
- `pdbs` is all the pdbs associated to the sequences present in your fasta file. 
- `templates` is a file with the explicit mapping of each sequence in your fasta file and each pdb you are providing.
  The template files should follow the corresponding syntax (mTM-align or 3D-Coffee correspondingly). You can find examples for both in the data folder.

<details markdown="1">
<summary>Output files</summary>

- `results/`
  - `msas/*.fa`: **alignment** files. 
  - `trees/`:  **trees** computed using your preferred sequence method (ME or ML) (trees/<ME|ML> folder) and the IMD trees (trees/IMD folder). **Tree replicates** are found in the replicates folder within the ME|ML|IMD folders respectively.
  - `multistrap_<ME|ML>_and_IMD/` the **Bootstrap support values** are stored as node labels in the trees found in multistrap_<ME|ML>_and_IMD folder. Here you will find one file with the tree with the <ME|ML> support values and one with the <IMD> bootstrap support values separatly and one with the multistrap support values.
  </details>




## Pipelines parameters 

You can modify the default pipeline parameters by using: 


<details markdown="1">
<summary>Parameters</summary>

  - Input parameters
      - `fasta` is a fasta file with the sequences you want to build the tree on. 
      - `pdbs` is all the pdbs associated to the sequences present in your fasta file. 
      - `templates` is a file with the explicit mapping of each sequence in your fasta file and each pdb you are providing.
        The template files should follow the corresponding syntax (mTM-align or 3D-Coffee correspondingly). You can find examples for both in the data folder.
  - Parameters for tree computation:
      - `seq_tree` determines the type of sequence based tree to be computed: either ME or ML. Default: ML. 
      - `gammaRate` that determines the gamma rate for FastME tree reconstruction. Default: 1.0.
      - `seedValue` that is the random seed for FastME tree reconstruction. Default: 5. 
      - `replicatesNum` that determines the number of bootstrap replicates. Default: 100. 
      - `tree_mode` that determines the distance mode to run the IMD distance matrix computation. Default: 10.
  - Output parameter:
      - `output` that determines where to store the outputs that the pipeline publishes. Default: ./results.
</details>



## Overview of the repository
For a more detailed overiview of the content of the repository please refer to [overview](https://github.com/l-mansouri/Phylo-IMD/blob/main/Overview.md)


## Analysis

In the paper we perform an extensive benchmark and produce accessory analyses to assess the robustness and validity of Multistrap. 
For more information on how to reproduce this please refer to [analysis](https://github.com/l-mansouri/Phylo-IMD/blob/main/Analysis.md)
