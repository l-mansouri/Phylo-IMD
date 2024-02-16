# Phylo-IMD
This repo contains a nextflow pipeline that reproduces the main results of the Phylo-IMD paper. It allows computing phylogenetic trees based on the **IMD** distance and **TMscore-based trees**, ME trees (using **FastME**), and ML trees (using **iQtree**). It implements 3 different alignment methods (**mTMalign**, **3d-coffee** and **t-coffee**). Other options include the trimming of the alignment (with **trimal**).

The pipeline is fully containerised and the docker files used to generate the containers used by the pipeline are present in `Dockerfiles`.

The different processes are implemented as `modules` and the pipeline is composed by `subworkflows` and `workflows` called by `main.nf`.

An example of input files can be found in the `data` directory.

## Pipeline profiles and parameters

The pipeline includes different profiles that set default parameters values: 
- `untrimmed_mTM` that runs the standard analysis with the **mTMalign untrimmed** alignment
- `trimmed_mTM` that runs the standard analysis with the **mTMalign** alignment trimmed by using **trimal**
- `untrimmed_3dc` that runs the standard analysis with the **3d-coffee untrimmed** alignment
- `trimmed_3dc` that runs the standard analysis with the **3d-coffee** alignment trimmed by using **trimal**
- `untrimmed_tc` that runs the standard analysis with the **t-coffee untrimmed** alignment
- `trimmed_tc` that runs the standard analysis with the **t-coffee** alignment trimmed by using **trimal**
-  `titration` that runs the titration analysis present in the paper; by default, it computes and uses **untrimmed mTMalign** alignment; if you wish to change the alignment method, you would have to change the `align` parameter and the input paths `fasta` and `template`
-  `titration_bootstrap` that runs the bootstrap for the titration-based trees; by default, it computes and uses **untrimmed mTMalign** alignment; if you wish to change the alignment method, you would have to change the `align` and the input paths `fasta` and `template`

But in general, the parameters are:
- General parameters:
    - `mode` that determines the type of analysis. It can either be `standard`, `titration`, or `titration_bootstrap`
    - `align` that determines the alignment method to use. It can be `mTMalign`, `3dcoffee`, `tcoffee`
    - `trimmer` that determined the type of trimming of the alignment. It can be `untrimmed` or `trimmal`
- Inputs parameters:
    - `fasta` that is the path to the input sequences (or the list of pdbs, in the case of mTMalign)
    - `templates` that is the path to the template files needed by t-coffee to compute the IMD distance matrices (as well as the alignment, in the case of 3d-coffee)
    - `pdb` that is the path to the pdb files
- Parameters for tree computation:
    - `gammaRate` that determines the gamma rate for FastME tree reconstruction
    - `seedValue` that is the random seed for FastME tree reconstruction
    - `replicatesNum` that determines the number of bootstrap replicates
    -  `tree_mode` that determines the distance mode to run the IMD distance matrix computation, default is 10.
- Output parameter:
    - `output` that determines where to store the outputs that the pipeline publishes

## How to run the pipeline

To run the main analysis present in the paper, you should run it like:
`nextflow run main.nf -profile untrimmed_mTM`

To run the titration, you should run the pipeline like:
`nextflow run main.nf -profile titration`

To run bootstrap on titration trees:
`nextflow run main.nf -profile titration_bootstrap`

## MULTISTRAP

To obtain the combined bootstrap support values in any dataset, the multistrap profile should be used. 
To see how to properly prepare the input files, please look into the example dataset in the ./data folder. 

The command line should be: 

`nextflow run main.nf -profile multistrap -fasta <id.fasta> -templates <id.template> -pdbs <id.seq1.pdb, id.seq2.pdb .. >`

This will: 
- compute the mTMalign MSA
- compute the sequence based tree and bootstrap replicates (Seq-ME or Seq-ML)
- compute the IMD-ME-based tree and bootstrap replicates
- return the tree with the combined (multistrap) bootstrap support values


