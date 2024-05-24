# Analysis
### How to reproduce the results reported in the paper

In the paper we perform an extensive benchmark and produce accessory analyses to assess the robustness and validity of Multistrap. 

You can get the input data for the full dataset at: [https://zenodo.org/records/7447443](https://zenodo.org/records/7447443)

We hereby provide the command lines necessary to reproduce the reported results. To run it on the real dataset, please download it and update the path of fasta, templates and pdb with the path you download the data to. 

To run the main analysis present in the paper (uses untrimmed mTM-align MSAs to produce ME, ML and IMD trees):
```nextflow run main.nf -profile untrimmed_mTM```

To modify the aligner or trimming to use you can do it with the following parameters: 
  - `align` that determines the alignment method to use. It can be `mTMalign`, `sap_tmalign`, `tcoffee`
  - `trimmer` that determined the type of trimming of the alignment. It can be `untrimmed` or `trimmal`

To run the titration, you should run the pipeline like:

```nextflow run main.nf -profile titration```

To run bootstrap on titration trees:

```nextflow run main.nf -profile titration_bootstrap```

To compute extra statistics needed for the complete analysis shown in the paper, such as the calculation of NiRMSD, you can use the profile analysis.

```nextflow run main.nf -profile analysis -msas <id.aln> -templates <id.template> -pdbs my_pdbs/*```

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

