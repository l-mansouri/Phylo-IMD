include { evaluate_nirmsd                    } from '../modules/evaluate_nirmsd.nf'
include { split_analysis                     } from '../subworkflows/split_analysis.nf'

workflow ANALYSIS{

    take:
        msas
        templates 
        pdb

    main: 

        msas
        .combine(templates, by:0)
        .combine(pdb, by:0)
        .set{aln_ch}


        // ------------------
        //     iRMSD
        // ------------------

        evaluate_nirmsd( aln_ch )


        // ------------------
        //     SPLITS- AUC
        // ------------------



}