include { evaluate_nirmsd                    } from '../modules/evaluate_nirmsd.nf'
include { split_analysis                     } from '../subworkflows/split_analysis.nf'
include { collecting_replicates              } from '../modules/collecting_replicates.nf'
include { computing_splitfiles               } from '../modules/computing_splitfiles.nf'

workflow ANALYSIS{

    take:
        msas
        templates 
        pdb
        trees
        replicates

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


        trees_ch.combine(replicates_ch, by:0).set{replicate_trees}.view()

        collecting_replicates(replicate_trees)
        computing_splitfiles(collecting_replicates.out)       



}