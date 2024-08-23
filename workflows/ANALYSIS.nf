include { evaluate_nirmsd                    } from '../modules/evaluate_nirmsd.nf'
include { split_analysis                     } from '../subworkflows/split_analysis.nf'
include { compute_splits                     } from '../subworkflows/compute_splits.nf'
include { TCOFFEE_SEQREFORMAT                } from '../modules/TCOFFEE_SEQREFORMAT.nf'


workflow ANALYSIS{

    take:
        msas
        templates 
        pdb
        trees_ch
        replicates_25_ch
        replicates_200_ch


    main: 

        msas
        .combine(templates, by:0)
        .combine(pdb, by:0)
        .set{aln_ch}


        // --------------------------------
        //     COMPUTE PERC SIMILARITY 
        // --------------------------------
        TCOFFEE_SEQREFORMAT(msas)


        // ------------------
        //     iRMSD
        // ------------------

        evaluate_nirmsd( aln_ch )


        // ------------------
        //     SPLITS- AUC
        // ------------------
    
        // mix replicates
        replicates = replicates_25_ch.mix(replicates_200_ch)

        // combine trees and replicates
        trees_ch.combine(replicates, by:0)
                .set{replicate_trees}

        // // prepare the type 
        replicate_trees.map{ fam, splits, tree, bs, replicates -> [fam, tree, replicates, splits+"_"+bs] }.set{ch_for_splits}

        // Compute the splits
        compute_splits(ch_for_splits)






}