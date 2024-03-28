include { evaluate_nirmsd                    } from '../modules/evaluate_nirmsd.nf'
include { split_analysis                     } from '../subworkflows/split_analysis.nf'
include { compute_splits                     } from '../subworkflows/compute_splits.nf'
include { comparing_splitfiles               } from '../modules/comparing_splitfiles.nf'
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

        //evaluate_nirmsd( aln_ch )


        // ------------------
        //     SPLITS- AUC
        // ------------------

        
        if ( params.trees ) {
            Channel
                .fromPath(params.trees)
                .map { item -> [ item.baseName.split('_')[0], item.getParent().getBaseName().split('_')[0].strip(), item] }
                .map { fam, type, tree -> [fam, type == "1d" ? "MEsplits" : type, tree] }
                .map { fam, type, tree -> [fam, type == "3d" ? "IMDsplits" : type, tree] }
                .map { fam, type, tree -> [fam, type == "ML" ? "MLsplits" : type, tree] }
                .set{trees_ch}
        }
        
        if ( params.replicates_25 ) {
            Channel
                .fromPath(params.replicates_25)
                .map { item -> [ item.baseName.split('_')[0],item.getParent().getBaseName().split('_')[0].strip(), item] }
                .map { fam, type, tree -> [fam, type == "1d" ? "MEbs25" : type, tree] }
                .map { fam, type, tree -> [fam, type == "3d" ? "IMDbs25" : type, tree] }
                .map { fam, type, tree -> [fam, type == "ML" ? "MLbs25" : type, tree] }
                .set{replicates_25_ch}
        }

        if ( params.replicates_200 ) {
            Channel
                .fromPath(params.replicates_200)
                .map { item -> [ item.baseName.split('_')[0],item.getParent().getBaseName().split('_')[0].strip(), item] }
                .map { fam, type, tree -> [fam, type == "1d" ? "MEbs200" : type, tree] }
                .map { fam, type, tree -> [fam, type == "3d" ? "IMDbs200" : type, tree] }
                .map { fam, type, tree -> [fam, type == "ML" ? "MLbs200" : type, tree] }
                .set{replicates_200_ch}
        }

    
        // mix replicates
        replicates = replicates_25_ch.mix(replicates_200_ch)

        // combine trees and replicates
        trees_ch.combine(replicates, by:0)
                .set{replicate_trees}

        // // prepare the type 
        replicate_trees.map{ fam, splits, tree, bs, replicates -> [fam, tree, replicates, splits+"_"+bs] }.set{ch_for_splits}

        // Compute the splits
        compute_splits(ch_for_splits)

        compute_splits.out.splits.view()

        // IMD splits, ME splits, ML splits, IMB bs, ME bs, ML bs
        


}