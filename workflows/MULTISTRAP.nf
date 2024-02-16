include { MTM_ALIGNMENT                         } from '../subworkflows/MTM_alignment.nf'
include { GENERATE_SEQ_ME_TREES                 } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES                 } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TREES                 } from '../subworkflows/IMD_ME_TREES.nf'
include { COMBINE_BOOTSTRAP_SUPPORTS            } from '../modules/combine_bootstrap_supports.nf'
include { MERGE_REPLICATES                      } from '../modules/merge_replicates.nf'

workflow MULTISTRAP{
    take:
        input_fasta
        templates
        structures

    main:

        // Compute alignment 
        MTM_ALIGNMENT(input_fasta, templates, structures)

        // Compute IMD trees and bootstrap
        GENERATE_IMD_ME_TREES(MTM_ALIGNMENT.out.fasta_aln, templates, structures)
        MERGE_REPLICATES(GENERATE_IMD_ME_TREES.out.tree.groupTuple())
        MERGE_REPLICATES.out.replicates.map{
            id, tree, replicates  -> [id.replaceAll("_.*", ""), tree, replicates]
        }.set{IMD_TREE_AND_REPLICATES}

        // Compute sequence_based trees and bootsrap
        if (params.seq_tree=='ME'){
            GENERATE_SEQ_ME_TREES(MTM_ALIGNMENT.out.fasta_aln)
            SEQ_TREE_AND_REPLICATES = GENERATE_SEQ_ME_TREES.out.tree.combine(GENERATE_SEQ_ME_TREES.out.replicates, by:0)
        }
        else if (params.align=='ML'){
            GENERATE_SEQ_ML_TREES(MTM_ALIGNMENT.out.fasta_aln)
            SEQ_TREE_AND_REPLICATES = GENERATE_SEQ_ML_TREES.out.tree.combine(GENERATE_SEQ_ML_TREES.out.replicates, by:0)
        }

        // Combine bootstrap supports
        COMBINE_BOOTSTRAP_SUPPORTS(SEQ_TREE_AND_REPLICATES.combine(IMD_TREE_AND_REPLICATES, by:0))

}