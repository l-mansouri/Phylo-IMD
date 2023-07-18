include { computing_IMD_matrices } from '../modules/computing_IMD_matrices.nf'
include {computing_IMD_tr_matrices} from '../modules/computing_trimmed_imd_mats.nf'
include { extracting_matrices } from '../modules/extracting_IMD_matrices.nf'
include { computing_IMD-ME_trees } from '../modules/computing_IMD_trees.nf'

workflow GENERATE_IMD_ME_TREES{
    take:
        fasta_aln
        templates
        structures
    main:

        fasta_aln.combine(templates, by:0).combine(structures, by:0).set{for_IMD_mat}

        computing_IMD_matrices(for_IMD_mat)
        extracting_matrices(computing_IMD_matrices.out)
        extracting_matrices.flatten().set{for_IMD_trees}
        computing_IMD-ME_trees(for_IMD_trees)

    emit:
    computing_IMD-ME_trees.out
}

workflow GENERATE_IMD_ME_TR_TREES{
    take:
        fasta_aln
        mapped_pos
        templates
        structures
    main:

        fasta_aln.combine(templates, by:0).combine(mapped_pos, by:0).combine(structures, by:0).set{for_IMD_mat_tr}

        computing_IMD_tr_matrices(for_IMD_mat_tr)
        extracting_matrices(computing_IMD_tr_matrices.out)
        extracting_matrices.flatten().set{for_IMD_trees}
        computing_IMD-ME_trees(for_IMD_trees)

    emit:
    computing_IMD-ME_trees.out
}