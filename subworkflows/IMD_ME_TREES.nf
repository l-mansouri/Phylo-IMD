include { computing_IMD_matrices } from '../modules/computing_IMD_matrices.nf'
include { computing_IMD_tr_matrices } from '../modules/computing_trimmed_imd_mats.nf'
include { computing_IMD_titr_matrices } from '../modules/computing_titration_imd_mats.nf'
include { computing_IMD_titr_matrices_bs } from '../modules/computing_titration_imd_bs_mats.nf'
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
        extracting_matrices.out.flatten().set{for_IMD_trees}
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
        extracting_matrices.out.flatten().set{for_IMD_trees}
        computing_IMD-ME_trees(for_IMD_trees)

    emit:
    computing_IMD-ME_trees.out
}

workflow GENERATE_IMD_ME_TITR_TREES{
    take:
        fasta_aln
        sampled_pos
        templates
        structures
    main:

        fasta_aln.combine(templates, by:0).combine(sampled_pos, by:0).combine(structures, by:0).set{for_IMD_mat_titr}

        computing_IMD_titr_matrices(for_IMD_mat_titr)
        computing_IMD-ME_trees(computing_IMD_titr_matrices.out)

    emit:
    computing_IMD-ME_trees.out
}

workflow GENERATE_IMD_ME_TITR_BS_TREES{
    take:
        fasta_aln
        sampled_pos
        templates
        structures
    main:

        fasta_aln.combine(templates, by:0).combine(sampled_pos, by:0).combine(structures, by:0).set{for_IMD_mat_titr}

        computing_IMD_titr_matrices_bs(for_IMD_mat_titr)
        extracting_matrices(computing_IMD_tr_matrices.out)
        extracting_matrices.out.flatten().set{for_IMD_trees}
        computing_IMD-ME_trees(for_IMD_trees.out)

    emit:
    computing_IMD-ME_trees.out
}