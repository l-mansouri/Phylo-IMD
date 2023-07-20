include { ALIGNMENT } from '../subworkflows/alignment.nf'
include { SELECTING_RANDOM_POSITIONS } from '../subworkflows/sampling_for_titration.nf'
include { GENERATE_SEQ_ME_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TITR_TREES } from '../subworkflows/IMD_ME_TREES.nf'


workflow titration_Phylo_IMD{
    take:
        input_fasta
        templates
        structures

    main:
        ALIGNMENT(input_fasta, templates, structures)
        SELECTING_RANDOM_POSITIONS(ALIGNMENT.out.fasta_aln)
        GENERATE_SEQ_ME_TREES(SELECTING_RANDOM_POSITIONS.out.randomized_msa)
        GENERATE_SEQ_ML_TREES(SELECTING_RANDOM_POSITIONS.out.randomized_msa)
        GENERATE_IMD_ME_TITR_TREES(ALIGNMENT.out.fasta_aln, SELECTING_RANDOM_POSITIONS.out.randomized_col, templates, structures)       

}