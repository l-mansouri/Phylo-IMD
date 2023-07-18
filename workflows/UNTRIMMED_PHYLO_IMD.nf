include { ALIGNMENT } from '../subworkflows/alignment.nf'
include { GENERATE_SEQ_ME_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TREES } from '../subworkflows/IMD_ME_TREES.nf'
include { GENERATE_TM_ME_TREES } from '../subworkflows/TM_ME_TREES.nf'


workflow untrimmed_Phylo_IMD{
    take:
        input_fasta
        templates
        structures

    main:
        ALIGNMENT(input_fasta, templates, structures)
        GENERATE_SEQ_ME_TREES(ALIGNMENT.out.fasta_aln)
        GENERATE_SEQ_ML_TREES(ALIGNMENT.out.fasta_aln)
        GENERATE_IMD_ME_TREES(ALIGNMENT.out.fasta_aln, templates, structures)
        if (params.align == 'mTMalign'){       
            GENERATE_TM_ME_TREES(ALIGNMENT.out.tmscore_matrix)
        }

}