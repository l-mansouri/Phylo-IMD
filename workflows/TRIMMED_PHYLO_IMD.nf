include { MTM_ALIGNMENT } from '../subworkflows/MTM_alignment.nf'
include { STM_ALIGNMENT } from '../subworkflows/STM_alignment.nf'
include { TC_ALIGNMENT } from '../subworkflows/TC_alignment.nf'
include { TRIMMING_ALIGNMENT } from '../subworkflows/trimming_aln.nf'
include { GENERATE_SEQ_ME_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TR_TREES } from '../subworkflows/IMD_ME_TREES.nf'


workflow trimmed_Phylo_IMD{
    take:
        input_fasta
        templates
        structures

    main:
        if (params.align=='mTMalign'){
            MTM_ALIGNMENT(input_fasta, templates, structures)
            TRIMMING_ALIGNMENT(MTM_ALIGNMENT.out.fasta_aln)
            }
        else if (params.align=='3dcoffee'){
            STM_ALIGNMENT(input_fasta, templates, structures)
            TRIMMING_ALIGNMENT(STM_ALIGNMENT.out.fasta_aln)

        }
        else if (params.align=='tcoffee'){
            TC_ALIGNMENT(input_fasta, templates, structures)
            TRIMMING_ALIGNMENT(TC_ALIGNMENT.out.fasta_aln)
        }
        TRIMMING_ALIGNMENT(ALIGNMENT.out.fasta_aln)
        GENERATE_SEQ_ME_TREES(TRIMMING_ALIGNMENT.out.trimmed_aln)
        GENERATE_SEQ_ML_TREES(TRIMMING_ALIGNMENT.out.trimmed_aln)
        GENERATE_IMD_ME_TR_TREES(ALIGNMENT.out.fasta_aln, TRIMMING_ALIGNMENT.out.mapped_pos, templates, structures)       

}