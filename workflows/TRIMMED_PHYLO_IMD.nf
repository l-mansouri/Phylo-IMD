include { MTM_ALIGNMENT } from '../subworkflows/MTM_alignment.nf'
include { STM_ALIGNMENT } from '../subworkflows/STM_alignment.nf'
include { TC_ALIGNMENT } from '../subworkflows/TC_alignment.nf'
include { TRIMMING_ALIGNMENT as TRIMMING_ALIGNMENT1 } from '../subworkflows/trimming_aln.nf'
include { TRIMMING_ALIGNMENT as TRIMMING_ALIGNMENT2 } from '../subworkflows/trimming_aln.nf'
include { TRIMMING_ALIGNMENT as TRIMMING_ALIGNMENT3 } from '../subworkflows/trimming_aln.nf'
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
            TRIMMING_ALIGNMENT1(MTM_ALIGNMENT.out.fasta_aln)
            ALIGNMENT = MTM_ALIGNMENT.out.fasta_aln
            TRIMMED_ALIGNMENT = TRIMMING_ALIGNMENT1
        }
        else if (params.align=='sap_tmalign'){
            STM_ALIGNMENT(input_fasta, templates, structures)
            TRIMMING_ALIGNMENT2(STM_ALIGNMENT.out.fasta_aln)
            ALIGNMENT = STM_ALIGNMENT.out.fasta_aln
            TRIMMED_ALIGNMENT = TRIMMING_ALIGNMENT2
        }
        else if (params.align=='tcoffee'){
            TC_ALIGNMENT(input_fasta, templates, structures)
            TRIMMING_ALIGNMENT3(TC_ALIGNMENT.out.fasta_aln)
            ALIGNMENT = TC_ALIGNMENT.out.fasta_aln
            TRIMMED_ALIGNMENT = TRIMMING_ALIGNMENT3
        }

        GENERATE_SEQ_ME_TREES(TRIMMED_ALIGNMENT.out.trimmed_aln)
        GENERATE_SEQ_ML_TREES(TRIMMED_ALIGNMENT.out.trimmed_aln)
        GENERATE_IMD_ME_TR_TREES(ALIGNMENT, TRIMMED_ALIGNMENT.out.mapped_pos, templates, structures)       

}