include { MTM_ALIGNMENT } from '../subworkflows/MTM_alignment.nf'
include { STM_ALIGNMENT } from '../subworkflows/STM_alignment.nf'
include { TC_ALIGNMENT } from '../subworkflows/TC_alignment.nf'
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
        if (params.align=='mTMalign'){
            MTM_ALIGNMENT(input_fasta, templates, structures)
            GENERATE_SEQ_ME_TREES(MTM_ALIGNMENT.out.fasta_aln)
            GENERATE_SEQ_ML_TREES(MTM_ALIGNMENT.out.fasta_aln)
            GENERATE_IMD_ME_TREES(MTM_ALIGNMENT.out.fasta_aln, templates, structures)
            GENERATE_TM_ME_TREES(MTM_ALIGNMENT.out.tmscore_matrix)
            }
        else if (params.align=='3dcoffee'){
            STM_ALIGNMENT(input_fasta, templates, structures)
            GENERATE_SEQ_ME_TREES(STM_ALIGNMENT.out.fasta_aln)
            GENERATE_SEQ_ML_TREES(STM_ALIGNMENT.out.fasta_aln)
            GENERATE_IMD_ME_TREES(STM_ALIGNMENT.out.fasta_aln, templates, structures)
        }
        else if (params.align=='tcoffee'){
            TC_ALIGNMENT(input_fasta, templates, structures)
            GENERATE_SEQ_ME_TREES(TC_ALIGNMENT.out.fasta_aln)
            GENERATE_SEQ_ML_TREES(TC_ALIGNMENT.out.fasta_aln)
            GENERATE_IMD_ME_TREES(TC_ALIGNMENT.out.fasta_aln, templates, structures)
        }

}