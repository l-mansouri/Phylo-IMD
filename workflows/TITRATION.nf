include { MTM_ALIGNMENT } from '../subworkflows/MTM_alignment.nf'
include { STM_ALIGNMENT } from '../subworkflows/STM_alignment.nf'
include { TC_ALIGNMENT } from '../subworkflows/TC_alignment.nf'
include { SELECTING_RANDOM_POSITIONS } from '../subworkflows/sampling_for_titration.nf'
include { GENERATE_SEQ_ME_NO_BS_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_NO_BS_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TITR_TREES } from '../subworkflows/IMD_ME_TREES.nf'


workflow 'titration_Phylo_IMD' {
    take:
        input_fasta
        templates
        structures

    main:
        if (params.align=='mTMalign'){
            MTM_ALIGNMENT(input_fasta, templates, structures)
            alignment=MTM_ALIGNMENT.out.fasta_aln
            SELECTING_RANDOM_POSITIONS(MTM_ALIGNMENT.out.fasta_aln)
            }
        else if (params.align=='sap_tmalign'){
            STM_ALIGNMENT(input_fasta, templates, structures)
            alignment=STM_ALIGNMENT.out.fasta_aln
            SELECTING_RANDOM_POSITIONS(STM_ALIGNMENT.out.fasta_aln)
        }
        else if (params.align=='tcoffee'){
            TC_ALIGNMENT(input_fasta, templates, structures)
            alignment=TC_ALIGNMENT.out.fasta_aln
            SELECTING_RANDOM_POSITIONS(TC_ALIGNMENT.out.fasta_aln)
        }
        random_msa= SELECTING_RANDOM_POSITIONS.out.randomized_msa.flatten().map{ item -> [ item.simpleName.split('_')[0], item] }
        random_col= SELECTING_RANDOM_POSITIONS.out.randomized_col.flatten().map{ item -> [ item.simpleName.split('_')[0], item] }

        GENERATE_SEQ_ME_NO_BS_TREES(random_msa)
        GENERATE_SEQ_ML_NO_BS_TREES(random_msa)
        GENERATE_IMD_ME_TITR_TREES(alignment, random_col, templates, structures)       

}