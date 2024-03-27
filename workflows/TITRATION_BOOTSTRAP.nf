include { MTM_ALIGNMENT } from '../subworkflows/MTM_alignment.nf'
include { STM_ALIGNMENT } from '../subworkflows/STM_alignment.nf'
include { TC_ALIGNMENT } from '../subworkflows/TC_alignment.nf'
include { SELECTING_RANDOM_POSITIONS } from '../subworkflows/sampling_for_titration.nf'
include { GENERATE_SEQ_ME_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TITR_BS_TREES } from '../subworkflows/IMD_ME_TREES.nf'


workflow 'titration_bootstrap_Phylo_IMD' {
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
        
        def sel_replicate = { tuple -> tuple[1] == params.replicate }
        def sel_columns ={ tuple ->tuple[2] == params.nCOL }
        def removeinfo = { tuple -> [ tuple[0], tuple[3] ] }

        SELECTING_RANDOM_POSITIONS.out.randomized_msa
                        .flatten()
                        //changing the channel from file to [id, replicate_number, number_of_columns, file]
                        .map{ item -> [item.simpleName.split('_')[0], item.baseName.replace(item.simpleName,'').split("_")[0].replace('.',''), item.baseName.replace(item.simpleName,'').split("_")[2], item ] }
                        //filtering to select the replicate
                        .filter( sel_replicate )
                        //filtering to select the number of columns
                        .filter( sel_columns )
                        //removing the extra info
                        .map( removeinfo )
                        .set{selected_replicate_msa}

        SELECTING_RANDOM_POSITIONS.out.randomized_col
                        .flatten()
                        //changing the channel from file to [id, replicate_number, number_of_columns, file]
                        .map{item -> [item.simpleName.split('_')[0], item.baseName.replace(item.simpleName, '').split('_')[3].replace('replicate.',''), item.baseName.replace(item.simpleName, '').split("_")[5], item] }
                        //filtering to select the replicate
                        .filter( sel_replicate )
                        //filtering to select the number of columns
                        .filter( sel_columns )
                        //removing the extra info
                        .map( removeinfo )
                        .set{selected_replicate_col}

        GENERATE_SEQ_ME_TREES(selected_replicate_msa)
        GENERATE_SEQ_ML_TREES(selected_replicate_msa)
        GENERATE_IMD_ME_TITR_BS_TREES(alignment, selected_replicate_col, templates, structures)   

}