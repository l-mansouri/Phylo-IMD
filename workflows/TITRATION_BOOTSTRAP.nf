include { ALIGNMENT } from '../subworkflows/alignment.nf'
include { SELECTING_RANDOM_POSITIONS } from '../subworkflows/sampling_for_titration.nf'
include { GENERATE_SEQ_ME_TREES } from '../subworkflows/SEQ_ME_TREES.nf'
include { GENERATE_SEQ_ML_TREES } from '../subworkflows/SEQ_ML_TREES.nf'
include { GENERATE_IMD_ME_TITR_BS_TREES } from '../subworkflows/IMD_ME_TREES.nf'


workflow titration_Phylo_IMD{
    take:
        input_fasta
        templates
        structures

    main:
        ALIGNMENT(input_fasta, templates, structures)
        SELECTING_RANDOM_POSITIONS(ALIGNMENT.out.fasta_aln)
        
        SELECTING_RANDOM_POSITIONS.out.randomized_msa
                        //changing the channel from [id, file] to [id, replicate_number, number_of_columns, file]
                        .map{it.id, it.file.baseName.split('.')[1].split('_')[0], it.file.basName.split('.')[1].split('_')[2], it.file}
                        //filtering to select the replicate
                        .filter(tuple -> tuple[1] == ${params.replicate} )
                        //filtering to select the number of columns
                        .filter(tuple ->tuple[2] == ${params.nCOL})
                        //removing the extra info
                        .map{ tuple -> [tuple[0], tuple[3]] }
                        .set{selected_replicate_msa}

        SELECTING_RANDOM_POSITIONS.out.randomized_col
                        //changing the channel from [id, file] to [id, replicate_number, number_of_columns, file]
                        .map{it.id, it.file.baseName.split('.')[2].split('_')[0], it.file.basName.split('.')[2].split('_')[2], it.file}
                        //filtering to select the replicate
                        .filter(tuple -> tuple[1] == ${params.replicate} )
                        //filtering to select the number of columns
                        .filter(tuple ->tuple[2] == ${params.nCOL})
                        //removing the extra info
                        .map{ tuple -> [tuple[0], tuple[3]] }
                        .set{selected_replicate_col}

        GENERATE_SEQ_ME_TREES(selected_replicate_msa)
        GENERATE_SEQ_ML_TREES(selected_replicate_msa)
        GENERATE_IMD_ME_TITR_BS_TREES(ALIGNMENT.out.fasta_aln, selected_replicate_col, templates, structures)       

}