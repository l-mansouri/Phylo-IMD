include { converting_TM_to_fastme; computing_TM_ME_trees } from '../modules/computing_tm-me_trees.nf'

workflow GENERATE_TM_ME_TREES{
    take:
        tmscore_matrix

    main:
        converting_TM_to_fastme(tmscore_matrix)
        computing_TM_ME_trees(converting_TM_to_fastme.out)
    
    emit:
        computing_TM_ME_trees.out
}