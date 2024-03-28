include { collecting_replicates } from '../modules/collecting_replicates.nf'
include { computing_splitfiles } from '../modules/computing_splitfiles.nf'

workflow compute_splits{
    
    take: 
    replicate_trees

    main:

    collecting_replicates(replicate_trees)
    computing_splitfiles(collecting_replicates.out.trees)

    emit: 
    splits =  computing_splitfiles.out.splits        
}