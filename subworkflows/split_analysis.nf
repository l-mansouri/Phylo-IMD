
include { computing_splitfiles } from '../modules/computing_splitfiles.nf'
include { collecting_replicates } from '../modules/collecting_replicates.nf'

workflow split_analysis{

    take: 
    trees

    main:
    // collect the bootstrap support
    // compute the splits 
    println "Computing splits"
}