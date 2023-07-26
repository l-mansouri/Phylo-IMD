include { removing_gaps } from '../modules/selecting_ungapped_columns.nf'
include { generating_randomized_fractions } from '../modules/random_sampling_of_the_aln.nf'

workflow SELECTING_RANDOM_POSITIONS{
    take:
    fasta_aln

    main:

    removing_gaps(fasta_aln)
    
    def lenFilter = { tuple -> tuple[1].toInteger() >= 200 }
    def removelen = { tuple -> [tuple[0], tuple[2]] }
    // Use the filter method to keep only the tuples where lenght of the remaining sequences is >= 200
    removing_gaps.out.filter(lenFilter).map(removelen).set{ungapped_alns}

    ungapped_alns.combine(fasta_aln, by:0).set {to_sample}
    generating_randomized_fractions(to_sample)

    emit:
    ungapped_alns
    randomized_msa = generating_randomized_fractions.out.randomized_fa
    randomized_col = generating_randomized_fractions.out.randomized_col
}