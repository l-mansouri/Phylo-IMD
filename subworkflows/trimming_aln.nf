include { trimming_aln } from '../modules/trimming_alignment.nf'
include { mapping_pos } from '../modules/mapping_position.nf'

workflow TRIMMING_ALIGNMENT{
    take:
    fasta_aln

    main:

    trimming_aln(fasta_aln)
    fasta_aln.combine(trimmed_aln.out, by:0).set{for_mapping}
    mapping_pos(for_mapping)

    emit:
    trimmed_aln=trimming_aln.out
    mapped_pos=mapping_pos.out
}