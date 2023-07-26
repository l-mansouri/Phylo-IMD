include { tcoffee_alignment } from '../modules/tcoffee_alignment.nf'

workflow TC_ALIGNMENT{
    take:
        input_fasta
        templates
        structures

    main:

        tcoffee_alignment(input_fasta)
        fasta_aln=tcoffee_alignment.out.fasta_aln

    emit:
        fasta_aln

}
    
