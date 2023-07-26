include { sap_tmalign_alignment } from '../modules/3dcoffee_alignment.nf'


workflow STM_ALIGNMENT{
    take:
        input_fasta
        templates
        structures

    main:

        input_fasta.combine(templates, by:0).combine(structures, by:0).set{to_align}
        sap_tmalign_alignment(to_align)
        fasta_aln=sap_tmalign_alignment.out.fasta_aln

    emit:
        fasta_aln


}
    
