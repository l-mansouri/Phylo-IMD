include { mTMalign_alignment } from '../modules/mtmalign_alignment.nf'



workflow MTM_ALIGNMENT{
    take:
        input_fasta
        templates
        structures

    main:

        input_fasta.combine(structures, by:0).set{to_align}
        mTMalign_alignment(to_align)
        fasta_aln=mTMalign_alignment.out.fasta_aln
        tmscore_matrix=mTMalign_alignment.out.tmscore_matrix

    emit:
        fasta_aln
        tmscore_matrix

}
    
