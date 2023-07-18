include { mTM-align_alignment } from '../modules/mtmalign_alignment.nf'
include { 3dcoffee_alignment } from '../modules/3dcoffee_alignment.nf'
include { tcoffee_alignment } from '../modules/tcoffee_alignment.nf'


workflow ALIGNMENT{
    take:
        input_fasta
        templates
        structures

    main:
        if(params.align=='mTMalign'){
            input_fasta.combine(structures, by:0).set{to_align}
            mTM-align_alignment(to_align)
            fasta_aln=mTM-align_alignment.out.fasta_aln
            tmscore_matrix=mTM-align_alignment.out.tmscore_matrix
            }
        else if (params.align=='3dcoffee'){
            input_fasta.combine(templates, by:0).combine(structures, by:0).set{to_align}
            3dcoffee_alignment(to_align)
            fasta_aln=3dcoffee_alignment.out.fasta_aln
        }
        else if (params.align=='tcoffee'){
            tcoffee_alignment(input_fasta)
            fasta_aln=tcoffee_alignment.out.fasta_aln
        }
    
    emit:
    fasta_aln
    if (params.align=='mTMalign'){
        tmscore_matrix
    }
}
    
