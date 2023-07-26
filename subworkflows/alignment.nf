include { mTMalign_alignment } from '../modules/mtmalign_alignment.nf'
include { sap_tmalign_alignment } from '../modules/3dcoffee_alignment.nf'
include { tcoffee_alignment } from '../modules/tcoffee_alignment.nf'


workflow ALIGNMENT{
    take:
        input_fasta
        templates
        structures

    main:
        if (params.align=='mTMalign'){
            input_fasta.combine(structures, by:0).set{to_align}
            mTMalign_alignment(to_align)
            fasta_aln=mTMalign_alignment.out.fasta_aln
            tmscore_matrix=mTMalign_alignment.out.tmscore_matrix

            }
        else if (params.align=='3dcoffee'){
            input_fasta.combine(templates, by:0).combine(structures, by:0).set{to_align}
            sap_tmalign_alignment(to_align)
            fasta_aln=sap_tmalign_alignment.out.fasta_aln

        }
        else if (params.align=='tcoffee'){
            tcoffee_alignment(input_fasta)
            fasta_aln=tcoffee_alignment.out.fasta_aln

        }
    emit:
        fasta_aln
        //tmscore_matrix

}
    
