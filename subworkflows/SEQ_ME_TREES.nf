include { converting_to_phylip } from '../modules/converting_to_phylip.nf'
include { computing_Seq-ME_trees } from '../modules/computing_seq-me_trees.nf'
include { computing_Seq-ME_trees_no_bs } from '../modules/computing_seq-me_trees.nf'


workflow GENERATE_SEQ_ME_TREES{
    take:
        fasta_aln
    main:
    converting_to_phylip(fasta_aln)
    computing_Seq-ME_trees(converting_to_phylip.out)
    
    emit:
    computing_Seq-ME_trees.out
}

workflow GENERATE_SEQ_ME_NO_BS_TREES{
    take:
        fasta_aln
    main:
    converting_to_phylip(fasta_aln)
    computing_Seq-ME_trees_no_bs(converting_to_phylip.out)
    
    emit:
    computing_Seq-ME_trees_no_bs.out
}