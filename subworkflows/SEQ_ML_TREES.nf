include { converting_to_phylip } from '../modules/converting_to_phylip.nf'
include { computing_Seq_ML_trees } from '../modules/computing_seq-ml_trees.nf'
include { computing_Seq_ML_trees_no_bs } from '../modules/computing_seq-ml_trees.nf'


workflow GENERATE_SEQ_ML_TREES{
    take:
        fasta_aln
    main:
    converting_to_phylip(fasta_aln)
    computing_Seq_ML_trees(converting_to_phylip.out.phylip_aln)
    
    emit:
    tree       = computing_Seq_ML_trees.out.tr_Seq_ML
    replicates = computing_Seq_ML_trees.out.rep_Seq_ML
}

workflow GENERATE_SEQ_ML_NO_BS_TREES{
    take:
        fasta_aln
    main:
    converting_to_phylip(fasta_aln)
    computing_Seq_ML_trees_no_bs(converting_to_phylip.out.phylip_aln)
    
    // emit:
    // computing_Seq_ML_trees_no_bs.out
}

