include { converting_to_phylip } from '../modules/converting_to_phylip.nf'
include { computing_Seq-ML_trees } from '../modules/computing_seq-ml_trees.nf'

workflow GENERATE_SEQ_ML_TREES{
    take:
        fasta_aln
    main:
    converting_to_phylip(fasta_aln)
    computing_Seq-ML_trees(converting_to_phylip.out)
    
    emit:
    computing_Seq-ML_trees.out
}