process 'removing_gaps'{
    tag="${id}"
    publishDir "${params.output}/UNGAPPED_5_fasta/", mode: 'copy', overwrite: true, pattern: "*.fa"

    input:
    tuple val(id), path(fasta_aln)

    output:
    tuple val(id), env(l), path("*_ungapped_5_columns.fa")

    script:
    """
    t_coffee -other_pg seq_reformat -in ${fasta_aln} -action +rm_gap 5 -output statistics > ${id}_ungapped_5_columns.stat
    l=$(head -2 ${id}_ungapped_5_columns.stat| tail -1|sed 's/[[:blank:]]\\+/ /g'|cut -d ' ' -f 3)
    t_coffee -other_pg seq_reformat -in ${fasta_aln} -action +rm_gap 5 >  ${id}_ungapped_5_columns.fa
    """
}