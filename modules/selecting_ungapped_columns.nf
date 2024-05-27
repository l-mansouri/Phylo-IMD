process 'removing_gaps'{
    tag "${id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'quay.io/biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

    input:
    tuple val(id), path(fasta_aln)

    output:
    tuple val(id), env(l), path("*_ungapped_5_columns.fa")

    script:
    """
    t_coffee -other_pg seq_reformat -in ${fasta_aln} -action +rm_gap 5 -output statistics > ${id}_ungapped_5_columns.stat
    l=\$(head -2 ${id}_ungapped_5_columns.stat| tail -1|sed 's/[[:blank:]]\\+/ /g'|cut -d ' ' -f 3)
    t_coffee -other_pg seq_reformat -in ${fasta_aln} -action +rm_gap 5 >  ${id}_ungapped_5_columns.fa
    """
}
