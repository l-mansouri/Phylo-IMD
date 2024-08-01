process 'tcoffee_alignment' {
    
    tag"${id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'quay.io/biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"


    input:
        tuple val(id), path(fasta)

    output:
        tuple val(id), path("*.fa"), emit: fasta_aln
    
    script:
    """
    t_coffee -in ${fasta} -output=fasta_aln >${id}_${params.align}.clustal
    mv ${id}.fasta_aln ${id}_${params.align}.fa
	sed -i "/^\$/d" ${id}_${params.align}.fa
    """
}
