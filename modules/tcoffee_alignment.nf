process 'tcoffee_alignment' {
    
    tag"${id}"
    container 'lmansouri/phylo_imd_tcoffee:1.0'


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
