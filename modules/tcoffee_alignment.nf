process 'tcoffee_alignment' {
    
    tag"${id}"
    publishDir "${params.output}/msa_fasta", mode: 'copy', overwrite: true
    container 'lmansouri/phylo_imd_tcoffee:1.0'


    input:
        tuple val(id), path(fasta)

    output:
        tuple val(id), path("*.fa"), emit: fasta_aln
    
    script:
    """
        t_coffee -in ${fasta} -output=fasta_aln >${id}_${params.align}.clustal
        mv ${id}.fasta_aln ${id}_${params.align}.fa
    """
}