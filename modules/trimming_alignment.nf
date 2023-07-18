process 'trimming_aln' {
    tag"${id}"
    publishDir "${params.output}/TRIMMED_fasta/", mode: 'copy', overwrite: true, pattern: "*.fa"

    input:
        tuple val(id), path(fasta)

    output:
        set val(id), file("*.fa"), emit: trimmed_fasta
    
    script:
    """
        trimal -in ${fasta} -out ${id}_${params.align}_trimmal_aln.fa -automated1
    """
}