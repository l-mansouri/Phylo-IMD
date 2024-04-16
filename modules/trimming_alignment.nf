process 'trimming_aln' {
    tag"${id}"
    publishDir "${params.output}/${params.align}_trimmed_fasta/", mode: 'copy', overwrite: true, pattern: "*.fa"
    container 'lmansouri/phylo_imd_trimal:1.0'

    input:
        tuple val(id), path(fasta)

    output:
        tuple val(id), file("*.fa"), emit: trimmed_fasta
    
    script:
    """
    trimal -in ${fasta} -out ${id}_${params.align}_${params.trimmer}.fa -automated1
    """
}