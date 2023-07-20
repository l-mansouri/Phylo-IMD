process 'generating_randomized_fractions' {
  //generates the randomized column alignment and the randomized column pairs
    tag "${id}"
    publishDir "${params.output}/random_sampled_alignment", mode: 'copy', overwrite: true, pattern: "*.fa"
    label 'short'

    input:
        tuple val(id), path(ungapped_aln), path(original_aln) 

    output:
        tuple val(id), path("*.fa"), emit: randomized_fa //random aln
        tuple val(id), path("*columns.txt"), emit: randomized_col //random column pairs

    script:
    """
        python ${baseDir}/bin/randomizing_msa_and_fractions_by_5.py ${ungapped_aln} ${original_aln} ${id} ${params.align}_${params.trimmer}
    """
}