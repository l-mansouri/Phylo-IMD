process 'converting_to_phylip' {
  tag"${id}"
  publishDir "${params.output}/msa_ph", mode: 'copy', overwrite: true

  input:
    tuple val(id), path(fasta) 

  output:
    tuple val(id),path("*.ph") , emit: phylip_aln

  script:
  """
    t_coffee -other_pg seq_reformat -in ${fasta} -output phylip_aln > ${id}_${params.align}.ph
  """

}