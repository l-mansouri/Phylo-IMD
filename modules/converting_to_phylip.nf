process 'converting_to_phylip' {
  tag"${id}"
  container 'lmansouri/phylo_imd_tcoffee:1.0'

  input:
    tuple val(id), path(fasta) 

  output:
    tuple val(id),path("*.ph") , emit: phylip_aln

  script:
  """
  t_coffee -other_pg seq_reformat -in ${fasta} -output phylip_aln > ${fasta.baseName}.ph
  """

}