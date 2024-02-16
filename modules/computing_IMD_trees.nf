process 'computing_IMD_ME_trees' {
  tag"${id}"
  publishDir "${params.output}/IMD_trees/", mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  container 'lmansouri/phylo_imd_fastme:1.0'

  input:
  tuple val(id), path(mat) 

  output:
  tuple val(id), path("*.nwk"), emit: tr_IMD_ME

  script:
  """
  fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  """
}