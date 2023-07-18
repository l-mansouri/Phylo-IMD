process 'computing_IMD-ME_trees' {
  tag"${id}"
  publishDir "${params.output}/IMD_trees/", mode: 'copy', overwrite: true
  errorStrategy 'ignore'

  input:
  tuple val(id), path(mat) 

  output:
  path("*.nwk"), emit tr_IMD-ME

  script:
  """
  fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  """
}