process 'computing_IMD_ME_trees' {
  tag"${id}"
  publishDir "${params.output}/${params.align}_3d_ME_${params.trimmer}_trees/$replicate", mode: 'copy', overwrite: true
  errorStrategy 'ignore'
  container 'lmansouri/phylo_imd_fastme:1.0'

  input:
  tuple val(id), path(mat)
  val(replicate) 

  output:
  tuple val(id), path("*.nwk"), emit: tr_IMD_ME

  script:
  """
  fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  """
}