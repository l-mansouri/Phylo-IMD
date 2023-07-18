process 'converting_TM_to_fastme'{
  input:
  tuple val(id), path(tm_mat) 

  output:
  tuple val(id), path("*4_fastme.matrix"), emit: matTM_4fastme

  script:
  """
    sed -E 's/^([^[:space:]]{10})[^[:space:]]*/\\1 /' ${tm_mat} >${id}_tmscore_4_fastme.matrix
  """
}

process 'computing_TM-ME_trees' {
  tag"${id}"
  publishDir "${params.output}/TM-ME_trees", mode: 'copy', overwrite: true
  //errorStrategy 'ignore'

  input:
  tuple val(id), path(mat-tm) 

  output:
  path("*.nwk"), emit tr_TM-ME 

  script:
  """
  fastme -i ${mat-tm} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  """
}
