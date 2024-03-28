process 'converting_TM_to_fastme'{
  container 'lmansouri/phylo_imd_base:1.0'

  input:
  tuple val(id), path(tm_mat) 

  output:
  tuple val(id), path("*4_fastme.matrix"), emit: matTM_4fastme

  script:
  """
  sed -E 's/^([^[:space:]]{10})*/\\1 /' ${tm_mat} >${id}_tmscore_4_fastme.matrix
  """
}

process 'computing_TM_ME_trees' {
  tag"${id}"
  publishDir "${params.output}/${params.align}_TM_${params.trimmer}_trees", mode: 'copy', overwrite: true
  container 'lmansouri/phylo_imd_fastme:1.0'

  //errorStrategy 'ignore'

  input:
  tuple val(id), path(mat_tm) 

  output:
  path("*.nwk"), emit: tr_TM_ME 

  script:
  """
  fastme -i ${mat_tm} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  mv ${id}*.nwk ${id}_${params.align}_TM_ME_${params.trimmer}.nwk
  """
}
