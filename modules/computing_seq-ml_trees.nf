process 'computing_Seq-ML_trees'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees/", mode: 'copy', overwrite: true
  
  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.treefile"), emit: tr_Seq-ML
  tuple val(id), path("*.boottrees"), emit: rep_Seq-ML

  script:
  """
  iqtree2 -s ${phylip} -b ${params.replicatesNum} 
  """
}

process 'computing_Seq-ML_trees_no_bs'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees/", mode: 'copy', overwrite: true
  
  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.treefile"), emit: tr_Seq-ML
  tuple val(id), path("*.boottrees"), emit: rep_Seq-ML

  script:
  """
  iqtree2 -s ${phylip} 
  """
}