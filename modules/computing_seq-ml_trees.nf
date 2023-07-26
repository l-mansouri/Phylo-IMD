process 'computing_Seq_ML_trees'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees/", mode: 'copy', overwrite: true
  container 'lmansouri/phylo_imd_iqtree:1.0'

  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.treefile"), emit: tr_Seq_ML
  tuple val(id), path("*.boottrees"), emit: rep_Seq_ML

  script:
  """
  iqtree -s ${phylip} -b ${params.replicatesNum} 
  """
}

process 'computing_Seq_ML_trees_no_bs'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees/", mode: 'copy', overwrite: true
  container 'lmansouri/phylo_imd_iqtree:1.0'

  
  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.treefile"), emit: tr_Seq_ML


  script:
  """
  iqtree -s ${phylip} 
  """
}