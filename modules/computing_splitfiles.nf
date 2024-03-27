process computing_splitfiles {
  tag "${id}_${params.type}"
  publishDir "${params.output}/split_files/${params.type}/", mode: 'copy', overwrite: true

  input:
  tuple val(id), path(tree)
  
  output:
  tuple val(id), path("*_split"), emit: splits

  script:
  """
  t_coffee -other_pg seq_reformat -in ${tree} -input treelist -action +treelist2splits |grep SPLIT1 |grep original > ${id}_${params.type}_split
  """
}

