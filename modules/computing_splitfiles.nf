process computing_splitfiles {
  tag "${id}_${type}"
  publishDir "${params.output}/split_files/${type}/", mode: 'copy', overwrite: true
  container 'lmansouri/phylo_imd_tcoffee:1.0'

  input:
  tuple val(id), val(type), path(tree)
  
  
  output:
  tuple val(id), val(type), path("*_split"), emit: splits

  script:
  """
  t_coffee -other_pg seq_reformat -in ${tree} -input treelist -action +treelist2splits |grep SPLIT1 |grep original > ${id}_${type}_split
  """
}

