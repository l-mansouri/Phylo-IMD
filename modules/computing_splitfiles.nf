process computing_splitfiles {
  tag "${id}_${type}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
      'quay.io/biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

  input:
  tuple val(id), val(type), path(tree)
  
  
  output:
  tuple val(id), val(type), path("*_split"), emit: splits

  script:
  """
  t_coffee -other_pg seq_reformat -in ${tree} -input treelist -action +treelist2splits |grep SPLIT1 |grep original > ${id}_${type}_split
  """
}

