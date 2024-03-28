process "collecting_replicates"{
  tag "${id}_${type}"
  publishDir "${params.output}/all_tree_files/", mode: 'copy', overwrite: true

  input:
  tuple val(id), path(tree), path (replicates), val(type)
  
  
  output:
  tuple val(id), val(type), path("*.trees"), emit: trees

  script:
  """
    cat ${tree} > ${id}_${type}.trees
    echo '' >> ${id}_${type}.trees
    cat ${replicates} >> ${id}_${type}.trees
    sed -i '/^\$/d' ${id}_${type}.trees
  """
}
