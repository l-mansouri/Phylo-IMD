process "collecting_replicates"{
  tag "${id}_${params.type}"
  publishDir "${params.output}/all_tree_files/", mode: 'copy', overwrite: true

  input:
  tuple val(id), path(tree), path (replicates)
  
  output:
  tuple val(id), path("*.trees"), emit: trees

  script:
  """
    cat ${tree} > ${id}_${params.type}.trees
    echo '' >> ${id}_${params.type}.trees
    cat ${replicates} >> ${id}_${params.type}.trees
    sed -i '/^\$/d' ${id}_${params.type}.trees
  """
}
