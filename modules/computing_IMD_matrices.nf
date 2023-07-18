process 'computing_IMD_matrices' {
 //errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/IMD_matrices/", mode: 'copy', overwrite: true

  input:
    tuple val(id), path(fasta), path(template), path(pdb)


  output:
    tuple val(id), path("*.matrices") emit: matrixOut

  script:
  """
    export THREED_TREE_MODE=${params.tree_mode}
    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +replicates  ${params.replicatesNum} +phylo3d +print_replicates -output dm > ${id}_${params.align}_phylo_IMD.matrices
    sed -i -E  's/^([^[:space:]]{1,10})[^[:space:]]*/\\1/' ${id}_${params.align}_phylo_IMD.matrices
  """
}