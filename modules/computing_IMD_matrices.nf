process 'computing_IMD_matrices' {
  tag"${id}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
      'quay.io/biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

  input:
    tuple val(id), path(fasta), path(template), path(pdb)


  output:
    tuple val(id), path("*.matrices"), emit: matrixOut

  script:
  """
  export THREED_TREE_MODE=${params.tree_mode}
  t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +replicates  ${params.replicatesNum} +phylo3d +print_replicates -output dm > ${fasta.baseName}.matrices
  sed -i -E  's/^([^[:space:]]{1,10})[^[:space:]]*/\\1/' ${fasta.baseName}.matrices

  if [ -f *columns* ]; then
    mv *.matrices ${id}_${params.align}_phylo_IMD.matrices
  fi
  """
}