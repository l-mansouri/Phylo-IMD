process 'computing_IMD_matrices' {
  tag"${id}"
  container 'lmansouri/phylo_imd_tcoffee:1.0'

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