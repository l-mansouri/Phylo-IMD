process 'computing_IMD_titr_matrices_bs' {
 //errorStrategy 'ignore'
  tag"${id}"
  container 'lmansouri/phylo_imd_tcoffee:1.0'

  input:
    tuple val(id), path(fasta), path(template), path(pairs), path(pdb)


  output:
    tuple val(id), path("*.matrices"), emit: matrixOut

  script:
  """
    export THREED_TREE_MODE=${params.tree_mode}
    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pairs} +replicates ${params.replicatesNum} +print_replicates +phylo3d  -output dm > ${pairs}.matrices 
    sed -i -E  's/^([^[:space:]]{1,10})[^[:space:]]*/\\1/' ${pairs}.matrices 
  """
}
