process 'computing_IMD_tr_matrices' {
 //errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/${params.align}_3d_ME_${params.trimmer}_matrix/", mode: 'copy', overwrite: true
  container 'lmansouri/phylo_imd_tcoffee:1.0'

  input:
    tuple val(id), path(fasta), path(template), path(pairs), path(pdb)


  output:
    tuple val(id), path("*.matrices"), emit: matrixOut

  script:
  """
    export THREED_TREE_MODE=${params.tree_mode}
    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pairs} +replicates ${params.replicatesNum} +print_replicates +phylo3d  -output dm > ${id}_${params.align}_${params.trimmer}.matrices 
    sed -i -E  's/^([^[:space:]]{1,10})[^[:space:]]*/\\1/' ${id}_${params.align}_${params.trimmer}.matrices 
  """
}