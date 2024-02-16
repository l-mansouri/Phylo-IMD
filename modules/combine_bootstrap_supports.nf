process 'COMBINE_BOOTSTRAP_SUPPORTS' {
    
  tag "${id}"
  publishDir "${params.output}/multistrap/$id" , mode: 'copy', overwrite: true
  container 'luisas/r_multistrap:1.0'
  
  input:
    tuple val(id), path(tree_1), path(replicates_1), path(tree_2), path(replicates_2)
  
  output:
    tuple val(id), file("*_multistrap_bs.nwk"), emit: bootstrap_supports

  script:
  """
  combine_bootstrap_supports.R --t1 ${tree_1}\
                                       --r1 ${replicates_1}\
                                       --t2 ${tree_2}\
                                       --r2 ${replicates_2}\
                                       --bs1 ${id}_1_bs.nwk\
                                       --bs2 ${id}_2_bs.nwk\
                                       --o ${id}_multistrap_bs.nwk
  """
}
