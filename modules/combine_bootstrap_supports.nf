process 'COMBINE_BOOTSTRAP_SUPPORTS' {
    
  tag "${id}"
  container 'luisas/r_multistrap:1.0'
  
  input:
    tuple val(id), path(tree_1), path(replicates_1), path(tree_2), path(replicates_2)
  
  output:
    tuple val(id), file("*_bs.nwk"), emit: bootstrap_supports

  script:
  """
  combine_bootstrap_supports.R --t1 ${tree_1}\
                                       --r1 ${replicates_1}\
                                       --t2 ${tree_2}\
                                       --r2 ${replicates_2}\
                                       --bs1 ${id}_${params.seq_tree}_tree_${params.seq_tree}bs.nwk\
                                       --bs2 ${id}_IMD_bs.nwk\
                                       --bs3 ${id}_${params.seq_tree}_multistrap_bs.nwk \
                                       --o ${id}_multistrap_bs.nwk
  """
}
