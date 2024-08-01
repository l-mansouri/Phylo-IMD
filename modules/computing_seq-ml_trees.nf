process 'computing_Seq_ML_trees'{
  errorStrategy 'ignore'
  tag"${id}"
  container 'lmansouri/phylo_imd_iqtree:1.0'

  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.nwk"), emit: tr_Seq_ML
  tuple val(id), path("*.replicates"), emit: rep_Seq_ML

  script:
  """
  iqtree -s ${phylip} -b ${params.replicatesNum}
  if [ ! -f *_columns.ph ]; then
    mv *.treefile ${id}_${params.align}_${params.trimmer}_ML.nwk
    mv *.boottrees ${id}_${params.align}_${params.trimmer}_ML.replicates
  else
    for file in *treefile; do
      new_file="\${file%.treefile}.nwk"
      mv \${file} \${new_file}
    done

    for file in *boottrees; do
      new_file="\${file%.boottrees}.replicates"
      mv \${file} \${new_file}
    done
  fi
  """
}

process 'computing_Seq_ML_trees_no_bs'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/${params.align}_ML_${params.trimmer}_trees/", mode: 'copy', overwrite: true, pattern: "*.nwk"
  container 'lmansouri/phylo_imd_iqtree:1.0'

  
  input:
  tuple val(id), path(phylip)

  output:
  tuple val(id), path("*.treefile"), emit: tr_Seq_ML


  script:
  """
  iqtree -s ${phylip} 
  if [ ! -f *_columns.ph ]; then
    mv *.treefile ${id}_${params.align}_${params.trimmer}_ML.treefile
  fi
  """
}