process 'mTMalign_alignment' {

  container 'lmansouri/phylo_imd_mtmalign:1.0'
  tag "${id}"
  publishDir "${params.output}/${params.align}_fasta/$id", pattern: "*.fa", overwrite: false
  publishDir "${params.output}/mTMalign_matrix/$id", pattern: "*.matrix", overwrite: false

  input:
    tuple val(id), path(inputs), path(pdb)

  output:
    tuple val(id), path("*.fa"), emit: fasta_aln
    tuple val(id), path("*.matrix"), emit: tmscore_matrix,  optional: true
  
  script:
  """
  mTM-align -i ${inputs}
  mv mTM_result/result.fasta ./${id}_${params.align}.fa
  mv mTM_result/infile ./${id}_${params.align}.mat
  sed -i "s/${id}.//g" ${id}_${params.align}.fa
  sed -i "/^\$/d" ${id}_${params.align}.fa
  python ${baseDir}/bin/mat_modification.py ${id}_${params.align}.fa ${id}_${params.align}.mat ${id}_${params.align}.matrix
  """

}
