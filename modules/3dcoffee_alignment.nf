process 'sap_tmalign_alignment' {
    
  tag "${id}"
  storeDir "${params.output}/sap_tmalign_${params.trimmer}_fasta/$id" 
  container 'lmansouri/phylo_imd_tcoffee:1.0'
  
  input:
    tuple val(id), path(fasta), path(template), path(pdb)
  
  output:
    tuple val(id), file("*.fa"), emit: fasta_aln

  script:
  """
  t_coffee -seq ${fasta} -template_file ${template} -method sap_pair TMalign_pair -output fasta_aln -outfile ${id}_${params.align}.fa
  sed -i "/^\$/d" ${id}_${params.align}.fa
  """
}
