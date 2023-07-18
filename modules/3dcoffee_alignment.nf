process '3dcoffee_alignment' {
    
  tag "${id}"
  publishDir "${params.output}/msa_fasta" , mode: 'copy', overwrite: true
  
  input:
    tuple val(id), path(fasta), path(template), path(pdb)
  
  output:
    tuple val(id), file("*.fa"), emit: fasta_aln

  script:
  """
  t_coffee -seq ${fasta} -template_file ${template} -method sap_pair TMalign_pair -output fasta_aln -outfile ${id}_${params.align}.fa
  """
}