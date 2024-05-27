process 'sap_tmalign_alignment' {
    
  tag "${id}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
      'quay.io/biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"
  
  input:
    tuple val(id), path(fasta), path(template), path(pdb)
  
  output:
    tuple val(id), file("*.fa"), emit: fasta_aln

  script:
  """
  t_coffee -seq ${fasta} -template_file ${template} -method sap_pair TMalign_pair -output fasta_aln -outfile ${id}_${params.align}.fa
  sed -i "/^\$/d" ${id}_${params.align}.fa
  # remove everything in each header after first space 
  sed -i 's/ .*//g' ${id}_${params.align}.fa
  """
}
