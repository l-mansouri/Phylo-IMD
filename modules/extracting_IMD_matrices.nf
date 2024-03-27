process 'extracting_matrices' {
   tag"${id}"
   publishDir "${params.output}/${params.align}_3d_ME_${params.trimmer}_matrix/", mode: 'copy', overwrite: true, pattern: "*.matrix"
   publishDir "${params.output}/${params.align}_3d_ME_${params.trimmer}_matrix/replicates", mode: 'copy', overwrite: true, pattern: "*.txt"
   container 'lmansouri/phylo_imd_base:1.0'


   input:
      tuple val(id), path(matrices) 

   output:
      path("*.matrix"), emit: main
      path("*.txt"), emit: replicates

   script:
   """
   awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}

   mv ${id}*matrices.1.txt ${id}_${params.align}_3d_${params.trimmer}.matrix
   """
}