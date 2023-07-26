process 'extracting_matrices' {
   tag"${id}"
   publishDir "${params.output}/IMD_matrices/single_matrix", mode: 'copy', overwrite: true
   container 'lmansouri/phylo_imd_base:1.0'


   input:
      tuple val(id), path(matrices) 

   output:
      path("*.txt"), emit: splitMatrix

   script:
   """
   awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
   """
}