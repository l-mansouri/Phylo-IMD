process 'extracting_matrices' {
   tag"${id}"
   publishDir "${params.output}/IMD_matrices/single_matrix", mode: 'copy', overwrite: true

   input:
      tuple val(id), path(matrices) 

   output:
      path("*.txt"), emit: splitMatrix

   script:
   """
   awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
   """
}