process 'extracting_matrices' {
   tag"${id}"
   container 'lmansouri/phylo_imd_base:1.0'


   input:
      tuple val(id), path(matrices) 

   output:
      path("*.matrix"), emit: main
      path("*.txt"), emit: replicates

   script:
   """
   awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
   if [ -f *columns*.txt ]; then
      mv ${id}*matrices.1.txt ${id}_${params.align}_3d_${params.trimmer}.matrix
   else 
      for file in *matrices.1.txt; do
         new_file="\${file%.txt}.matrix"
         # Rename the file
         mv \${file} \${new_file}
         echo "File \${file} renamed to \${new_file}" 
      done
   fi
   """
}