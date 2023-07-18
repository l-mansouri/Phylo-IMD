leilaDir = "${baseDir}/data"


params.mode = "10"   //1-9
params.replicatesNum="1"
params.align='mTMalign'

// params.fastafile="${leilaDir}/*.fasta"
// params.templates="${leilaDir}/*.templ*"
params.templates="${leilaDir}/mtmalign_template_lists/*.templ*"
params.inputlists="${leilaDir}/mtmalign_input_lists/PF*"
params.pdb="${leilaDir}/pdbs/*.pdb"

params.gammaRate="1.0"
params.seedValue="5"

params.output = "${baseDir}/results"

if ( params.inputlists ) {
  Channel
  .fromPath(params.inputlists)
  .map { item -> [ item.baseName.split('_')[0] , item] }
  .into { input_MTM_Ch ; input_MTM}
}
// 
// if ( params.fastafile ) {
//   Channel
//   .fromPath(params.fastafile)
//   .map { item -> [ item.baseName , item] }
//   .into { fastaSeqs ; fastas}
// }
 
if ( params.templates ) {
  Channel
  .fromPath(params.templates)
  .map { item -> [ item.baseName.split('_')[0] , item] }
  .into { templates ; templ_ch; templates2}
}

if ( params.pdb ) {
  Channel
  .fromPath(params.pdb)
  .map { item -> [ item.simpleName , item] }
  .groupTuple()
  //.collect()
  .into { pdbFiles ; pdb_ch; pdbFiles2}
}

input_MTM_Ch
        //.combine(templates, by:0)
        .combine(pdbFiles, by:0)
        .set{aln_ch}  

process 'running alignment' {

   publishDir "${params.output}/msa_fasta", mode: 'copy', overwrite: true, pattern: "*.fa"
   publishDir "${params.output}/mTMalign_matrix", mode: 'copy', overwrite: true, pattern: "*.matrix"
    
   input:
     set val(id), file(inputs), file(pdb) from aln_ch


   output:
    set val(id), file("*.fa") into fasta_aln, fasta_aln2, fasta_aln3, original_msa_ch, oriSeqs
    set val(id), file("*.matrix") into mtmalign_matrix_ch
   script:
   """
    /users/cn/lmansouri/bin/mTM-align -i ${inputs}
    mv mTM_result/result.fasta ./${id}_mTMalign.fa
    mv mTM_result/infile ./${id}_mTMalign.matrix
    sed -i "s/${id}//g" ${id}_mTMalign.fa
   """

}

process 'converting 1d format' {
  tag"${id}"
  publishDir "${params.output}/msa_ph", mode: 'copy', overwrite: true

  input:
    set val(id), file(fasta) from  fasta_aln

  output:
    set val(id),file("*.ph") into phylip_aln, phy_aln

  script:
  """
    t_coffee -other_pg seq_reformat -in ${fasta} -output phylip_aln > ${id}_mTMalign.ph
  """

}


process 'computing trees on 1d' {
    tag"${id}"
    publishDir "${params.output}/ME_trees/untrimmed", mode: 'copy', overwrite: true

    input:
    set val(id),file(phylip) from phylip_aln

    output:
    set val(id),file("*.nwk") into unttree1d
    set val(id), file("*.replicates") into untreplicate1d

    script:
    """
    /users/cn/lmansouri/bin/fastme -i ${phylip} -o ${id}_tmalign.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b 100 -B ${id}_tmalign.1dtree.replicates
    """
}

process 'computing_ML_trees_untrimmed'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/ML_trees/untrimmed/", mode: 'copy', overwrite: true
  
  input:
  set val(id), file(phylip) from phy_aln

  output:
  set val(id), file("*.treefile"), file("*.boottrees") into unt_iqtree_ch

  script:
  """
  /users/cn/abaltzis/bin/iqtree2 -s ${phylip} -b ${params.replicatesNum} 
  """
}

process 'computing_NJ_tree'{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/NJ_trees/untrimmed/", mode: 'copy', overwrite: true
  
  input:
  set val(id), file(matrix) from mtmalign_matrix_ch

  output:
  set val(id), file("*.nwk") into unt_NJtree_ch

  script:
  """
  echo "${id}_mTMalign.matrix" > NJ_in
  echo "J" >> NJ_in
  echo "5" >> NJ_in
  echo "Y" >> NJ_in  
  /nfs/users/cn/lmansouri/tcoffee_extras/plugins/linux/neighbor < NJ_in
  mv outtree ${id}_NJ_tree.nwk
  """
}

oriSeqs
       .combine(templates, by:0)
       .combine(pdb_ch, by:0)
       .set{tc_mat_ch}


process 'computing_3d_matrices_untrimmed' {
 //errorStrategy 'ignore'
 tag"${id}"
 publishDir "${params.output}/3d_matrices/untrimmed/", mode: 'copy', overwrite: true

 input:
   set val(id), file(fasta), file(template), file(pdb) from tc_mat_ch


 output:
   set val(id),file("*.matrices") into untmatrixOut

 script:
 """
   export THREED_TREE_MODE=${params.mode}
   t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +replicates  ${params.replicatesNum} +phylo3d +print_replicates -output dm > ${id}_phylo3d.matrices 2>${id}_phylo3d_mat_err_log
    sed -i -E  's/^([^[:space:]]{1,10})[^[:space:]]*/\\1/'  ${id}_phylo3d.matrices
 """
}


process 'extr_matrices_per_family_untrimmed' {
   tag"${id}"
   publishDir "${params.output}/3d_matrices/untrimmed/single_matrices", mode: 'copy', overwrite: true

   input:
    set val(id), file(matrices) from untmatrixOut

   output:
    file("*.txt") into untsplitMatrix

   script:
   """
   awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
   """
}

untsplitMatrix
  .flatten()
  .map { item -> [ item.simpleName, item] }
  .set{untsplitMatrix2}


process 'fastme on 3d matrices_untrimmed' {
  tag"${id}"
  publishDir "${params.output}/3d_trees/untrimmed", mode: 'copy', overwrite: true
  errorStrategy 'ignore'

  input:
  set val(id), file(mat) from untsplitMatrix2

  output:
  file("*.nwk") into unttrees3d

  script:
  """
  /users/cn/lmansouri/bin/fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
  """
}

workflow.onComplete {
println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}

