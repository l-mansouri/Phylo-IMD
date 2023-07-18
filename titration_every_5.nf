leilaDir = "/users/cn/lmansouri/PROJECTS/Phylo3D/NEW_Phylo3D/NF_TMalign/mTMalign"

params.msa = "mTMalign"
params.mode = "10"   //1-11
params.columns="200"

params.msafile="${leilaDir}/ungapped_5/fasta_aln_titration/*_gapped_5_columns.fa"
params.originalmsa="${leilaDir}/msa_fasta/PF*.fa"
params.templates="${leilaDir}/template_lists/PF*"
params.pdb="${leilaDir}/pdb/*"

params.replicatesNum='1'
params.gammaRate="1.0"
params.seedValue="5"

params.output = "${leilaDir}/titration_every_5/"


if ( params.msafile ) {
  Channel
  .fromPath(params.msafile)
  .map { item -> [ item.baseName.replace('_gapped_5_columns','') , item] }
  .into { fastaAln; fastaSeqs ; trmsa}
}

if ( params.originalmsa ) {
  Channel
  .fromPath(params.originalmsa)
  .map { item -> [ item.baseName.replace('_mTMalign','') , item] }
  .into { ori_aln; oriSeqs ; orifa}
}

fastaAln
       .combine(ori_aln, by: 0)
       .set{fastaAln_ch}

if ( params.templates ) {
  Channel
  .fromPath(params.templates)
  .map { item -> [ item.baseName, item] }
  .into { templates ; tmpl}
}

//tmpl.println()

if ( params.pdb ) {
  Channel
  .fromPath(params.pdb)
  .map { item -> [ item.simpleName, item] }
  .groupTuple()
  .into { pdbFiles ; pdbs}
}
//pdbs.println()


process 'generating_randomized_fractions' {
  //generates the randomized column alignment and the randomized column pairs
    tag "${id}"
    publishDir "${params.output}/${id}/raw_files/", mode: 'copy', overwrite: true
    label 'short'

    input:
      set val(id), file(trimmed), file(original) from fastaAln_ch

    output:
      file("*.fa") into rnd1dFra //random aln
      file("*columns.txt") into rndCOL //random column pairs

    script:
    """
      python ${leilaDir}/randomizing_msa_and_fractions_by_5.py ${trimmed} ${original} ${id}
    """
}

rnd1dFra
  .flatten()
  .map{ item -> [ item.simpleName.replace('_random_msa_replicate','') , item] }
  .into { rndmsa;rndmsa2;check}


//PF*.random_column_pairs_replicate.9_with_99_columns.txt
rndCOL
  .flatten()
  .map{ item -> [ item.simpleName, item.baseName.replace(item.simpleName,'').replace('.random_column_pairs_replicate.','').replaceAll('_with_([0-9]+)_columns', '') , item]}
  .into { rndPcol; check1}

//check1.println()


process 'converting_1d_format'{
  tag "${id}"
  publishDir "${params.output}/${id}/phylip_aln", mode: 'copy', overwrite: true
  label 'long'
  input:
    set val(id), file(rfasta) from  rndmsa

  output:
    set val(id), file("*.ph") into rnd1dFrac, rndMLFrac

  script:
  """
    t_coffee -other_pg seq_reformat -in ${rfasta} -output phylip_aln -out ${rfasta.baseName}.ph
  """
}


process computing_trees_on_1d_fractions{
    tag"${id}"
    publishDir "${params.output}/${id}/1d_trees/", mode: 'copy', overwrite: true
    label 'short'

    input:
    set val(id),file(fraphy) from rnd1dFrac

    output:
    set val(id), file("*.nwk") into rnd1dFrT

    script:
    """
    /users/cn/lmansouri/bin/fastme -i ${fraphy} -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}

process computing_ML_trees{
  errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/${id}/ML_trees/", mode: 'copy', overwrite: true
  label 'long'

  input:
  set val(id),file(phylip) from rndMLFrac

  output:
  set val(id), file("*.treefile") into iqtree_ch

  script:
  """
  /users/cn/abaltzis/bin/iqtree2 -s ${phylip}
  """
}

// params.raw_pairs_4mat="${params.output}/PF**/raw_files/*_columns.txt"
//
// if ( params.raw_pairs_4mat ) {
//   Channel
//   .fromPath(params.raw_pairs_4mat)
//   .flatten()
//   .map{ item -> [ item.simpleName, item.baseName.replace(item.simpleName,'').replace('.random_column_pairs_replicate.','').replaceAll('_with_([0-9]+)_columns', '') , item]}
//   .into { rndPcol; check1 }
// }

//check1.view()

rndPcol
  .combine(oriSeqs, by:0)
  .combine(templates, by: 0)
  .combine(pdbs, by:0)
  .into{ ori_3d_mat_ch; maxd_mat_ch }


  //ori_3d_mat_ch.println()
//maxd_mat_ch.view()

process computing_3d_matrices{
  // errorStrategy 'ignore'
  tag"${id}"
  publishDir "${params.output}/${id}/3d_matrices/", mode: 'copy', overwrite: true
  label 'long'

  input:
    set val(id), val(rep),  file(pair), file(fasta), file(template), file(pdb) from ori_3d_mat_ch

  output:
    set val(id),file("*.matrices") into orimatOut

  script:
  """
    export THREED_TREE_MODE=${params.mode};
    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pair} +replicates ${params.replicatesNum} +phylo3d  -output dm > ${pair}.matrices 2> ${pair}_3dtrees_err_log
    #+print_replicates
  """
}

process fastme_matrices_on_pairs{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${id}/3d_trees", mode: 'copy', overwrite: true
    label 'short'
    
    input:
    set val(id), file(mat) from orimatOut

    output:
    file("*.nwk") into matricesOut

    script:
    """
    /users/cn/lmansouri/bin/fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
