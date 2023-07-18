leilaDir = "/users/cn/lmansouri/PROJECTS/Phylo3D/NEW_Phylo3D/NF_TMalign/mTMalign"

params.msa = "mTMalign"
params.mode = "10"   //1-11
params.columns = "200"
params.replicatesNum = "100"
params.gammaRate="1.0"
params.seedValue="5"

params.output = "${leilaDir}/titration_every_5_bootstrap_200_columns/"

params.raw_phylip_ME="${leilaDir}/titration_every_5/PF**/phylip_aln/PF*_random_msa_replicate.0_with_200_columns.ph"
params.raw_phylip_ML="${leilaDir}/titration_every_5/PF**/phylip_aln/PF*_random_msa_replicate.0_with_200_columns.ph"
params.raw_pairs_4mat="${leilaDir}/titration_every_5/PF**/raw_files/*random_column_pairs_replicate.0_with_200_columns.txt"
params.originalmsa="${leilaDir}/msa_fasta/PF*.fa"
params.templates="${leilaDir}/template_lists/PF*"
params.pdb="${leilaDir}/pdb/*"


if ( params.raw_phylip_ME ) {
    Channel
    .fromPath(params.raw_phylip_ME)
    .flatten()
    .map{ item -> [item.baseName.replace('_random_msa_replicate.0_with_200_columns','') , item]}
    .set { rnd1dFrac }
}

if ( params.raw_phylip_ML ) {
    Channel
    .fromPath(params.raw_phylip_ML)
    .flatten()
    .map{ item -> [item.baseName.replace('_random_msa_replicate.0_with_200_columns','') , item]}
    .set {rndMLFrac }
}

if ( params.raw_pairs_4mat ) {
    Channel
    .fromPath(params.raw_pairs_4mat)
    .flatten()
    .map{ item -> [ item.simpleName, item]}
    .set { rndPcol }
}


if ( params.originalmsa ) {
    Channel
    .fromPath(params.originalmsa)
    .map { item -> [ item.baseName.replace('_mTMalign','') , item] }
    .set { oriSeqs }
}


if ( params.templates ) {
    Channel
    .fromPath(params.templates)
    .map { item -> [ item.baseName, item] }
    .set { templates }
}


if ( params.pdb ) {
    Channel
    .fromPath(params.pdb)
    .map { item -> [ item.simpleName, item] }
    .groupTuple()
    .set { pdbs }
}


rndPcol
    .combine(oriSeqs, by:0)
    .combine(templates, by: 0)
    .combine(pdbs, by:0)
    .set{ ori_3d_mat_ch }


process computing_trees_on_1d_fractions{
    tag"${id}"
    publishDir "${params.output}/${id}/1d_trees/", mode: 'copy', overwrite: true
    label 'process_low'

    input:
        set val(id),file(fraphy) from rnd1dFrac

    output:
        set val(id), file("*.nwk"), file("*_fastme_boot.txt") into rnd1dFrT

    script:
    """
        /users/cn/lmansouri/bin/fastme -i ${fraphy} -m BioNJ -p LG -g ${params.gammaRate} -b ${params.replicatesNum} -s -n -z ${params.seedValue}
    """
}

process computing_ML_trees{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${id}/ML_trees/", mode: 'copy', overwrite: true
    label 'process_high'

    input:
        set val(id),file(phylip) from rndMLFrac

    output:
        set val(id), file("*.treefile"), file("*.boottrees") into iqtree_ch

    script:
    """
        /users/cn/abaltzis/bin/iqtree2 -s ${phylip} -b ${params.replicatesNum} 
    """
}


process computing_3d_matrices{
  // errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${id}/3d_matrices/", mode: 'copy', overwrite: true
    label 'process_medium'

    input:
        set val(id), file(pair), file(fasta), file(template), file(pdb) from ori_3d_mat_ch

    output:
        set val(id),file("*.matrices") into orimatOut

    script:
    """
        export THREED_TREE_MODE=${params.mode};
        t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pair} +replicates ${params.replicatesNum} +print_replicates +phylo3d  -output dm > ${pair}.matrices 2> ${pair}_3dtrees_err_log
    """
}



process extr_trmatrices_per_family{
    tag "${id}"
    publishDir "${params.output}/${id}/3d_matrices/", mode: 'copy', overwrite: true
    label 'process_low'

    input:
     set val(id), file(matrices) from orimatOut

    output:
     file("*.txt") into splitMatrix
     //set val(id), file ("*.1.*.txt")into input_tr_3dM_Ch

    script:
    """
    awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
    """
}

splitMatrix
    .flatten()
    .map { item -> [ item.simpleName, item] }
    .set{splitMatrix2}


process fastme_IMD_matrices{
    tag"${id}"
    publishDir "${params.output}/${id}/3d_trees/", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    label 'process_low'

    input:
    set val(id), file(mat) from splitMatrix2

    output:
    file("*.nwk") into matricesOut

    script:
    """
    fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}


workflow.onComplete {
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
