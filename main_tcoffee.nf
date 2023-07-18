leilaDir = "/users/cn/lmansouri/PROJECTS/Phylo3D/NEW_Phylo3D/NF_draft"


params.mode = "10"   //1-11
params.replicatesNum="100"
params.alnr = 'tcoffee'


params.templates="${leilaDir}/template_lists/*.templ*"
params.fasta="${leilaDir}/fasta/PF*"
params.pdb="${leilaDir}/pdb/*.pdb"

params.gammaRate="1.0"
params.seedValue="5"

params.output = "${leilaDir}"

if ( params.fasta ) {
    Channel
        .fromPath(params.fasta)
        .map { item -> [ item.baseName , item] }
        .into { input_fasta_Ch ; input_fasta; in_fasta}
}

if ( params.templates ) {
    Channel
        .fromPath(params.templates)
        .map { item -> [ item.baseName.replace('_ref', '') , item] }
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

input_fasta_Ch
        .combine(templates, by:0)
        .combine(pdbFiles, by:0)
        .set{aln_ch}

process 'running_alignment' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_fasta", mode: 'copy', overwrite: true, pattern: "*.fasta"
    publishDir "${params.output}/${params.alnr}_ph", mode: 'copy', overwrite: true, pattern: "*.ph"

    input:
        set val(id), file(fasta) from in_fasta

    output:
        set val(id), file("*.fasta") into fasta_aln, fasta_aln2, fasta_aln3, original_msa_ch, oriSeqs, oriSeqs2
        set val(id), file("*.ph") into phylip_aln, phy_aln, phy_aln2
    
    script:
    """
        t_coffee -in ${fasta} -output=fasta_aln >${id}_${params.alnr}.fasta
        mv ${id}_${params.alnr}.fasta ${id}_${params.alnr}.clustal
        mv ${id}.fasta_aln ${id}_${params.alnr}.fasta
        t_coffee -other_pg seq_reformat -in ${fasta} -output phylip_aln > ${id}_${params.alnr}.ph

    """

}

process 'computing_untrimmed_ME_trees' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_trees/", mode: 'copy', overwrite: true, pattern: "*_untrimmed_ME_tree.nwk"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_trees/replicates", mode: 'copy', overwrite: true, pattern: "*.replicates"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_matrix/", mode: 'copy', overwrite: true, pattern: "*.matrix"


    input:
    set val(id),file(phylip) from phy_aln

    output:
    set val(id),file("*.nwk") into unttreeng1d
    set val(id), file("*.replicates") into untreplicateng1d
    set val(id), file ("*.mat*") into untrmat1d

    script:
    """
    fastme -i ${phylip} -o ${id}_${params.alnr}_untrimmed_ME_tree.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b 100 -B ${id}_${params.alnr}_untrimmed_ME_tree.replicates 
    fastme -i ${phylip} -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -O ${id}_${params.alnr}_untrimmed_ME.matrix
    """
}

process 'computing_untrimmed_ME_trees_NO_GAMMA' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_trees_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: "*_untrimmed_ME_tree_NO_GAMMA.nwk"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_trees_NO_GAMMA/replicates", mode: 'copy', overwrite: true, pattern: "*.replicates"
    publishDir "${params.output}/${params.alnr}_ME_untrimmed_matrix_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: "*.matrix"


    input:
    set val(id),file(phylip) from phylip_aln

    output:
    set val(id),file("*.nwk") into unttree1dng
    set val(id), file("*.replicates") into untreplicate1dng
    set val(id), file ("*.mat*") into untrmat1dng

    script:
    """
    fastme -i ${phylip} -o ${id}_${params.alnr}_untrimmed_ME_tree_NO_GAMMA.nwk  -m BioNJ -p LG -s -n -z ${params.seedValue} -b 100 -B ${id}_${params.alnr}_untrimmed._ME_tree_NO_GAMMA.replicates 
    fastme -i ${phylip} -m BioNJ -p LG -s -n -z ${params.seedValue} -O ${id}_${params.alnr}_untrimmed_ME_NO_GAMMA.matrix
    """
}

process 'computing_untrimmed_ML_trees'{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ML_untrimmed_trees/", mode: 'copy', overwrite: true, pattern: "*.treefile"
    publishDir "${params.output}/${params.alnr}_ML_untrimmed_trees/replicates", mode: 'copy', overwrite: true, pattern: "*.boottrees"
    
    input:
        set val(id), file(phylip) from phy_aln2

    output:
        set val(id), file("*.treefile") into unt_iqtree_ch
        set val(id), file("*.boottrees") into unt_iq_rep

    script:
    """
        iqtree -s ${phylip} -b ${params.replicatesNum} 
    """
}

/*process 'computing_untrimmed_NJ_tree'{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_NJ_untrimmed_trees/", mode: 'copy', overwrite: true

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
        /users/cn/lmansouri/tcoffee_extras/plugins/linux/neighbor < NJ_in
        mv outtree ${id}_${params.alnr}_untrimmed_NJ_tree.nwk
    """
}
*/
oriSeqs
        .combine(templates2, by:0)
        .combine(pdb_ch, by:0)
        .set{tc_mat_ch}


process 'computing_untrimmed_3d_matrices' {
    // errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_matrices/", mode: 'copy', overwrite: true

    input:
        set val(id), file(fasta), file(template), file(pdb) from tc_mat_ch


    output:
        set val(id),file("*.matrices") into untmatrixOut

    script:
    """
        export THREED_TREE_MODE=${params.mode}
        t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +replicates  ${params.replicatesNum} +phylo3d +print_replicates -output dm > ${id}_${params.alnr}_untrimmed_phylo3d.matrices 2>${id}_${params.alnr}_phylo3d_mat_err_log
    """
}

process 'extr_untrimmed_matrices_per_family' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_matrices/single_matrices", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_matrices/", mode: 'copy', overwrite: true, pattern: '*.1.txt'

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
    .into{untsplitMatrix2; untsplitMatrix2NG}


process 'fastme_untrimmed_3d_matrices' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_trees/replicates", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_trees/", mode: 'copy', overwrite: true, pattern: '*.1.*'
    errorStrategy 'ignore'

    input:
        set val(id), file(mat) from untsplitMatrix2

    output:
        file("*.nwk") into unttrees3d

    script:
    """
        fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}

process 'fastme_untrimmed_3d_matrices_NO_GAMMA' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_trees_NO_GAMMA/replicates", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_untrimmed_3d_trees_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: '*.1.*'
    errorStrategy 'ignore'

    input:
        set val(id), file(mat) from untsplitMatrix2NG

    output:
        file("*.nwk") into unttreesng3d

    script:
    """
        fastme -i ${mat} -s -n -z ${params.seedValue}
    """
}

process 'trimming_aln' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_TRIMMED_fasta/", mode: 'copy', overwrite: true, pattern: "*.fa"
    publishDir "${params.output}/${params.alnr}_TRIMMED_ph/", mode: 'copy', overwrite: true, pattern: "*.ph"

    input:
        set val(id), file(fasta) from fasta_aln

    output:
        set val(id), file("*.ph") into trimmed_aln, trimmed_aln2, blocked_aln
        set val(id), file("*.fa") into trimmed_fasta, blocked_fasta
    
    script:
    """
        trimal -in ${fasta} -out ${id}_${params.alnr}_trimmal_aln.ph -phylip -automated1
        trimal -in ${fasta} -out ${id}_${params.alnr}_trimmal_aln.fa -automated1
    """
}

process 'computing_trimmed_ME_trees'{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_trees/", mode: 'copy', overwrite: true, pattern: "*_trimmed_ME_tree.nwk"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_trees/replicates", mode: 'copy', overwrite: true, pattern: "*.replicates"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_matrix/", mode: 'copy', overwrite: true, pattern: "*.matrix"

    input:
        set val(id), file(phylip) from trimmed_aln

    output:
        set val(id), file("*.nwk"), file("*.replicates") into fastme_ch
        set val(id), file ("*.mat*") into trmat1d


    script:
    """
        fastme -i ${phylip} -o ${id}_${params.alnr}_trimmed_ME_tree.nwk  -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b ${params.replicatesNum} -B ${id}_${params.alnr}_trimmed_ME_tree.replicates 
        fastme -i ${phylip} -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -O ${id}_${params.alnr}_trimmed_ME.matrix
    """
}

process 'computing_trimmed_ME_trees_NO_GAMMA'{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_trees_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: "*_trimmed_ME_tree.nwk"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_trees_NO_GAMMA/replicates", mode: 'copy', overwrite: true, pattern: "*.replicates"
    publishDir "${params.output}/${params.alnr}_ME_trimmed_matrix_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: "*.matrix"

    input:
        set val(id), file(phylip) from trimmed_aln2


    output:
        set val(id), file("*.nwk"), file("*.replicates") into fastmeNG_ch
        set val(id), file ("*.mat*") into trmat1dng

    script:
    """
        fastme -i ${phylip} -o ${id}_${params.alnr}_trimmed_ME_tree.nwk  -m BioNJ -p LG -s -n -z ${params.seedValue} -b ${params.replicatesNum} -B ${id}_${params.alnr}_trimmed_ME_tree.replicates
        fastme -i ${phylip} -m BioNJ -p LG -s -n -z ${params.seedValue} -O ${id}_${params.alnr}_trimmed_ME_NO_GAMMA.matrix
    """
}

process 'computing_trimmed_ML_trees'{
    errorStrategy 'ignore'
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_ML_trimmed_trees/", mode: 'copy', overwrite: true, pattern: "*.treefile"
    publishDir "${params.output}/${params.alnr}_ML_trimmed_trees/replicates", mode: 'copy', overwrite: true, pattern: "*.boottrees"
    input:
        set val(id), file(phylip) from blocked_aln

    output:
        set val(id), file("*.treefile"), file("*.boottrees") into iqtree_ch

    script:
    """
        iqtree -s ${phylip} -b ${params.replicatesNum} 
    """
}

original_msa_ch
                .combine(trimmed_fasta, by:0)
                .set{for_mapping_ch}

process mapping_position {
    tag"${id}"
    input:
        set val(id), file(original_msa), file(trimmed_msa) from for_mapping_ch
    
    output:
        set val(id), file("*_selected_columns*.txt") into pairs_ch

    script:
    """
        python ${leilaDir}/mapping_position.py ${trimmed_msa} ${original_msa} ${id} ${params.trimmer}
    """
}

oriSeqs2
        .combine(templ_ch, by: 0)
        .combine(pairs_ch, by:0)
        .combine(pdbFiles2, by: 0)
        .set{phylo_3d_mat_ch}

process 'computing_trimmed_3d_matrices'{
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_trimmed_3d_matrices/", mode: 'copy', overwrite: true

    input:
        set val(id), file(fasta), file(template), file(pair), file(pdb) from phylo_3d_mat_ch

    output:
        set val(id),file("*.matrices") into orimatOut

    script:
    """
        export THREED_TREE_MODE=${params.mode}; 
        t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pair} +replicates ${params.replicatesNum} +print_replicates +phylo3d  -output dm > ${pair}.matrices 2> ${pair}_3dtrees_err_log
    """
}

process 'extr_trimmed_matrices_per_family'{
    tag "${id}"
    publishDir "${params.output}/${params.alnr}_trimmed_3d_matrices/single_matrices", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_trimmed_3d_matrices/", mode: 'copy', overwrite: true, pattern: '*.1.txt'

    input:
        set val(id), file(matrices) from orimatOut

    output:
        file("*.txt") into splitMatrix

    script:
    """
        awk -v RS= '{print > ("${matrices}."NR".txt")}' ${matrices}
    """
}

splitMatrix
    .flatten()
    .map { item -> [ item.simpleName, item] }
    .into{splitMatrix2; splitMatrix2NG}

process 'fastme_trimmed_3d_matrices' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_trimmed_3d_trees/replicates", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_trimmed_3d_trees/", mode: 'copy', overwrite: true, pattern: '*.1.*'
    errorStrategy 'ignore'

    input:
        set val(id), file(mat) from splitMatrix2

    output:
        file("*.nwk") into trees3d

    script:
    """
        fastme -i ${mat} -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}

process 'fastme_trimmed_3d_matrices_NO_GAMMA' {
    tag"${id}"
    publishDir "${params.output}/${params.alnr}_trimmed_3d_trees_NO_GAMMA/replicates", mode: 'copy', overwrite: true
    publishDir "${params.output}/${params.alnr}_trimmed_3d_trees_NO_GAMMA/", mode: 'copy', overwrite: true, pattern: '*.1.*'
    errorStrategy 'ignore'

    input:
        set val(id), file(mat) from splitMatrix2NG

    output:
        file("*.nwk") into treesng3d

    script:
    """
        fastme -i ${mat} -s -n -z ${params.seedValue}
    """
}

workflow.onComplete {
println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}