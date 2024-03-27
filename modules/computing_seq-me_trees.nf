process 'computing_Seq_ME_trees' {

    tag"${id}"
    publishDir "${params.output}/${params.align}_ME_${params.trimmer}_trees/", mode: 'copy', overwrite: true, pattern: "*.nwk"
    publishDir "${params.output}/${params.align}_ME_${params.trimmer}_trees/replicates", mode: 'copy', overwrite: true, pattern: "*.replicates"
    container 'lmansouri/phylo_imd_fastme:1.0'

    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq_ME
    tuple val(id), path("*.replicates"), emit: rep_Seq_ME

    script:
    """
    fastme -i ${phylip} -o ${id}_${params.align}_${params.trimmer}_ME.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b ${params.replicatesNum} -B ${id}_${params.align}_${params.trimmer}_ME.replicates
    """
}

process 'computing_Seq_ME_trees_no_bs' {

    tag"${id}"
    publishDir "${params.output}/${params.align}_ME_${params.trimmer}_matrix/", mode: 'copy', overwrite: true, pattern: "*.matrix"
    container 'lmansouri/phylo_imd_fastme:1.0'


    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq_ME
    tuple val(id), path("*.matrix"), emit: mat_Seq_ME


    script:
    """
    fastme -i ${phylip} -o ${id}_${params.align}_${params.trimmer}_ME.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -O ${id}_${params.align}_${params.trimmer}_ME.matrix
    """
}