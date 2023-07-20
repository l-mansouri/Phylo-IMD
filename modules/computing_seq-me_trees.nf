process 'computing_Seq-ME_trees' {

    tag"${id}"
    publishDir "${params.output}/ME_trees/", mode: 'copy', overwrite: true

    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq-ME
    tuple val(id), path("*.replicates"), emit: rep_Seq-ME

    script:
    """
    fastme -i ${phylip} -o ${id}_${params.align}_${params.trimmer}.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b 100 -B ${id}_${params.align}_${params.trimmer}.1dtree.replicates
    """
}

process 'computing_Seq-ME_trees_no_bs' {

    tag"${id}"
    publishDir "${params.output}/ME_trees/", mode: 'copy', overwrite: true

    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq-ME
    tuple val(id), path("*.replicates"), emit: rep_Seq-ME

    script:
    """
    fastme -i ${phylip} -o ${id}_${params.align}_${params.trimmer}.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue}
    """
}