process 'computing_Seq_ME_trees' {

    tag"${id}"
    container 'lmansouri/phylo_imd_fastme:1.0'

    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq_ME
    tuple val(id), path("*.replicates"), emit: rep_Seq_ME

    script:
    """
    fastme -i ${phylip} -o ${phylip.baseName}.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -b ${params.replicatesNum} -B ${phylip.baseName}.replicates
    if [ ! -f *_columns.nwk ]; then
        mv *.nwk ${id}_${params.align}_${params.trimmer}_ME.nwk
        mv *.replicates ${id}_${params.align}_${params.trimmer}_ME.replicates
    fi
    """
}

process 'computing_Seq_ME_trees_no_bs' {

    tag"${id}"
    container 'lmansouri/phylo_imd_fastme:1.0'


    input:
    tuple val(id),path(phylip) 

    output:
    tuple val(id),path("*.nwk"), emit: tr_Seq_ME
    tuple val(id), path("*.matrix"), emit: mat_Seq_ME


    script:
    """
    fastme -i ${phylip} -o ${phylip.baseName}.nwk -m BioNJ -p LG -g ${params.gammaRate} -s -n -z ${params.seedValue} -O ${phylip.baseName}.matrix
    if [ ! -f *_columns.nwk ]; then
        mv *.nwk ${id}_${params.align}_${params.trimmer}_ME.nwk
        mv *.matrix ${id}_${params.align}_${params.trimmer}_ME.matrix
    fi
    """
}