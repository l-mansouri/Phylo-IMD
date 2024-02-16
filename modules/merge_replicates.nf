process MERGE_REPLICATES{
    tag"${id}"
    
    input:
    tuple val(id), path(reps) 
    
    output:
    tuple val(id), path("*.nwk"), path("*.replicates"), emit: replicates
    
    script:
    """
    cat ${reps} | grep . - > tmp
    # extract the tree
    head -n 1 tmp > ${id}.nwk
    tail -n +2 tmp > ${id}.replicates
    """
}