process MERGE_REPLICATES{
    tag"${id}"
    
    input:
    tuple val(id), path(reps) 
    
    output:
    tuple val(id), path("*.replicates"), emit: replicates
    
    script:
    """
    cat ${reps} | grep . - > ${id}.replicates
    """
}