process 'mapping_pos' {
    tag"${id}"
    container 'lmansouri/phylo_imd_iqtree:1.0'
    
    input:
        tuple val(id), path(original_msa), path(trimmed_msa) 
    
    output:
        tuple val(id), path("*_selected_columns*.txt"), emit: mapped_pairs

    script:
    """
        python ${baseDir}/bin/mapping_position.py ${trimmed_msa} ${original_msa} ${id} ${params.trimmer}
    """
}