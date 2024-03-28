process comparing_splitfiles {
    tag "${id}_${type}"
    publishDir "${params.output}/split_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(imd_tree), path(mE_tree), path(mL_tree), path(ref_iMD), path(ref_mE), path(ref_mL)
    
    output:
    tuple val(id), path("*txt"), emit: splits

    script:
    """
    ./comparing_splits.py ${mE_tree} ${imd_tree} ${mL_tree} ${ref_iMD} ${ref_mE} ${ref_mL} ${id}_ME_splits_bs_3d_mL_mE_RefBR.txt
    """
}
