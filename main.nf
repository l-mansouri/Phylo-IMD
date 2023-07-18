include { untrimmed_Phylo_IMD } from 'workflows/UNTRIMMED_PHYLO_IMD.nf '
include { trimmed_Phylo_IMD } from 'workflows/TRIMMES_PHYLO_IMD.nf'

workflow PHYLO_IMD{
    //Prepare input channels

    if ( params.fasta ) {
    Channel
        .fromPath(params.fasta)
        .map { item -> [ item.baseName.split('_')[0] , item] }
        .set { input_fasta }
    }

    if ( params.templates ) {
        Channel
            .fromPath(params.templates)
            .map { item -> [ item.baseName.split('_')[0] , item] }
            .set { templates }
    }

    if ( params.pdb ) {
        Channel
            .fromPath(params.pdb)
            .map { item -> [ item.simpleName , item] }
            .groupTuple()
            .set { structures }
    }

    //log info

    log.info """\
    PHYLO-IMD - version 1.0
    ====================================
    General parameters:
    -------------------
    Type of analysis                        : ${params.mode}
    Alignment method                        : ${params.align}
    Trimming                                : ${params.trimmer}

    Parameters for tree computation:
    --------------------------------
    Gamma rate                              : ${params.gammaRate}
    Random seed value                       : ${params.seedValue}
    Number of bootstrap replicates          : ${params.replicatesNum}
    IMD_tree mode                           : ${params.tree_mode}

    Input data:
    -----------
    Input (FASTA or input list)		        : ${params.fasta}
    Input template file		                : ${params.templates}
    Input structures (PDBs)                 : ${params.pdb}

    Output folder:                          : ${params.output}
    """

    //code with conditions
    if (params.mode == 'standard' && params.trimmer == 'trimmal'){
        trimmed_Phylo_IMD( input_fasta, templates, structures)
    } 
    else if (params.mode == 'standard' && params.trimmer != 'trimmal'){
        untrimmed_Phylo_IMD( input_fasta, templates, structures)
    }
    // else if (params.mode == 'titration'){

    // }
    // else if (params.mode == 'titration_bootstrap'){

    // }
}

workflow{
    PHYLO_IMD()
}


workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        workDir      : ${workflow.workDir}
        exit status  : ${workflow.exitStatus}
        runName      : ${workflow.runName}
        """
        .stripIndent()

}