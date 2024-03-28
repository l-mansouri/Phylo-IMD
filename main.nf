include { untrimmed_Phylo_IMD           } from './workflows/UNTRIMMED_PHYLO_IMD.nf'
include { trimmed_Phylo_IMD             } from './workflows/TRIMMED_PHYLO_IMD.nf'
include { titration_Phylo_IMD           } from './workflows/TITRATION.nf'
include { titration_bootstrap_Phylo_IMD } from './workflows/TITRATION_BOOTSTRAP.nf'
include { MULTISTRAP                    } from './workflows/MULTISTRAP.nf'
include { ANALYSIS                      } from './workflows/ANALYSIS.nf'

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
    if (params.mode == 'standard' && params.trimmer == 'trimmed'){
        trimmed_Phylo_IMD( input_fasta, templates, structures )
    } 
    else if (params.mode == 'standard' && params.trimmer == 'untrimmed'){
        untrimmed_Phylo_IMD( input_fasta, templates, structures )
    }
    else if (params.mode == 'titration'){
        titration_Phylo_IMD( input_fasta, templates, structures )
    }
    else if (params.mode == 'titration_bootstrap'){
        titration_bootstrap_Phylo_IMD( input_fasta, templates, structures )
    }
    else if (params.mode == "multistrap"){
        MULTISTRAP( input_fasta, templates, structures )
    }
    else if (params.mode == "analysis"){

        if ( params.msas ) {
            Channel
                .fromPath(params.msas)
                .map { item -> [ item.simpleName.split("_")[0], item.simpleName.split("_")[1..-1].join("_").replaceAll(".fa", ""), item] }
                .filter { id, method, file -> method in ["tcoffee", "sap_tmalign", "mTMalign"] }
                .set { msas }
        }

        if ( params.trees ) {
            Channel
                .fromPath(trees)
                .map { item -> [ item.baseName.split('_')[0], item] }
                .set{trees_ch}
        }
        

        if ( params.replicates ) {
            Channel
                .fromPath(replicates)
                .map { item -> [ item.baseName.split('_')[0], item] }
                .set{replicates_ch}
        }


        ANALYSIS( msas, templates, structures )
    }
    else{
        error "Mode not recognized"
    }

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