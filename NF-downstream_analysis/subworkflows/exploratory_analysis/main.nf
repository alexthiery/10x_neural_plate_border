/*
 * Downstream exploratory analysis on full Seurat filtered dataset
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                            = [:]
analysis_scripts.scatterplot3d                  = file("$baseDir/bin/other/scatterplot3d.R", checkIfExists: true)
analysis_scripts.gene_modules                   = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification           = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)

params.scatterplot3d_options                    = [:]
params.gene_module_options                      = [:]
params.state_classification_options        = [:]

// Include R processes

include {R as SCATTERPLOT3D} from "$baseDir/modules/local/r/main"                addParams(     options: params.scatterplot3d_options,
                                                                                                script: analysis_scripts.scatterplot3d )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                addParams(      options: params.gene_module_options,
                                                                                                script: analysis_scripts.gene_modules )

include {R as STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"        addParams(      options: params.state_classification_options,
                                                                                                script: analysis_scripts.state_classification )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow EXPLORATORY_ANALYSIS {
    take:
    seurat_out //Channel: [[meta], [output]]

    main:
    // Run processes on full filtered dataset
    SCATTERPLOT3D(seurat_out)
    GENE_MODULES( seurat_out )
    STATE_CLASSIFICATION( seurat_out )
    
    emit:
    gene_modules_out = GENE_MODULES.out //Channel: [[meta], [output]]
    state_classification_out = STATE_CLASSIFICATION.out //Channel: [[meta], [output]]
}

