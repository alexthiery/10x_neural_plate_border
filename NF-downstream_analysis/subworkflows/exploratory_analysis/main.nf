/*
 * Downstream exploratory analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                            = [:]
analysis_scripts.gene_modules                   = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.cell_state_classification      = file("$baseDir/bin/seurat/cell_state_classification.R", checkIfExists: true)
analysis_scripts.scatterplot3                   = file("$baseDir/bin/seurat/scatterplot3d.R", checkIfExists: true)

params.gene_module_options                      = [:]
params.cell_state_classification_options        = [:]
params.scatterplot3d_options                    = [:]

// Include R processes
include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                addParams(      options: params.gene_module_options,
                                                                                                script: analysis_scripts.gene_modules )

include {R as CELL_STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"   addParams(      options: params.cell_state_classification_options,
                                                                                                script: analysis_scripts.cell_state_classification )

include {R as SCATTERPLOT3D} from "$baseDir/modules/local/r/main"               addParams(      options: params.scatterplot3d_options,
                                                                                                script: analysis_scripts.scatterplot3d)

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
    // Run Seurat pipeline
    GENE_MODULES( seurat_out )
    CELL_STATE_CLASSIFICATION( seurat_out )
    SCATTERPLOT3D(CELL_STATE_CLASSIFICATION.out)
    
    emit:
    gene_modules_out = GENE_MODULES.out //Channel: [[meta], [output]]
    cell_state_classification_out = CELL_STATE_CLASSIFICATION.out //Channel: [[meta], [output]]
}

