/*
 * Downstream exploratory analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                            = [:]
analysis_scripts.scatterplot3d                  = file("$baseDir/bin/other/scatterplot3d.R", checkIfExists: true)
analysis_scripts.gene_modules                   = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.cell_state_classification      = file("$baseDir/bin/seurat/cell_state_classification.R", checkIfExists: true)
analysis_scripts.subset_npb                   = file("$baseDir/bin/seurat/subset_cells.R", checkIfExists: true)

params.scatterplot3d_options                    = [:]
params.gene_module_options                      = [:]
params.cell_state_classification_options        = [:]
params.subset_npb_options                    = [:]

// Include R processes

include {R as SCATTERPLOT3D} from "$baseDir/modules/local/r/main"                addParams(     options: params.scatterplot3d_options,
                                                                                                script: analysis_scripts.scatterplot3d )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                addParams(      options: params.gene_module_options,
                                                                                                script: analysis_scripts.gene_modules )

include {R as CELL_STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"   addParams(      options: params.cell_state_classification_options,
                                                                                                script: analysis_scripts.cell_state_classification )

include {R as SUBSET_NPB} from "$baseDir/modules/local/r/main"                  addParams(      options: params.subset_npb_options,
                                                                                                script: analysis_scripts.subset_npb )


/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(params, analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow EXPLORATORY_ANALYSIS {
    take:
    seurat_out //Channel: [[meta], [output]]

    main:
    // Run Seurat pipeline
    SCATTERPLOT3D(seurat_out)
    GENE_MODULES( seurat_out )

    CELL_STATE_CLASSIFICATION( seurat_out )
    SUBSET_NPB( CELL_STATE_CLASSIFICATION.out )
    
    emit:
    gene_modules_out = GENE_MODULES.out //Channel: [[meta], [output]]
    cell_state_classification_out = CELL_STATE_CLASSIFICATION.out //Channel: [[meta], [output]]
}

