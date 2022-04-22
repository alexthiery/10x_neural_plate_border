/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts = [:]
analysis_scripts.preprocess     = file("$baseDir/bin/seurat/1_preprocess.R", checkIfExists: true)
analysis_scripts.integration        = file("$baseDir/bin/seurat/2_integration.R", checkIfExists: true)
analysis_scripts.integration_qc     = file("$baseDir/bin/seurat/3_integration_qc.R", checkIfExists: true)
analysis_scripts.sex_filt           = file("$baseDir/bin/seurat/4_sex_filt.R", checkIfExists: true)
analysis_scripts.cell_cycle         = file("$baseDir/bin/seurat/5_cell_cycle.R", checkIfExists: true)
analysis_scripts.contamination_filt = file("$baseDir/bin/seurat/6_contamination_filt.R", checkIfExists: true)

// Include Seurat R processes
include {R as PREPROCESS} from "$baseDir/modules/local/r/main"             addParams( script: analysis_scripts.preprocess )

include {R as INTEGRATION} from "$baseDir/modules/local/r/main"               addParams( script: analysis_scripts.integration )

include {R as INTEGRATION_QC} from "$baseDir/modules/local/r/main"            addParams( script: analysis_scripts.integration_qc )

include {R as SEX_FILT} from "$baseDir/modules/local/r/main"                  addParams( script: analysis_scripts.sex_filt )

include {R as CELL_CYCLE} from "$baseDir/modules/local/r/main"                addParams( script: analysis_scripts.cell_cycle )

include {R as CONTAMINATION_FILT} from "$baseDir/modules/local/r/main"        addParams( script: analysis_scripts.contamination_filt )


/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_FILTERING {
    take:
    cellranger_filtered_feature_bc_matrix

    main:
    // Run Seurat pipeline
    PREPROCESS( cellranger_filtered_feature_bc_matrix )
    INTEGRATION( PREPROCESS.out )
    INTEGRATION_QC( INTEGRATION.out )
    SEX_FILT( INTEGRATION_QC.out )
    CELL_CYCLE( SEX_FILT.out )
    CONTAMINATION_FILT( CELL_CYCLE.out )

    
    emit:
    preprocess_out              = PREPROCESS.out //Channel: [[meta], [output]]
    annotations                 = PREPROCESS.out.map{[it[0], it[1].findAll{it =~ /seurat_annotations.csv/}[0]]} //Channel: [[meta], annotations]
    integration_out             = INTEGRATION.out //Channel: [[meta], [output]]
    integration_qc_out          = INTEGRATION_QC.out //Channel: [[meta], [output]]
    sex_filt_out                = SEX_FILT.out //Channel: [[meta], [output]]
    cell_cycle_out              = CELL_CYCLE.out //Channel: [[meta], [output]]
    contamination_filt_out      = CONTAMINATION_FILT.out //Channel: [[meta], [output]]
}

