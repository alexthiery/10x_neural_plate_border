/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

if(params.integration && params.integration != 'Seurat' && params.integration != 'STACAS'){
    log.error "'${params.integration}' is not a valid option for integration. Please use 'Seurat' or 'STACAS'."
    System.exit(1)
}

def analysis_scripts = [:]
analysis_scripts.preprocessing     = file("$baseDir/bin/seurat/1_preprocessing.R", checkIfExists: true)
analysis_scripts.integration        = params.integration == 'Seurat' ? file("$baseDir/bin/seurat/2_integration.R", checkIfExists: true) : file("$baseDir/bin/seurat/2_integration_STACAS.R", checkIfExists: true)
analysis_scripts.integration_qc     = file("$baseDir/bin/seurat/3_integration_qc.R", checkIfExists: true)
analysis_scripts.sex_filt           = file("$baseDir/bin/seurat/4_sex_filt.R", checkIfExists: true)
analysis_scripts.cell_cycle         = file("$baseDir/bin/seurat/5_cell_cycle.R", checkIfExists: true)
analysis_scripts.contamination_filt = file("$baseDir/bin/seurat/6_contamination_filt.R", checkIfExists: true)

params.preprocessing_options        = [:]
params.integration_options          = [:]
params.integration_qc_options       = [:]
params.sex_filt_options             = [:]
params.cell_cycle_options           = [:]
params.contamination_filt_options   = [:]

// Include Seurat R processes
include {R as PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(      options: params.preprocessing_options,
                                                                                                script: analysis_scripts.preprocessing )

include {R as INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(        options: params.integration_options,
                                                                                                script: analysis_scripts.integration )

include {R as INTEGRATION_QC} from "$baseDir/modules/local/r/main"            addParams(        options: params.integration_qc_options,
                                                                                                script: analysis_scripts.integration_qc )

include {R as SEX_FILT} from "$baseDir/modules/local/r/main"                  addParams(        options: params.sex_filt_options,
                                                                                                script: analysis_scripts.sex_filt )

include {R as CELL_CYCLE} from "$baseDir/modules/local/r/main"                addParams(        options: params.cell_cycle_options,
                                                                                                script: analysis_scripts.cell_cycle )

include {R as CONTAMINATION_FILT} from "$baseDir/modules/local/r/main"        addParams(        options: params.contamination_filt_options,
                                                                                                script: analysis_scripts.contamination_filt )


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
    PREPROCESSING( cellranger_filtered_feature_bc_matrix )
    INTEGRATION( PREPROCESSING.out )
    INTEGRATION_QC( INTEGRATION.out )
    SEX_FILT( INTEGRATION_QC.out )
    CELL_CYCLE( SEX_FILT.out )
    CONTAMINATION_FILT( CELL_CYCLE.out )

    
    emit:
    preprocessing_out           = PREPROCESSING.out //Channel: [[meta], [output]]
    annotations                 = PREPROCESSING.out.map{[it[0], it[1].findAll{it =~ /seurat_annotations.csv/}[0]]} //Channel: [[meta], annotations]
    integration_out             = INTEGRATION.out //Channel: [[meta], [output]]
    integration_qc_out          = INTEGRATION_QC.out //Channel: [[meta], [output]]
    sex_filt_out                = SEX_FILT.out //Channel: [[meta], [output]]
    cell_cycle_out              = CELL_CYCLE.out //Channel: [[meta], [output]]
    contamination_filt_out      = CONTAMINATION_FILT.out //Channel: [[meta], [output]]
}

