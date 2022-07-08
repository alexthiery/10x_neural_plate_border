/*
 * Generate HCR intensity plots and model expression in scRNAseq
 */


def analysis_scripts = [:]
analysis_scripts.hcr_intensity_plot             = file("$baseDir/bin/other/hcr_intensity_plot.R", checkIfExists: true)
analysis_scripts.spatial_expression_modelling   = file("$baseDir/bin/other/spatial_expression_modelling.R", checkIfExists: true)

// Include subworkflows
include {METADATA} from "$baseDir/subworkflows/metadata/main"

// Include modules
include {R as HCR_INTENSITY_PLOT} from "$baseDir/modules/local/r/main"                      addParams( script: analysis_scripts.hcr_intensity_plot )
include {R as SPATIAL_EXPRESSION_MODELLING} from "$baseDir/modules/local/r/main"            addParams( script: analysis_scripts.spatial_expression_modelling )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow HCR {
    take:
    hcr_intensity_samplesheet
    hh7_seurat

    main:
    METADATA( hcr_intensity_samplesheet )
    HCR_INTENSITY_PLOT( METADATA.out )
    SPATIAL_EXPRESSION_MODELLING( hh7_seurat )
}

