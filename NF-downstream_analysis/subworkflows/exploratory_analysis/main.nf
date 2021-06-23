/*
 * Downstream exploratory analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts              = [:]
analysis_scripts.gene_modules     = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)

params.gene_module_options       = [:]

// Include R processes
include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"          addParams(        options: params.gene_module_options,
                                                                                            script: analysis_scripts.gene_modules )
/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow EXPLORATORY_ANALYSIS {
    take:
    seurat_out //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // Run Seurat pipeline
    GENE_MODULES( seurat_out )
    
    emit:
    gene_modules_out = GENE_MODULES.out //Channel: [[meta], [output]]
}

