/*
 * Downstream exploratory analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts              = [:]
analysis_scripts.gene_modules_latent_time     = file("$baseDir/bin/other/gene_modules_latent_time.R", checkIfExists: true)

params.gene_modules_latent_time_options       = [:]

// Include R processes
include {R as GENE_MODULES_LATENT_TIME} from "$baseDir/modules/local/r/main"                  addParams(    options: params.gene_modules_latent_time_options,
                                                                                                            script: analysis_scripts.gene_modules_latent_time )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow GENE_MODULE_LATENT_TIME {
    take:
    seurat_out //Channel: [[meta], *.rds_file]
    antler_out //Channel: [[meta], *.rds_file]
    scvelo_out //Channel: [[meta], *.csv]

    main:
    
    ch_gene_module_latent_time = seurat_out.combine(antler_out, by: 0).combine(scvelo_out, by: 0)
    ch_gene_module_latent_time = ch_gene_module_latent_time.map{[it[0], [it[1], it[2], it[3]]]}

    // Run gene modules across latent time
    GENE_MODULES_LATENT_TIME( ch_gene_module_latent_time )
}

