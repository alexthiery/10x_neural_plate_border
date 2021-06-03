/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                    = [:]
analysis_scripts.stage_split            = file("$baseDir/bin/seurat/stage_split.R", checkIfExists: true)
analysis_scripts.stage_cluster          = file("$baseDir/bin/seurat/stage_cluster.R", checkIfExists: true)
analysis_scripts.stage_gene_modules     = file("$baseDir/bin/other/stage_gene_modules.R", checkIfExists: true)

params.stage_split_options              = [:]
params.stage_cluster_options            = [:]
params.stage_gene_modules_options       = [:]

// Include Seurat R processes
include {R as STAGE_SPLIT} from "$baseDir/modules/local/r/main"                 addParams(  options: params.stage_split_options,
                                                                                            script: analysis_scripts.stage_split )

include {R as STAGE_CLUSTER} from "$baseDir/modules/local/r/main"               addParams(  options: params.stage_cluster_options,
                                                                                            script: analysis_scripts.stage_cluster )

include {R as STAGE_GENE_MODULES} from "$baseDir/modules/local/r/main"          addParams(  options: params.stage_gene_modules_options,
                                                                                            script: analysis_scripts.stage_gene_modules )
/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_STAGE_PROCESS {
    take:
    seurat_out //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // Run Seurat pipeline
    STAGE_SPLIT( seurat_out )

    STAGE_SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_stage }

    STAGE_CLUSTER( ch_split_stage )

    STAGE_GENE_MODULES( STAGE_CLUSTER.out )

    emit:
    stage_cluster_out   = STAGE_CLUSTER.out //Channel: [[meta], [output]]
}

