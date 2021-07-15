/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                    = [:]
analysis_scripts.run_split            = file("$baseDir/bin/seurat/seurat_split.R", checkIfExists: true)
// analysis_scripts.run_cluster          = file("$baseDir/bin/seurat/run_cluster.R", checkIfExists: true)
// analysis_scripts.run_gene_modules     = file("$baseDir/bin/other/run_gene_modules.R", checkIfExists: true)

params.run_split_options              = [:]
// params.run_cluster_options            = [:]
// params.run_gene_modules_options       = [:]

// Include Seurat R processes
include {R as RUN_SPLIT} from "$baseDir/modules/local/r/main"                 addParams(  options: params.run_split_options,
                                                                                            script: analysis_scripts.run_split )

// include {R as RUN_CLUSTER} from "$baseDir/modules/local/r/main"               addParams(  options: params.run_cluster_options,
//                                                                                             script: analysis_scripts.run_cluster )

// include {R as RUN_GENE_MODULES} from "$baseDir/modules/local/r/main"          addParams(  options: params.run_gene_modules_options,
//                                                                                             script: analysis_scripts.run_gene_modules )
/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_RUN_PROCESS {
    take:
    seurat_out //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // Run Seurat pipeline
    RUN_SPLIT( seurat_out )

    RUN_SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }

    // RUN_CLUSTER( ch_split_run )

    // RUN_GENE_MODULES( RUN_CLUSTER.out )

    // emit:
    // run_cluster_out   = RUN_CLUSTER.out //Channel: [[meta], [output]]
    // run_gene_modules_out = RUN_GENE_MODULES.out //Channel: [[meta], [output]]
}

