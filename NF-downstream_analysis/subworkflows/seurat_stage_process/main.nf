/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.split                              = file("$baseDir/bin/seurat/split_seurat.R", checkIfExists: true)
analysis_scripts.cluster                            = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.gene_modules                       = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification               = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)

params.split_options                                = [:]
params.cluster_options                              = [:]
params.gene_modules_options                         = [:]
params.state_classification_options                 = [:]

// Include Seurat R processes
include {R as SPLIT} from "$baseDir/modules/local/r/main"                           addParams(  options:                        params.split_options,
                                                                                                script:                         analysis_scripts.split )

include {R as CLUSTER} from "$baseDir/modules/local/r/main"                         addParams(  options:                        params.cluster_options,
                                                                                                script:                         analysis_scripts.cluster )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                    addParams(  options:                        params.gene_modules_options,
                                                                                                script:                         analysis_scripts.gene_modules )

include {R as STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"            addParams(  options:                        params.state_classification_options,
                                                                                                script:                         analysis_scripts.state_classification )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_STAGE_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // Run Seurat pipeline
    SPLIT( seurat_out )

    SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    CLUSTER( ch_split_run )
    STATE_CLASSIFICATION( CLUSTER.out )
    GENE_MODULES( STATE_CLASSIFICATION.out )

    emit:
    cluster_out                     = CLUSTER.out                               //Channel: [[meta], [output]]
    state_classification_out        = STATE_CLASSIFICATION.out                  //Channel: [[meta], [output]]
    gene_modules_out                = GENE_MODULES.out                          //Channel: [[meta], [output]]
}

