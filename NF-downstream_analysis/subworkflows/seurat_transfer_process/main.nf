/*
 * Subset clusters from filtered Seurat object and run downstream analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.subset                             = file("$baseDir/bin/seurat/subset_cells.R", checkIfExists: true)
analysis_scripts.cluster                            = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.gene_modules                       = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification               = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)
analysis_scripts.phate                              = file("$baseDir/bin/other/phateR.R", checkIfExists: true)

params.subset_options                               = [:]
params.cluster_options                              = [:]
params.gene_modules_options                         = [:]
params.state_classification_options                 = [:]
params.phate_options                                = [:]

// Include Seurat R processes
include {R as SUBSET} from "$baseDir/modules/local/r/main"                          addParams(  options:                        params.subset_options,
                                                                                                script:                         analysis_scripts.subset )

include {R as CLUSTER} from "$baseDir/modules/local/r/main"                         addParams(  options:                        params.cluster_options,
                                                                                                script:                         analysis_scripts.cluster )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                    addParams(  options:                        params.gene_modules_options,
                                                                                                script:                         analysis_scripts.gene_modules )
/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_TRANSFER_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // Run Seurat pipeline
    SUBSET( seurat_out )

    SUBSET.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    CLUSTER( ch_split_run )
    GENE_MODULES( CLUSTER.out )

    emit:
    cluster_out                     = CLUSTER.out                               //Channel: [[meta], [output]]
    gene_modules_out                = GENE_MODULES.out                          //Channel: [[meta], [output]]
}

