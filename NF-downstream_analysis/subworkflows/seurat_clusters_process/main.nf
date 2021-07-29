/*
 * Subset clusters from filtered Seurat object and run downstream analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.clusters_subset                    = file("$baseDir/bin/seurat/subset_cells.R", checkIfExists: true)
analysis_scripts.clusters_cluster                   = file("$baseDir/bin/seurat/cluster_cells.R", checkIfExists: true)
analysis_scripts.clusters_gene_modules              = file("$baseDir/bin/other/clusters_gene_modules.R", checkIfExists: true)
analysis_scripts.clusters_state_classification      = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)

params.clusters_subset_options                      = [:]
params.clusters_cluster_options                     = [:]
params.clusters_gene_modules_options                = [:]
params.clusters_state_classification_options        = [:]
params.seurat_h5ad_options                          = [:]
params.seurat_intersect_loom_options                = [:]
params.scvelo_run_options                           = [:]



// Include Seurat R processes
include {R as CLUSTERS_SUBSET} from "$baseDir/modules/local/r/main"                 addParams(  options:                        params.clusters_subset_options,
                                                                                                script:                         analysis_scripts.clusters_subset )

include {R as CLUSTERS_CLUSTER} from "$baseDir/modules/local/r/main"                addParams(  options:                        params.clusters_cluster_options,
                                                                                                script:                         analysis_scripts.clusters_cluster )

include {R as CLUSTERS_GENE_MODULES} from "$baseDir/modules/local/r/main"           addParams(  options:                        params.clusters_gene_modules_options,
                                                                                                script:                         analysis_scripts.clusters_gene_modules )

include {R as CLUSTERS_STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"   addParams(  options:                        params.clusters_state_classification_options,
                                                                                                script:                         analysis_scripts.clusters_state_classification )

include {SEURAT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"                addParams(  options:                        params.seurat_h5ad_options )

include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"             addParams(  seurat_intersect_loom_options:  params.seurat_intersect_loom_options,
                                                                                                scvelo_run_options:             params.scvelo_run_options)

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_CLUSTERS_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]
    loom            //Channel: merged.loom
    annotations     //Channel: seurat_annotations.csv

    main:
    // Run Seurat pipeline
    CLUSTERS_SUBSET( seurat_out )
    CLUSTERS_SUBSET.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_clusters_subset }

    CLUSTERS_CLUSTER( ch_clusters_subset )
    CLUSTERS_GENE_MODULES( CLUSTERS_CLUSTER.out )
    CLUSTERS_STATE_CLASSIFICATION( CLUSTERS_CLUSTER.out )

    SEURAT_H5AD( CLUSTERS_CLUSTER.out )

    // Run scVelo
    SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv

    emit:
    clusters_cluster_out                = CLUSTERS_CLUSTER.out              //Channel: [[meta], [output]]
    clusters_gene_modules_out           = CLUSTERS_GENE_MODULES.out         //Channel: [[meta], [output]]
    clusters_state_classification_out   = CLUSTERS_STATE_CLASSIFICATION.out //Channel: [[meta], [output]]
    clusters_scvelo_out                 = SEURAT_SCVELO.out                 //Channel:
}

