/*
 * Split Seurat object by batch and run downstream analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.run_split                          = file("$baseDir/bin/seurat/split_seurat.R", checkIfExists: true)
analysis_scripts.run_cluster                        = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.run_gene_modules                   = file("$baseDir/bin/other/subset_gene_modules.R", checkIfExists: true)
analysis_scripts.run_state_classification           = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)

params.run_split_options                            = [:]
params.run_cluster_options                          = [:]
params.run_gene_modules_options                     = [:]
params.run_state_classification_options             = [:]
params.seurat_h5ad_options                          = [:]
params.seurat_intersect_loom_options                = [:]
params.scvelo_run_options                           = [:]



// Include Seurat R processes
include {R as RUN_SPLIT} from "$baseDir/modules/local/r/main"                       addParams(  options: params.run_split_options,
                                                                                                script: analysis_scripts.run_split )

include {R as RUN_CLUSTER} from "$baseDir/modules/local/r/main"                     addParams(  options: params.run_cluster_options,
                                                                                                script: analysis_scripts.run_cluster )

include {R as RUN_GENE_MODULES} from "$baseDir/modules/local/r/main"                addParams(  options: params.run_gene_modules_options,
                                                                                                script: analysis_scripts.run_gene_modules )

include {R as RUN_STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"        addParams(  options: params.run_state_classification_options,
                                                                                                script: analysis_scripts.run_state_classification )

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

workflow SEURAT_RUN_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]
    loom            //Channel: merged.loom
    annotations     //Channel: seurat_annotations.csv

    main:
    // Run Seurat pipeline
    RUN_SPLIT( seurat_out )

    RUN_SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    RUN_CLUSTER( ch_split_run )
    RUN_GENE_MODULES( RUN_CLUSTER.out )
    RUN_STATE_CLASSIFICATION( RUN_CLUSTER.out )

    SEURAT_H5AD( RUN_CLUSTER.out )

    // Run scVelo
    SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv

    emit:
    run_cluster_out                     = RUN_CLUSTER.out                               //Channel: [[meta], [output]]
    run_gene_modules_out                = RUN_GENE_MODULES.out                          //Channel: [[meta], [output]]
    run_state_classification_out        = RUN_STATE_CLASSIFICATION.out                  //Channel: [[meta], [output]]   

    run_scvelo_run_out_metadata         = SEURAT_SCVELO.out.scvelo_run_out_metadata     //Channel: [[meta], csv]
    run_scvelo_run_out_5had             = SEURAT_SCVELO.out .scvelo_run_out_h5ad        //Channel: [[meta], h5ad]
}
