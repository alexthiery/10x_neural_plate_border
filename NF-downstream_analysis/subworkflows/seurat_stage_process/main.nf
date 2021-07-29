/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.stage_split                        = file("$baseDir/bin/seurat/split_seurat.R", checkIfExists: true)
analysis_scripts.stage_cluster                      = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.stage_gene_modules                 = file("$baseDir/bin/other/subset_gene_modules.R", checkIfExists: true)
analysis_scripts.stage_state_classification         = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)

params.stage_split_options                          = [:]
params.stage_cluster_options                        = [:]
params.stage_gene_modules_options                   = [:]
params.stage_state_classification_options           = [:]
params.seurat_h5ad_options                          = [:]
params.seurat_intersect_loom_options                = [:]
params.scvelo_run_options                           = [:]



// Include Seurat R processes
include {R as STAGE_SPLIT} from "$baseDir/modules/local/r/main"                     addParams(  options:                        params.stage_split_options,
                                                                                                script:                         analysis_scripts.stage_split )

include {R as STAGE_CLUSTER} from "$baseDir/modules/local/r/main"                   addParams(  options:                        params.stage_cluster_options,
                                                                                                script:                         analysis_scripts.stage_cluster )

include {R as STAGE_GENE_MODULES} from "$baseDir/modules/local/r/main"              addParams(  options:                        params.stage_gene_modules_options,
                                                                                                script:                         analysis_scripts.stage_gene_modules )

include {R as STAGE_STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"      addParams(  options:                        params.stage_state_classification_options,
                                                                                                script:                         analysis_scripts.stage_state_classification )

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

workflow SEURAT_STAGE_PROCESS {
    take:
    seurat_out //Channel: [[meta], [plot_dir, rds_dir]]
    loom            //Channel: merged.loom
    annotations     //Channel: seurat_annotations.csv

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
    STAGE_STATE_CLASSIFICATION( STAGE_CLUSTER.out )

    SEURAT_H5AD( STAGE_CLUSTER.out )

    // Run scVelo
    SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv

    emit:
    stage_cluster_out                   = STAGE_CLUSTER.out                             //Channel: [[meta], [output]]
    stage_gene_modules_out              = STAGE_GENE_MODULES.out                        //Channel: [[meta], [output]]
    stage_state_classification_out      = STAGE_STATE_CLASSIFICATION.out                //Channel: [[meta], [output]]

    stage_scvelo_run_out_metadata       = SEURAT_SCVELO.out.scvelo_run_out_metadata     //Channel: [[meta], csv]
    stage_scvelo_run_out_5had           = SEURAT_SCVELO.out .scvelo_run_out_h5ad        //Channel: [[meta], h5ad]
}

