/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts = [:]
analysis_scripts.stage_split        = file("$baseDir/bin/seurat/stage_split.R", checkIfExists: true)

params.stage_split_options          = [:]

// Include Seurat R processes
include {R as STAGE_SPLIT} from "$baseDir/modules/local/r/main"               addParams(    options: params.stage_split_options,
                                                                                            script: analysis_scripts.stage_split )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_STAGE_PROCESS {
    take:
    seurat_out //Channel: [[meta], rds_dir_path]

    main:
    // Run Seurat pipeline
    SEURAT_STAGE_SPLIT( seurat_out )

    SEURAT_STAGE_SPLIT.out
        .flatMap {it[1].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_stage }
    
    emit:
    test                 = ch_split_stage //Channel: [[meta], annotations]
}

