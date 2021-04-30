/*
 * Merge loom and run scVelo
 */

params.merge_loom_options         = [:]
params.seurat_intersect_loom_options = [:]

include { MERGE_LOOM } from '../../modules/local/scvelo/merge_loom/main' addParams( options: params.merge_loom_options )
include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom//main' addParams( options: params.seurat_intersect_loom_options ) 

workflow SCVELO {
    take:
    loom_dir
    seurat_input
    annotations

    main:
    /*
     * Merge loom files
     */
    // MERGE_LOOM ( loom_dir ) //channel: [ loom_directory ]
    SEURAT_INTERSECT_LOOM ( loom_dir, seurat_input, annotations ) //channel: [ loom_directory ]

    emit:
    // loom            = MERGE_LOOM.out.loom                  // channel: [ loom ]
    loom = SEURAT_INTERSECT_LOOM.out.loom
}