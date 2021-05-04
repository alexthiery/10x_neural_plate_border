/*
 * Merge loom and run scVelo
 */

params.merge_loom_options         = [:]
params.seurat_intersect_loom_options = [:]
params.scvelo_options = [:]

include { MERGE_LOOM } from '../../modules/local/scvelo/merge_loom/main' addParams( options: params.merge_loom_options )
include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom/main' addParams( options: params.seurat_intersect_loom_options ) 
include { SCVELO } from '../../modules/local/scvelo/scvelo/main' addParams( options: params.scvelo_options ) 

workflow SEURAT_SCVELO {
    take:
    loom_dir
    seurat_input
    annotations

    main:
    // println loom_dir

    MERGE_LOOM ( loom_dir ) //channel: [ loom_directory ]
    SEURAT_INTERSECT_LOOM ( MERGE_LOOM.out.loom, seurat_input, annotations ) //channel: [ loom_directory ]
    SCVELO (SEURAT_INTERSECT_LOOM.out.loom)
}