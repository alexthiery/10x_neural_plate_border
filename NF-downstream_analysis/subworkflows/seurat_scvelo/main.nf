/*
run scVelo
 */

params.seurat_intersect_loom_options    = [:]
params.scvelo_options                   = [:]

include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom/main'  addParams( options: params.seurat_intersect_loom_options ) 
include { SCVELO } from '../../modules/local/scvelo/scvelo/main'                                addParams( options: params.scvelo_options ) 

workflow SEURAT_SCVELO {
    take:
    looms //channel: [loom]
    seurat_input //channel: [hd5a]
    annotations //channel: [annotations]

    main:

    SEURAT_INTERSECT_LOOM ( loom, seurat_input, annotations )
    SCVELO (SEURAT_INTERSECT_LOOM.out.loom)
}