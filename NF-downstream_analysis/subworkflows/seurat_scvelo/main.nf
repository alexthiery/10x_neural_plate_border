/*
 * Merge loom and run scVelo
 */

params.merge_loom_options               = [:]
params.seurat_intersect_loom_options    = [:]
params.scvelo_options                   = [:]

include { MERGE_LOOM } from '../../modules/local/scvelo/merge_loom/main'                        addParams( options: params.merge_loom_options )
include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom/main'  addParams( options: params.seurat_intersect_loom_options ) 
include { SCVELO } from '../../modules/local/scvelo/scvelo/main'                                addParams( options: params.scvelo_options ) 

workflow SEURAT_SCVELO {
    take:
    looms //channel: [loom]
    seurat_input //channel: [hd5a]
    annotations //channel: [annotations]

    main:
    // println loom_dir

    if(!params.merge_loom_options.skip_process){
        MERGE_LOOM ( looms ) //channel: [ looms ]
        loom = MERGE_LOOM.out
    } else {
        loom = looms
    }
    
    SEURAT_INTERSECT_LOOM ( loom, seurat_input, annotations )
    SCVELO (SEURAT_INTERSECT_LOOM.out.loom)
}