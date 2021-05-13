/*
run scVelo
 */

params.seurat_intersect_loom_options    = [:]
params.scvelo_run_options                   = [:]

include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom/main'  addParams( options: params.seurat_intersect_loom_options ) 
include { SCVELO_RUN } from '../../modules/local/scvelo/scvelo_run/main'                                addParams( options: params.scvelo_run_options ) 

workflow SEURAT_SCVELO {
    take:
    loom //channel: [[meta], loom]
    seurat_input //channel: [[meta], hd5a]
    annotations //channel: [[meta], annotations]

    main:

    // Combine loom, seurat and annotations for running scvelo
    loom.concat(seurat_input, annotations)
        .groupTuple(by: 0, size: 3)
        .map{[it[0], it[1][0], it[1][1], it[1][2]]}
        .set{ch_scveloInput}

    SEURAT_INTERSECT_LOOM ( ch_scveloInput )
    SCVELO_RUN (SEURAT_INTERSECT_LOOM.out.loom)
}