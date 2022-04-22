/*
run scVelo
 */

include { SEURAT_INTERSECT_LOOM } from '../../modules/local/scvelo/seurat_intersect_loom/main'
include { SCVELO_RUN } from '../../modules/local/scvelo/scvelo_run/main'

workflow SEURAT_SCVELO {
    take:
    seurat_input //channel: [[meta], hd5a]
    loom //channel: loom
    annotations //channel: annotations

    main:

    // Combine loom, seurat and annotations for running scvelo
    seurat_input
        .combine(loom)
        .combine(annotations)
        .set{ch_scveloInput}
        
    SEURAT_INTERSECT_LOOM ( ch_scveloInput ) //channel: [[meta], seurat, loom, annotations]
    SCVELO_RUN (SEURAT_INTERSECT_LOOM.out.loom) //channel: [[meta], loom]

    emit:
    // scvelo_run_out_metadata     = SCVELO_RUN.out.csv //Channel: [[meta], csv]
    scvelo_run_out_h5ad         = SCVELO_RUN.out.h5ad //Channel: [[meta], h5ad]
}