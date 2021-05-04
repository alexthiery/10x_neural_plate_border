#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// This is to conditionally skip the merging process IF loomInput is not a directory... Ask chris the best way to do this....
file(params.loomInput).isFile()? :


Channel
    .fromPath( params.loomInput )
    .set { ch_loomInput }

Channel
    .fromPath( params.seuratInput )
    .set { ch_seuratInput }

Channel
    .fromPath( params.annotations )
    .set { ch_annotations }

/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/
include {SEURAT_SCVELO} from "../scvelo/main" addParams(    merge_loom_options:             modules['merge_loom'],
                                                            seurat_intersect_loom_options:  modules['seurat_intersect_loom'],
                                                            scvelo_options:                 modules['scvelo'] )

workflow {

    SEURAT_SCVELO(ch_loomInput, ch_seuratInput, ch_annotations)
}