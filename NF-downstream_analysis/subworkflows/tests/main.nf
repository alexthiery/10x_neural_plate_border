#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

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
include {SCVELO} from "../scvelo/main" addParams( merge_loom_options: modules['integration'] )

workflow {
    SCVELO(ch_loomInput, ch_seuratInput, ch_annotations)
}