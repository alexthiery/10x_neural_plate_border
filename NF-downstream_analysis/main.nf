#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"

include {r as integration} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration'], script: modules['integration'].script)
include {r as integration_qc} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_qc'], script: modules['integration_qc'].script)
include {r as sexfilt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['sexfilt'], script: modules['sexfilt'].script)

// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    metadata(params.input)

    // Run pipeline on seurat integration
    integration( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_qc( integration.out )
    sexfilt( integration_qc.out )

}