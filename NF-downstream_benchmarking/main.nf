#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"

include {r as integration_min} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_min'], script: modules['integration_min'].script)
include {r as integration_low} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_low'], script: modules['integration_low'].script)
include {r as integration_med} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_med'], script: modules['integration_med'].script)
include {r as integration_high} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_high'], script: modules['integration_high'].script)
include {r as integration_seurat_low} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat_low'], script: modules['integration_seurat_low'].script)
include {r as integration_seurat_med} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat_med'], script: modules['integration_seurat_med'].script)
include {r as integration_seurat_high} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat_high'], script: modules['integration_seurat_high'].script)

// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    metadata(params.input)

    // Run pipeline on seurat integration
    integration_min( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_low( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_med( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_high( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_seurat_low( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_seurat_med( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_seurat_high( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
}