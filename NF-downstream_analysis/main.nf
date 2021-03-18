#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"

include {r as integration_seurat} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat'], script: modules['integration_seurat'].script)
include {r as integration_seurat_qc} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat_qc'], script: modules['integration_seurat_qc'].script)
include {r as integration_seurat_sexfilt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_seurat_sexfilt'], script: modules['integration_seurat_sexfilt'].script)


include {r as integration_STACAS} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_STACAS'], script: modules['integration_STACAS'].script)
include {r as integration_STACAS_qc} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_STACAS_qc'], script: modules['integration_STACAS_qc'].script)
include {r as integration_STACAS_sexfilt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_STACAS_sexfilt'], script: modules['integration_STACAS_sexfilt'].script)


// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    metadata(params.input)

    // Run pipeline on seurat integration
    integration_seurat( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_seurat_qc( integration_seurat.out )
    integration_seurat_sexfilt( integration_seurat_qc.out )

    // Run pipeline on STACAS integration
    integration_STACAS( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_STACAS_qc( integration_STACAS.out )
    integration_STACAS_sexfilt( integration_STACAS_qc.out )
}