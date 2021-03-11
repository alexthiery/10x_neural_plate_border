#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"
include {r_analysis as seurat_integrate} from "$baseDir/../modules/tools/r_analysis/main.nf" addParams(options: modules['seurat_integrate'], script: modules['seurat_integrate'].script)
include {r_analysis as seurat_2} from "$baseDir/../modules/tools/r_analysis/main.nf" addParams(options: modules['seurat_2'], script: modules['seurat_2'].script)


// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {

    metadata(params.input)

    //  Run seurat_integrate
    seurat_integrate( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )

    //  Run seurat_2
    seurat_2( seurat_integrate.out )
}