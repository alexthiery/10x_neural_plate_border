#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"
include {r as seurat_integrate} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['seurat_integrate'], script: modules['seurat_integrate'].script)
include {r as seurat_sexfilt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['seurat_sexfilt'], script: modules['seurat_sexfilt'].script)
include {r as seurat_sexfilt_log} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['seurat_sexfilt_log'], script: modules['seurat_sexfilt_log'].script)


// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {

    metadata(params.input)

    //  Run seurat_integrate
    seurat_integrate( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )


    //  Run seurat_integrate
    seurat_integrate_2( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )

    // //  Run seurat_2
    // seurat_sexfilt( seurat_integrate.out )

    // //  Run seurat_2
    // seurat_sexfilt_log( seurat_integrate.out )
}