#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"

// Include Seurat R processes
include {r as integration} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration'], script: modules['integration'].script)
include {r as integration_qc} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['integration_qc'], script: modules['integration_qc'].script)
include {r as sex_filt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['sex_filt'], script: modules['sex_filt'].script)
include {r as cell_cycle} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['cell_cycle'], script: modules['cell_cycle'].script)
include {r as contamination_filt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['contamination_filt'], script: modules['contamination_filt'].script)
include {r as poor_cluster_filt} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['poor_cluster_filt'], script: modules['poor_cluster_filt'].script)

// Include other downstream processes
include {r as gene_modules} from "$baseDir/../modules/tools/r/main.nf" addParams(options: modules['gene_modules'], script: modules['gene_modules'].script)


// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    metadata(params.input)

    // Run Seurat pipeline
    integration( metadata.out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    integration_qc( integration.out )
    poor_cluster_filt( integration_qc.out )
    sex_filt( poor_cluster_filt.out )
    cell_cycle( sex_filt.out )
    contamination_filt( cell_cycle.out )

    // Run downstream analyses
    gene_modules( contamination_filt.out )
}