#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

if(params.integration && params.integration != 'Seurat' && params.integration != 'STACAS'){
    log.error "'${params.integration}' is not a valid option for integration. Please use 'Seurat' or 'STACAS'."
    System.exit(1)
}

def analysis_scripts = [:]
analysis_scripts.integration = params.integration == 'STACAS' ? file("$baseDir/bin/seurat/1_integration_STACAS.R", checkIfExists: true) : file("$baseDir/bin/seurat/1_integration.R", checkIfExists: true)
analysis_scripts.integration_qc = file("$baseDir/bin/seurat/2_integration_qc.R", checkIfExists: true)
analysis_scripts.poor_cluster_filt = file("$baseDir/bin/seurat/3_poor_cluster_filt.R", checkIfExists: true)
analysis_scripts.sex_filt = file("$baseDir/bin/seurat/4_sex_filt.R", checkIfExists: true)
analysis_scripts.cell_cycle = file("$baseDir/bin/seurat/5_cell_cycle.R", checkIfExists: true)
analysis_scripts.contamination_filt = file("$baseDir/bin/seurat/6_contamination_filt.R", checkIfExists: true)
analysis_scripts.seurat_to_h5ad = file("$baseDir/bin/seurat/seurat_to_h5ad.R", checkIfExists: true)
analysis_scripts.gene_modules = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {metadata} from "$baseDir/../modules/tools/metadata/main.nf"

// Include Seurat R processes
include {r as integration} from "$baseDir/../modules/tools/r/main.nf"           addParams(  options: modules['integration'],
                                                                                            script: analysis_scripts.integration )

include {r as integration_qc} from "$baseDir/../modules/tools/r/main.nf"        addParams(  options: modules['integration_qc'],
                                                                                            script: analysis_scripts.integration_qc )

include {r as poor_cluster_filt} from "$baseDir/../modules/tools/r/main.nf"     addParams(  options: modules['poor_cluster_filt'],
                                                                                            script: analysis_scripts.poor_cluster_filt )

include {r as sex_filt} from "$baseDir/../modules/tools/r/main.nf"              addParams(  options: modules['sex_filt'],
                                                                                            script: analysis_scripts.sex_filt )

include {r as cell_cycle} from "$baseDir/../modules/tools/r/main.nf"            addParams(  options: modules['cell_cycle'],
                                                                                            script: analysis_scripts.cell_cycle )

include {r as contamination_filt} from "$baseDir/../modules/tools/r/main.nf"    addParams(  options: modules['contamination_filt'],
                                                                                            script: analysis_scripts.contamination_filt )

// Include other downstream processes
include {r as gene_modules} from "$baseDir/../modules/tools/r/main.nf"          addParams(  options: modules['gene_modules'],
                                                                                            script: analysis_scripts.gene_modules )

include {r as contamination_filt_h5ad} from "$baseDir/../modules/tools/r/main.nf" addParams(  options: modules['contamination_filt_h5ad'],
                                                                                                script: analysis_scripts.seurat_to_h5ad )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {
    log.info Headers.build_debug_param_summary(params, params.monochrome_logs)
    log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)
}


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
    contamination_filt_h5ad( contamination_filt.out )

    // Run downstream analyses
    gene_modules( contamination_filt.out )
}