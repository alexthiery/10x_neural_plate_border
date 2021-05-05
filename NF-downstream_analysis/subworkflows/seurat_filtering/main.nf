/*
 * Merge loom and run scVelo
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

if(params.integration && params.integration != 'Seurat' && params.integration != 'STACAS'){
    log.error "'${params.integration}' is not a valid option for integration. Please use 'Seurat' or 'STACAS'."
    System.exit(1)
}

def analysis_scripts = [:]
analysis_scripts.integration        = params.integration == 'STACAS' ? file("$baseDir/bin/seurat/1_integration_STACAS.R", checkIfExists: true) : file("$baseDir/bin/seurat/1_integration.R", checkIfExists: true)
analysis_scripts.integration_qc     = file("$baseDir/bin/seurat/2_integration_qc.R", checkIfExists: true)
analysis_scripts.poor_cluster_filt  = file("$baseDir/bin/seurat/3_poor_cluster_filt.R", checkIfExists: true)
analysis_scripts.sex_filt           = file("$baseDir/bin/seurat/4_sex_filt.R", checkIfExists: true)
analysis_scripts.cell_cycle         = file("$baseDir/bin/seurat/5_cell_cycle.R", checkIfExists: true)
analysis_scripts.contamination_filt = file("$baseDir/bin/seurat/6_contamination_filt.R", checkIfExists: true)
// analysis_scripts.gene_modules       = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.seurat_h5ad     = file("$baseDir/bin/seurat/seurat_h5ad.R", checkIfExists: true)



params.integration_options          = [:]
params.integration_qc_options       = [:]
params.poor_cluster_filt_options    = [:]
params.sex_filt_options             = [:]
params.cell_cycle_options           = [:]
params.contamination_filt_options   = [:]
// params.gene_modules                 = [:]
params.seurat_h5ad_options          = [:]

// Include Seurat R processes
include {r as INTEGRATION} from "$baseDir/../modules/tools/r/main.nf"               addParams(  options: params.integration_options,
                                                                                                script: analysis_scripts.integration )

include {r as INTEGRATION_QC} from "$baseDir/../modules/tools/r/main.nf"            addParams(  options: params.integration_qc_options,
                                                                                                script: analysis_scripts.integration_qc )

include {r as POOR_CLUSTER_FILT} from "$baseDir/../modules/tools/r/main.nf"         addParams(  options: params.poor_cluster_filt_options,
                                                                                                script: analysis_scripts.poor_cluster_filt )

include {r as SEX_FILT} from "$baseDir/../modules/tools/r/main.nf"                  addParams(  options: params.sex_filt_options,
                                                                                                script: analysis_scripts.sex_filt )

include {r as CELL_CYCLE} from "$baseDir/../modules/tools/r/main.nf"                addParams(  options: params.cell_cycle_options,
                                                                                                script: analysis_scripts.cell_cycle )

include {r as CONTAMINATION_FILT} from "$baseDir/../modules/tools/r/main.nf"        addParams(  options: params.contamination_filt_options,
                                                                                                script: analysis_scripts.contamination_filt )

include {r as CONTAMINATION_FILT_H5AD} from "$baseDir/../modules/tools/r/main.nf"   addParams(  options: params.seurat_h5ad_options,
                                                                                                script: analysis_scripts.seurat_h5ad )

// include {r as GENE_MODULES} from "$baseDir/../modules/tools/r/main.nf"          addParams(  options: params.gene_modules_options,
//                                                                                             script: analysis_scripts.gene_modules )



workflow SEURAT_FILTERING {
    take:
    metadata_out

    main:
    // Run Seurat pipeline
    INTEGRATION( metadata_out.filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' } )
    INTEGRATION_QC( INTEGRATION.out )
    POOR_CLUSTER_FILT( INTEGRATION_QC.out )
    SEX_FILT( POOR_CLUSTER_FILT.out )
    CELL_CYCLE( SEX_FILT.out )
    CONTAMINATION_FILT( CELL_CYCLE.out )

    emit:
    integration_out             = INTEGRATION.out
    integration_qc_out          = INTEGRATION_QC.out
    poor_cluster_filt_out       = POOR_CLUSTER_FILT.out
    sex_filt_out                = SEX_FILT.out
    cell_cycle_out              = CELL_CYCLE.out
    contamination_filt_out      = CONTAMINATION_FILT.out
    contamination_filt_h5ad_out = CONTAMINATION_FILT_H5AD.out
}