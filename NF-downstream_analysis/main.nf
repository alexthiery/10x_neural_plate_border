#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

def analysis_scripts                                = [:]
analysis_scripts.transfer_labels                    = file("$baseDir/bin/seurat/transfer_labels.R", checkIfExists: true)
analysis_scripts.gene_modules_subset_latent_time    = file("$baseDir/bin/other/gene_modules_subset_latent_time.R", checkIfExists: true)
analysis_scripts.gene_modules_npb_latent_time       = file("$baseDir/bin/other/gene_modules_npb_latent_time.R", checkIfExists: true)
analysis_scripts.coexpression_analysis_npb          = file("$baseDir/bin/other/coexpression_analysis_npb.R", checkIfExists: true)
analysis_scripts.coexpression_nc_ppr_modules_npb    = file("$baseDir/bin/other/coexpression_nc_ppr_modules_npb.R", checkIfExists: true)
analysis_scripts.identify_lineage_drivers           = file("$baseDir/bin/other/identify_lineage_drivers.R", checkIfExists: true)

/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/

include {METADATA} from "$baseDir/subworkflows/metadata/main"
include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"

// Modules and subworkflows for running scVelo/cellrank                                             
include {MERGE_LOOM} from "$baseDir/modules/local/merge_loom/main"
include {SEURAT_FILTERED_PROCESS} from "$baseDir/subworkflows/seurat_filtered_process/main"

// Subworkflows for split stage and run and run downstream analysis
include {SEURAT_SPLIT_PROCESS as SEURAT_STAGE_PROCESS} from "$baseDir/subworkflows/seurat_split_process/main"

// Subworkflows and modules for label transfer and subsequent cluster subsets and run downstream analysis
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                                                            addParams( script: analysis_scripts.transfer_labels )
include {SEURAT_TRANSFER_FULL_PROCESS} from "$baseDir/subworkflows/seurat_transfer_full_process/main"
include {SEURAT_TRANSFER_PROCESS as SEURAT_TRANSFER_PPR_NC_PROCESS} from "$baseDir/subworkflows/seurat_transfer_process/main"
include {R as GENE_MODULES_SUBSET_LATENT_TIME} from "$baseDir/modules/local/r/main"                                            addParams( script: analysis_scripts.gene_modules_subset_latent_time )
include {R as GENE_MODULES_NPB_LATENT_TIME} from "$baseDir/modules/local/r/main"                                               addParams( script: analysis_scripts.gene_modules_npb_latent_time )
include {R as COEXPRESSION_ANALYSIS_NPB} from "$baseDir/modules/local/r/main"                                                  addParams( script: analysis_scripts.coexpression_analysis_npb )
include {R as COEXPRESSION_NC_PPR_MODULES_NPB} from "$baseDir/modules/local/r/main"                                            addParams( script: analysis_scripts.coexpression_nc_ppr_modules_npb )
include {R as IDENTIFY_LINEAGE_DRIVERS_FULL} from  "$baseDir/modules/local/r/main"                                             addParams( script: analysis_scripts.identify_lineage_drivers )
include {R as IDENTIFY_LINEAGE_DRIVERS_NPB} from  "$baseDir/modules/local/r/main"                                              addParams( script: analysis_scripts.identify_lineage_drivers )

// HCR intensity subworkflow
include {HCR} from "$baseDir/subworkflows/hcr/main"

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}

/*------------------------------------------------------------------------------------
Set channels
--------------------------------------------------------------------------------------*/

// Set channel for binary knowledge matrix for cell state classification
Channel
    .value("$baseDir/binary_knowledge_matrix.csv")
    .set{ch_binary_knowledge_matrix}

hcr_intensity_samplesheet = "$baseDir/hcr_data/samplesheet.csv"

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow {
    METADATA( params.input )

    /*------------------------------------------------------------------------------------*/
    /* Run inital seurat pipeline
    --------------------------------------------------------------------------------------*/
    // Set channel for cellranger counts
    METADATA.out
        .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
        .map {[it[0], it[1].collect{ file(it+"/cellranger/count/filtered_feature_bc_matrix", checkIfExists: true) }]}
        .set {ch_scRNAseq_counts}

    SEURAT_FILTERING( ch_scRNAseq_counts )
        
    /*------------------------------------------------------------------------------------*/
    /* Prepare inputs for scVelo
    --------------------------------------------------------------------------------------*/
   
    // Set channel for input looms
    METADATA.out
        .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
        .map {[it[0], it[1].collect{ file(it+"/velocyto", checkIfExists: true) }]}
        .set {ch_loomInput}

    MERGE_LOOM( ch_loomInput )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on full filtered seurat object
    --------------------------------------------------------------------------------------*/

    SEURAT_FILTERED_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]}, ch_binary_knowledge_matrix )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on stage split
    --------------------------------------------------------------------------------------*/ 
    SEURAT_STAGE_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]}, ch_binary_knowledge_matrix )


    /*------------------------------------------------------------------------------------*/
    /* Transfer cell type labels from stage to full dataset
    --------------------------------------------------------------------------------------*/     

    // Collect rds files from all stages
    ch_combined = SEURAT_STAGE_PROCESS.out.state_classification_out
        .concat(SEURAT_FILTERING.out.contamination_filt_out)
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'all_stages_filtered'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    // Transfer labels from stage subsets to full data
    TRANSFER_LABELS( ch_combined )

    SEURAT_TRANSFER_FULL_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on cluster subsets after transfer labels
    --------------------------------------------------------------------------------------*/
    SEURAT_TRANSFER_PPR_NC_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    ch_full_state_classification    = TRANSFER_LABELS.out.map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
    ch_full_cellrank                = SEURAT_TRANSFER_FULL_PROCESS.out.cellrank_run_out_metadata.map{it[1]}

    ch_full_latent_time             = SEURAT_TRANSFER_FULL_PROCESS.out.gene_modules_out
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
                                        .combine(ch_full_state_classification)
                                        .combine(ch_full_cellrank)
                                        .map{[[sample_id:'full_gm_latent_time'], it]}
    
    GENE_MODULES_SUBSET_LATENT_TIME( ch_full_latent_time )

    ch_seurat_npb_subset            = SEURAT_TRANSFER_PPR_NC_PROCESS.out.cluster_out.map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}

    ch_npb_cellrank                = SEURAT_TRANSFER_PPR_NC_PROCESS.out.cellrank_run_out_metadata.map{it[1]}

    ch_npb_latent_time              = SEURAT_TRANSFER_PPR_NC_PROCESS.out.gene_modules_out
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
                                        .combine(ch_seurat_npb_subset)
                                        .combine(ch_npb_cellrank)
                                        .combine(ch_binary_knowledge_matrix)
                                        .map{[[sample_id:'npb_gm_latent_time'], it]}

    GENE_MODULES_NPB_LATENT_TIME( ch_npb_latent_time )


    // Run coexpression analysis on modules calculated from NPB subset

    COEXPRESSION_NC_PPR_MODULES_NPB( ch_npb_latent_time )


    // Set up channels for running co-expression analysis on npb subset with gms from ss8 //[[meta], ['HH5.rds', 'HH6.rds' … 'ss8.rds', 'npb_subset.rds', 'ss8_antler.rds']]
    ch_stage_data                   = SEURAT_STAGE_PROCESS.out
                                        .state_classification_out
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]} // it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]
                                        .collect()
    
    ch_ss8_gms                      = SEURAT_STAGE_PROCESS.out
                                        .gene_modules_out
                                        .filter{ it[0].sample_id == 'ss8_splitstage_data' }
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}

    ch_coexpression_analysis_npb    = ch_ss8_gms
                                        .combine(ch_stage_data)
                                        .combine(ch_seurat_npb_subset)
                                        .map{[[sample_id:'coexpression_analysis_npb'], it]}


    COEXPRESSION_ANALYSIS_NPB(ch_coexpression_analysis_npb)

    /*------------------------------------------------------------------------------------*/
    /* Run HCR intensity subworkflow
    --------------------------------------------------------------------------------------*/

    hh7_seurat                      = SEURAT_STAGE_PROCESS.out
                                        .state_classification_out
                                        .filter{ it[0].sample_id == 'HH7_splitstage_data' }
                                        .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]}

    HCR(hcr_intensity_samplesheet, hh7_seurat)

    /*------------------------------------------------------------------------------------*/
    /* Identify lineage drivers
    --------------------------------------------------------------------------------------*/
    
    GENE_MODULES_SUBSET_LATENT_TIME.out
        .combine(SEURAT_FILTERING.out.annotations.map{it[1]})
        .map{meta, output, annotations_csv -> [meta, [output.findAll{it =~ /rds_files/}[0].listFiles()[0], annotations_csv]]}
        .set{ch_identify_lineage_drivers_full}
        
    IDENTIFY_LINEAGE_DRIVERS_FULL( ch_identify_lineage_drivers_full )

    // GENE_MODULES_NPB_LATENT_TIME.out
    //     .combine(SEURAT_FILTERING.out.annotations.map{it[1]})
    //     .map{meta, output, annotations_csv -> [meta, [output.findAll{it =~ /rds_files/}[0].listFiles()[0], annotations_csv]]}
    //     .set{ch_identify_lineage_drivers_npb}

    // IDENTIFY_LINEAGE_DRIVERS_NPB( ch_identify_lineage_drivers_npb )
}