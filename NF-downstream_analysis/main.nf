#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def skip_seurat_filtering = params.skip_seurat_filtering ? true : false
// def skip_scvelo = params.skip_scvelo ? true : false

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

def analysis_scripts                                = [:]
analysis_scripts.transfer_labels                    = file("$baseDir/bin/seurat/transfer_labels.R", checkIfExists: true)
analysis_scripts.gene_modules_latent_time           = file("$baseDir/bin/other/gene_modules_latent_time.R", checkIfExists: true)
analysis_scripts.refined_gene_modules_latent_time   = file("$baseDir/bin/other/refined_gene_modules_latent_time.R", checkIfExists: true)


/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/

include {METADATA} from "$baseDir/subworkflows/metadata/main"

include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"                                                   addParams(  preprocessing_options:                  modules['preprocessing'],
                                                                                                                                            integration_options:                    modules['integration'],
                                                                                                                                            integration_qc_options:                 modules['integration_qc'],
                                                                                                                                            sex_filt_options:                       modules['sex_filt'],
                                                                                                                                            cell_cycle_options:                     modules['cell_cycle'],
                                                                                                                                            contamination_filt_options:             modules['contamination_filt'])

// Modules and subworkflows for running scVelo/cellrank                                             
include {MERGE_LOOM} from "$baseDir/modules/local/merge_loom/main"                                                              addParams(  options:                                modules['merge_loom'])

include {SEURAT_FILTERED_PROCESS} from "$baseDir/subworkflows/seurat_filtered_process/main"                                     addParams(  scatterplot3d_options:                  modules['scatterplot3d'],
                                                                                                                                            gene_module_options:                    modules['gene_modules'],
                                                                                                                                            state_classification_options:           modules['state_classification'],
                                                                                                                                            phate_options:                          modules['phate'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['scvelo_run'],
                                                                                                                                            cellrank_run_options:                   modules['cellrank_run'])

// Subworkflows for subset stage, run and clusters from seurat object and run downstream analysis pipelines
include {SEURAT_SPLIT_PROCESS as SEURAT_STAGE_PROCESS} from "$baseDir/subworkflows/seurat_split_process/main"                   addParams(  split_options:                          modules['stage_split'],
                                                                                                                                            cluster_options:                        modules['stage_cluster'],
                                                                                                                                            gene_modules_options:                   modules['stage_gene_modules'],
                                                                                                                                            state_classification_options:           modules['stage_state_classification'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['stage_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['stage_scvelo_run'])

include {SEURAT_SPLIT_PROCESS as SEURAT_RUN_PROCESS} from "$baseDir/subworkflows/seurat_split_process/main"                     addParams(  split_options:                          modules['run_split'],
                                                                                                                                            cluster_options:                        modules['run_cluster'],
                                                                                                                                            gene_modules_options:                   modules['run_gene_modules'],
                                                                                                                                            state_classification_options:           modules['run_state_classification'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['run_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['run_scvelo_run'])

include {SEURAT_SUBSET_PROCESS as SEURAT_NPB_PROCESS} from "$baseDir/subworkflows/seurat_subset_process/main"                   addParams(  subset_options:                         modules['npb_subset'],
                                                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                                                            phate_options:                          modules['clusters_phate'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['clusters_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['clusters_scvelo_run'])

include {SEURAT_SUBSET_PROCESS as SEURAT_PLACODAL1_PROCESS} from "$baseDir/subworkflows/seurat_subset_process/main"             addParams(  subset_options:                         modules['placodal1_subset'],
                                                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                                                            phate_options:                          modules['clusters_phate'],
                                                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['clusters_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['clusters_scvelo_run'])

include {SEURAT_SUBSET_PROCESS as SEURAT_PLACODAL2_PROCESS} from "$baseDir/subworkflows/seurat_subset_process/main"             addParams(  subset_options:                         modules['placodal2_subset'],
                                                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                                                            phate_options:                          modules['clusters_phate'],
                                                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['clusters_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['clusters_scvelo_run'])

include {SEURAT_SUBSET_CELLRANK_PROCESS as SEURAT_FILTER_CONTAM_PROCESS} from "$baseDir/subworkflows/seurat_subset_cellrank_process/main"         addParams(  subset_options:                         modules['filter_contam_subset'],
                                                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                                                            phate_options:                          modules['clusters_phate'],
                                                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['clusters_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['filter_contam_scvelo_run'],
                                                                                                                                            cellrank_run_options:                   modules['clusters_cellrank_run'])


include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                                                             addParams(  options:                                modules['transfer_labels'],
                                                                                                                                            script:                                 analysis_scripts.transfer_labels )

include {SEURAT_TRANSFER_PROCESS as SEURAT_TRANSFER_NPB_PROCESS} from "$baseDir/subworkflows/seurat_transfer_process/main"      addParams(  subset_options:                         modules['transfer_subset_npb'],
                                                                                                                                            cluster_options:                        modules['transfer_subset_cluster'],
                                                                                                                                            gene_modules_options:                   modules['transfer_subset_gene_modules'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['transfer_subset_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['transfer_subset_scvelo_run'])

include {SEURAT_TRANSFER_PROCESS as SEURAT_TRANSFER_FILTER_PROCESS} from "$baseDir/subworkflows/seurat_transfer_process/main"   addParams(  subset_options:                         modules['transfer_subset_filter'],
                                                                                                                                            cluster_options:                        modules['transfer_subset_cluster'],
                                                                                                                                            gene_modules_options:                   modules['transfer_subset_gene_modules'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['transfer_subset_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['transfer_subset_scvelo_run'])

include {SEURAT_TRANSFER_FULL_PROCESS} from "$baseDir/subworkflows/seurat_transfer_full_process/main"                           addParams(  gene_modules_options:                   modules['transfer_labels_gene_modules'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['transfer_labels_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['transfer_labels_scvelo_run'],
                                                                                                                                            cellrank_run_options:                   modules['transfer_labels_cellrank_run'])

include {SEURAT_TRANSFER_FULL_PROCESS as REFINED_FULL_PROCESS} from "$baseDir/subworkflows/seurat_transfer_full_process/main"   addParams(  gene_modules_options:                   modules['transfer_labels_gene_modules'],
                                                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                                                            seurat_intersect_loom_options:          modules['transfer_labels_seurat_intersect_loom'],
                                                                                                                                            scvelo_run_options:                     modules['transfer_labels_scvelo_run'],
                                                                                                                                            cellrank_run_options:                   modules['refined_cellrank_run'])

include {R as GENE_MODULES_LATENT_TIME} from "$baseDir/modules/local/r/main"                                                    addParams( options:                                modules['gene_modules_latent_time'],
                                                                                                                                            script:                                 analysis_scripts.gene_modules_latent_time )

include {R as REFINED_GENE_MODULES_LATENT_TIME} from "$baseDir/modules/local/r/main"                                            addParams(  options:                                modules['refined_gene_modules_latent_time'],
                                                                                                                                            script:                                 analysis_scripts.refined_gene_modules_latent_time )


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
    SEURAT_FILTERED_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on stage and run split
    --------------------------------------------------------------------------------------*/ 
    SEURAT_STAGE_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )
    SEURAT_RUN_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on cluster subsets
    --------------------------------------------------------------------------------------*/     
    SEURAT_NPB_PROCESS( SEURAT_FILTERED_PROCESS.out.state_classification_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )
    SEURAT_PLACODAL1_PROCESS( SEURAT_FILTERED_PROCESS.out.state_classification_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )
    SEURAT_PLACODAL2_PROCESS( SEURAT_FILTERED_PROCESS.out.state_classification_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )
    // SEURAT_FILTER_CONTAM_PROCESS( SEURAT_FILTERED_PROCESS.out.state_classification_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    // Collect rds files from all stages
    ch_combined = SEURAT_STAGE_PROCESS.out.state_classification_out
        .concat(SEURAT_FILTERING.out.contamination_filt_out)
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'all_stages_filtered'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    // Transfer labels from stage subsets to full data
    TRANSFER_LABELS( ch_combined )

    SEURAT_TRANSFER_NPB_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    SEURAT_TRANSFER_FILTER_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    SEURAT_TRANSFER_FULL_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )
    

    ch_full_state_classification    = SEURAT_FILTERED_PROCESS.out.state_classification_out.map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
    ch_full_cellrank                = SEURAT_FILTERED_PROCESS.out.cellrank_run_out_metadata.map{it[1]}

    ch_full_latent_time             = SEURAT_FILTERED_PROCESS.out.gene_modules_out
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
                                        .combine(ch_full_state_classification)
                                        .combine(ch_full_cellrank)
                                        .map{[[sample_id:'full_gm_latent_time'], it]}

    ch_stage_latent_time            = SEURAT_STAGE_PROCESS.out.gene_modules_out
                                        .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]}
                                        .combine(ch_full_state_classification)
                                        .combine(ch_full_cellrank)
                                        .map{[[sample_id:it[0].sample_id.split("_")[0]+'_gm_latent_time'], [it[1], it[2], it[3]]]}
    
    GENE_MODULES_LATENT_TIME( ch_full_latent_time.concat(ch_stage_latent_time) )




    // Run cellrank and gm dynamics with refined terminal states
    REFINED_FULL_PROCESS( TRANSFER_LABELS.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} )

    ch_transfer_state_classification    = TRANSFER_LABELS.out.map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
    ch_transfer_cellrank                = REFINED_FULL_PROCESS.out.cellrank_run_out_metadata.map{it[1]}

    ch_transfer_latent_time             = REFINED_FULL_PROCESS.out.gene_modules_out
                                        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
                                        .combine(ch_transfer_state_classification)
                                        .combine(ch_transfer_cellrank)
                                        .map{[[sample_id:'refined_gm_latent_time'], it]}

    REFINED_GENE_MODULES_LATENT_TIME( ch_transfer_latent_time )
}