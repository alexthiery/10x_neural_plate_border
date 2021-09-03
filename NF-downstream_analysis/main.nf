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
analysis_scripts.gene_modules_latent_time           = file("$baseDir/bin/other/gene_modules_latent_time.R", checkIfExists: true)


/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/

include {METADATA} from "$baseDir/subworkflows/metadata/main"

include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"                   addParams(  preprocessing_options:                  modules['preprocessing'],
                                                                                                            integration_options:                    modules['integration'],
                                                                                                            integration_qc_options:                 modules['integration_qc'],
                                                                                                            sex_filt_options:                       modules['sex_filt'],
                                                                                                            cell_cycle_options:                     modules['cell_cycle'],
                                                                                                            contamination_filt_options:             modules['contamination_filt'] )


// Modules and subworkflows for running scVelo/cellrank                                             
include {MERGE_LOOM} from "$baseDir/modules/local/merge_loom/main"                              addParams(  options:                                modules['merge_loom'] )

include {SEURAT_FILTERED_PROCESS} from "$baseDir/subworkflows/seurat_filtered_process/main"     addParams(  scatterplot3d_options:                  modules['scatterplot3d'],
                                                                                                            gene_module_options:                    modules['gene_modules'],
                                                                                                            state_classification_options:           modules['state_classification'],
                                                                                                            phate_options:                          modules['phate'])

// Subworkflows for subset stage, run and clusters from seurat object and run downstream analysis pipelines
include {SEURAT_STAGE_PROCESS} from "$baseDir/subworkflows/seurat_stage_process/main"           addParams(  split_options:                          modules['stage_split'],
                                                                                                            cluster_options:                        modules['stage_cluster'],
                                                                                                            gene_modules_options:                   modules['stage_gene_modules'],
                                                                                                            state_classification_options:           modules['stage_state_classification'])

include {SEURAT_RUN_PROCESS} from "$baseDir/subworkflows/seurat_run_process/main"               addParams(  split_options:                          modules['run_split'],
                                                                                                            cluster_options:                        modules['run_cluster'],
                                                                                                            gene_modules_options:                   modules['run_gene_modules'],
                                                                                                            state_classification_options:           modules['run_state_classification'])

include {SEURAT_CLUSTERS_PROCESS as SEURAT_NPB_PROCESS} from "$baseDir/subworkflows/seurat_clusters_process/main"     addParams(  subset_options:                         modules['npb_subset'],
                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                            state_classification_options:           modules['clusters_state_classification'])



include {SEURAT_CLUSTERS_PROCESS as SEURAT_HH4_PROCESS} from "$baseDir/subworkflows/seurat_clusters_process/main"     addParams(  subset_options:                         modules['hh4_subset'],
                                                                                                            cluster_options:                        modules['hh4_cluster'],
                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                            phate_options:                          modules['phate'])

include {SEURAT_CLUSTERS_PROCESS as SEURAT_NPB_HH4_PROCESS} from "$baseDir/subworkflows/seurat_clusters_process/main"     addParams(  subset_options:                         modules['npb_subset'],
                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                            phate_options:                          modules['phate'])

include {SEURAT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"                            addParams(  options:                                modules['seurat_h5ad'] )

include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"                         addParams(  seurat_intersect_loom_options:          modules['seurat_intersect_loom'],
                                                                                                            scvelo_run_options:                     modules['scvelo_run'])

include {R as GENE_MODULES_LATENT_TIME} from "$baseDir/modules/local/r/main"                    addParams(  options:                                modules['gene_modules_latent_time'],
                                                                                                            script:                                 analysis_scripts.gene_modules_latent_time )


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
    SEURAT_FILTERED_PROCESS( SEURAT_FILTERING.out.contamination_filt_out)

    /*------------------------------------------------------------------------------------*/
    /* Split data and cluster batches and stages separately (+GMs)
    --------------------------------------------------------------------------------------*/ 
    SEURAT_STAGE_PROCESS( SEURAT_FILTERING.out.contamination_filt_out)  
    SEURAT_RUN_PROCESS( SEURAT_FILTERING.out.contamination_filt_out)
    SEURAT_HH4_PROCESS( SEURAT_FILTERING.out.contamination_filt_out)
    SEURAT_NPB_PROCESS( SEURAT_FILTERED_PROCESS.out.state_classification_out)
    // SEURAT_NPB_HH4_PROCESS( SEURAT_NPB_PROCESS.out.state_classification_out )

    // Prepare outputs for scVelo
    ch_seurat_concat =          SEURAT_FILTERED_PROCESS.out.state_classification_out
                                    .concat(SEURAT_STAGE_PROCESS.out.cluster_out)
                                    .concat(SEURAT_RUN_PROCESS.out.cluster_out)
                                    .concat(SEURAT_HH4_PROCESS.out.cluster_out)
                                    .concat(SEURAT_NPB_PROCESS.out.cluster_out)
                                    // .concat(SEURAT_NPB_HH4_PROCESS.out.cluster_out)


    ch_gene_modules_concat =    SEURAT_FILTERED_PROCESS.out.gene_modules_out
                                    .concat(SEURAT_STAGE_PROCESS.out.gene_modules_out)
                                    .concat(SEURAT_RUN_PROCESS.out.gene_modules_out)
                                    .concat(SEURAT_HH4_PROCESS.out.gene_modules_out)
                                    .concat(SEURAT_NPB_PROCESS.out.gene_modules_out)
                                    // .concat(SEURAT_NPB_HH4_PROCESS.out.gene_modules_out)


    // Run scVelo
    SEURAT_H5AD( ch_seurat_concat )
    SEURAT_SCVELO( SEURAT_H5AD.out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]} ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv
    
    // Run gene module analysis across latent time
    ch_cluster_rds              = ch_seurat_concat.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    ch_gene_modules_rds         = ch_gene_modules_concat.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    ch_gene_module_latent_time  = ch_cluster_rds.combine(ch_gene_modules_rds, by: 0).combine(SEURAT_SCVELO.out.scvelo_run_out_metadata, by: 0)
    ch_gene_module_latent_time  = ch_gene_module_latent_time.map{[it[0], [it[1], it[2], it[3]]]}
    
    GENE_MODULES_LATENT_TIME( ch_gene_module_latent_time )
}