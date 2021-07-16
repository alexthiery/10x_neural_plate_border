#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def skip_seurat_filtering = params.skip_seurat_filtering ? true : false
def skip_scvelo = params.skip_scvelo ? true : false

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/

include {METADATA} from "$baseDir/subworkflows/metadata/main"

include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"           addParams(  preprocessing_options:              modules['preprocessing'],
                                                                                                    integration_options:                modules['integration'],
                                                                                                    integration_qc_options:             modules['integration_qc'],
                                                                                                    sex_filt_options:                   modules['sex_filt'],
                                                                                                    cell_cycle_options:                 modules['cell_cycle'],
                                                                                                    contamination_filt_options:         modules['contamination_filt'] )

include {EXPLORATORY_ANALYSIS} from "$baseDir/subworkflows/exploratory_analysis/main"   addParams(  gene_module_options:                modules['gene_modules'],
                                                                                                    scatterplot3d_options: modules['scatterplot3'])

include {SEURAT_STAGE_PROCESS} from "$baseDir/subworkflows/seurat_stage_process/main"   addParams(  stage_split_options:                modules['stage_split'],
                                                                                                    stage_cluster_options:              modules['stage_cluster'],
                                                                                                    stage_gene_modules_options:         modules['stage_gene_modules'])

include {SEURAT_RUN_PROCESS} from "$baseDir/subworkflows/seurat_run_process/main"   addParams(  run_split_options:                modules['run_split'])
                                                                                                    // run_cluster_options:              modules['run_cluster'],
                                                                                                    // run_gene_modules_options:         modules['run_gene_modules'])
                                                                                                    
include {MERGE_LOOM} from "$baseDir/modules/local/merge_loom/main"                      addParams(  options:                            modules['merge_loom'] )

include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"                 addParams(  seurat_intersect_loom_options:      modules['seurat_intersect_loom'],
                                                                                                    scvelo_run_options:                 modules['scvelo_run'] )
                                                                                                    
include {SEURAT_SUBSET_H5AD} from "$baseDir/subworkflows/seurat_subset_h5ad/main"       addParams(  contamination_filt_h5ad_options:    modules['contamination_filt_h5ad'] )

include {EXPLORATORY_LATENT_TIME} from "$baseDir/subworkflows/exploratory_latent_time/main" addParams( gene_modules_latent_time_options: modules['gene_modules_latent_time'],
                                                                                                       cell_state_classification_options: modules['cell_state_classification'])

include {EXPLORATORY_LATENT_TIME as STAGE_EXPLORATORY_LATENT_TIME} from "$baseDir/subworkflows/exploratory_latent_time/main" addParams(  gene_modules_latent_time_options: modules['stage_gene_modules_latent_time'])


workflow {
    METADATA( params.input )

    /*------------------------------------------------------------------------------------*/
    /* Run inital seurat pipeline
    --------------------------------------------------------------------------------------*/
   if(!skip_seurat_filtering){
        // Set channel for cellranger counts
        METADATA.out
            .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
            .map {[it[0], it[1].collect{ file(it+"/cellranger/count/filtered_feature_bc_matrix", checkIfExists: true) }]}
            .set {ch_scRNAseq_counts}

        SEURAT_FILTERING( ch_scRNAseq_counts )
        
        SEURAT_STAGE_PROCESS( SEURAT_FILTERING.out.contamination_filt_out )

        SEURAT_RUN_PROCESS( SEURAT_FILTERING.out.contamination_filt_out )

        EXPLORATORY_ANALYSIS( SEURAT_FILTERING.out.contamination_filt_out )
        
        // Convert seurat to h5ad format
        SEURAT_SUBSET_H5AD( SEURAT_FILTERING.out.contamination_filt_out.concat(SEURAT_STAGE_PROCESS.out.stage_cluster_out) )
        
        ch_seurat_h5ad = SEURAT_SUBSET_H5AD.out.contamination_filt_h5ad_out
        ch_seurat_annotations = SEURAT_FILTERING.out.annotations


   } else {
       seurat_h5ad = [[[sample_id:'NF-scRNAseq_alignment_out'], file(params.seurat_h5ad, checkIfExists: true)]]
       ch_seurat_h5ad = Channel.from(seurat_h5ad)

       seurat_annotations = [[[sample_id:'NF-scRNAseq_alignment_out'], file(params.seurat_annotations, checkIfExists: true)]]
       ch_seurat_annotations = Channel.from(seurat_annotations)
   }

    /*------------------------------------------------------------------------------------*/
    /* Prepare inputs for scVelo
    --------------------------------------------------------------------------------------*/


    if(!skip_scvelo){
        // ch_seurat_h5ad = ch_seurat_h5ad.filter{it[0].sample_id == 'NF-scRNAseq_alignment_out'}
        // Set channel for input looms
        METADATA.out
            .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
            .map {[it[0], it[1].collect{ file(it+"/velocyto", checkIfExists: true) }]}
            .set {ch_loomInput}

        MERGE_LOOM( ch_loomInput )
        
        SEURAT_SCVELO( ch_seurat_h5ad, MERGE_LOOM.out.loom.map{it[1]}, ch_seurat_annotations.map{it[1]} ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv
    }

    ch_seurat_data = SEURAT_FILTERING.out.contamination_filt_out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    ch_antler_data = EXPLORATORY_ANALYSIS.out.gene_modules_out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    EXPLORATORY_LATENT_TIME(ch_seurat_data, ch_antler_data, SEURAT_SCVELO.out.scvelo_run_out_metadata)

    ch_seurat_stage_data = SEURAT_STAGE_PROCESS.out.stage_cluster_out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    ch_antler_stage_data = SEURAT_STAGE_PROCESS.out.stage_gene_modules_out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    STAGE_EXPLORATORY_LATENT_TIME(ch_seurat_stage_data, ch_antler_stage_data, SEURAT_SCVELO.out.scvelo_run_out_metadata)
}