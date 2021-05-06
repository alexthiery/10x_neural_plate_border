#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def merge_loom_options = modules['merge_loom']
merge_loom_options.skip_process = file(params.loomInput).isFile()

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

include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"       addParams(  integration_options:                modules['integration'],
                                                                                                integration_qc_options:             modules['integration_qc'],
                                                                                                poor_cluster_filt_options:          modules['poor_cluster_filt'],
                                                                                                sex_filt_options:                   modules['sex_filt'],
                                                                                                cell_cycle_options:                 modules['cell_cycle'],
                                                                                                contamination_filt_options:         modules['contamination_filt'] )


include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"             addParams(  merge_loom_options:                 merge_loom_options,
                                                                                                seurat_intersect_loom_options:      modules['seurat_intersect_loom'],
                                                                                                scvelo_options:                     modules['scvelo'] )

include {SEURAT_SUBSET_H5AD} from "$baseDir/subworkflows/seurat_subset_h5ad/main"   addParams(  contamination_filt_h5ad_options:    modules['contamination_filt_h5ad'])


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

        seurat_out = SEURAT_FILTERING.out.contamination_filt_out
        seurat_annotations = SEURAT_FILTERING.out.annotations
   } else {
       seurat_out = [[sample_id:'temp'], [params.seurat_out]]
       seurat_annotations = params.seurat_annotations
   }

    /*------------------------------------------------------------------------------------*/
    /* Prepare inputs for scVelo
    --------------------------------------------------------------------------------------*/
    // Convert seurat to h5ad format
    SEURAT_SUBSET_H5AD( seurat_out )

    if(!skip_scvelo){
        // Set channel for input looms
        METADATA.out
            .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
            .map {it[1].collect{ file(it+"/velocyto", checkIfExists: true) }}
            .set {ch_loomInput}

        SEURAT_SCVELO( ch_loomInput, SEURAT_SUBSET_H5AD.out.contamination_filt_h5ad_out.map{it[1]}, seurat_annotations )
    }

}