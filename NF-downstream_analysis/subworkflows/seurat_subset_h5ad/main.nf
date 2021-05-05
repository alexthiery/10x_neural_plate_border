/*
 * Convert seurat to h5ad format for python
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/
params.contamination_filt_h5ad_options   = [:]

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/
include {SEURAT_H5AD as CONTAMINATION_FILT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"   addParams(  options: params.contamination_filt_h5ad_options)

/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_SUBSET_H5AD {
    take:
    contamination_filt_out

    main:
    // Run Seurat pipeline
    CONTAMINATION_FILT_H5AD( contamination_filt_out )

    emit:
    contamination_filt_h5ad_out = CONTAMINATION_FILT_H5AD.out
}