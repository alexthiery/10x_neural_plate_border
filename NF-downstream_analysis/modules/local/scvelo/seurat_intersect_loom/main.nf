process SEURAT_INTERSECT_LOOM {
    tag "$meta.sample_id"
    label 'process_medium'
    // maxForks 1

    container "alexthiery/10x-npb-scvelo:base-2.0.0"

    input:
    tuple val(meta), path(seurat), path(loom), path(annotations)

    output:
    tuple val(meta), path("*.loom"), emit: loom

    script:
    def prefix   = task.ext.prefix ?: "${meta.sample_id}"
    
    """
    export HDF5_USE_FILE_LOCKING=FALSE
    $moduleDir/bin/seurat_intersect_loom.py --loomInput ${loom} --seuratInput ${seurat} --annotations ${annotations} --output ${prefix}_seurat_intersect.loom
    """
}