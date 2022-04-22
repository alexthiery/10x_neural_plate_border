process SEURAT_H5AD {
    tag "$meta.sample_id"
    label 'process_medium'

    container "alexthiery/10x-npb-schelper:latest"

    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*h5ad')

    script:
    def args = task.ext.args  ?: ''

    """
    Rscript $moduleDir/bin/seurat_h5ad.R ${args}
    """
}
