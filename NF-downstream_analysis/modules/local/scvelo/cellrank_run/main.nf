process CELLRANK_RUN {
    tag "$meta.sample_id"
    label 'process_high'

    container "alexthiery/10x-npb-scvelo:base-2.0.1"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_cellrank.h5ad"), emit: h5ad
    tuple val(meta), path("*_metadata.csv"), emit: csv
    path "figures", emit: plots

    script:
    def args = task.ext.args  ?: ''
    def prefix   = task.ext.prefix ?: "${meta.sample_id}"

    """
    export HDF5_USE_FILE_LOCKING=FALSE
    $moduleDir/bin/cellrank_run.py --input ${h5ad} --output ${prefix} --ncores ${task.cpus} ${args} 
    """
}