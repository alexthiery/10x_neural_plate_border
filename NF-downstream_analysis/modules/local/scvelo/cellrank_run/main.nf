process CELLRANK_RUN {
    tag "$meta.sample_id"
    label 'process_high'

    container "quay.io/biocontainers/cellrank:1.5.1--pyhdfd78af_0"

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