#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SCVELO_RUN {
    tag "$meta.sample_id"
    label 'process_high'

    container "quay.io/biocontainers/scvelo:0.2.4--pyhdfd78af_0"

    input:
    tuple val(meta), path(loom)

    output:
    tuple val(meta), path("*_scvelo.h5ad"), emit: h5ad
    path "figures", emit: plots

    script:
    def args = task.ext.args  ?: ''
    def prefix   = task.ext.prefix ?: "${meta.sample_id}"

    """
    export HDF5_USE_FILE_LOCKING=FALSE
    $moduleDir/bin/scvelo_run.py --input ${loom} --output ${prefix} --ncores ${task.cpus} ${args} 
    """
}