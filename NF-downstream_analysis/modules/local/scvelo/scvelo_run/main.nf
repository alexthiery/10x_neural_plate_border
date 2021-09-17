#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCVELO_RUN {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    container "alexthiery/10x-npb-scvelo:base-1.8"

    input:
        tuple val(meta), path(loom)

    output:
        tuple val(meta), path("*_scvelo.h5ad"), emit: h5ad
        tuple val(meta), path("*_metadata.csv"), emit: csv
        path "figures", emit: plots

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "${meta.sample_id}"

        """
        export HDF5_USE_FILE_LOCKING=FALSE
        $moduleDir/bin/scvelo_run.py --input ${loom} --output ${prefix} --ncores ${task.cpus} ${options.args} 
        """
}