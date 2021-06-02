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
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    container "alexthiery/10x-npb-scvelo:base-1.5"

    input:
        tuple val(meta), path(loom)

    output:
        path "figures", emit: plots

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "seurat_merged"
        """
        $moduleDir/bin/scvelo_run.py --input ${loom} --ncores ${task.cpus} ${options.args}
        """
}