#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_LOOM {
    tag "$meta.sample_id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    container "alexthiery/10x-npb-scvelo:base-1.4"

    input:
        tuple val(meta), path(loom_dir)

    output:
        tuple val(meta), path("*.loom"), emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "merged"
        """
        $moduleDir/bin/merge_loom.py --input ${loom_dir} --output ${prefix}.loom
        """
}