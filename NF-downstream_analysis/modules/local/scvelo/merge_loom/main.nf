#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_LOOM {

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "alexthiery/10x-npb-scvelo:latest"

    input:
        path(loom_dir)

    output:
        path "${prefix}.loom", emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? ${options.prefix} : "merged"
        """
        merge_loom.py --input ${loom_dir} --output ${prefix}.loom
        """
}