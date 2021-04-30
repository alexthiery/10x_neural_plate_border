#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCVELO {

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "alexthiery/10x-npb-scvelo:latest"

    input:
        path(loom)

    output:
        path "${prefix}.loom", emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? ${options.prefix} : "seurat_merged"
        """
        scvelo.py --input ${loom} --output ${prefix}.loom
        """
}