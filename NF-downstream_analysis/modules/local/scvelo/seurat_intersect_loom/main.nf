#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEURAT_INTERSECT_LOOM {

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "alexthiery/10x-npb-scvelo:latest"

    input:
        tuple val(meta), path(loom), path(seurat), path(annotations)

    output:
        tuple val(meta), path("*.loom"), emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "seurat_merged"
        
        """
        $moduleDir/bin/seurat_intersect_loom.py --loomInput ${loom} --seuratInput ${seurat} --annotations ${annotations} --output ${prefix}.loom
        """
}