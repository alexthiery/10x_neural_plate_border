#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEURAT_INTERSECT_LOOM {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    container "alexthiery/10x-npb-scvelo:base-1.7"

    input:
        tuple val(meta), path('input/'seurat), path('input/'loom), path('input/'annotations)

    output:
        tuple val(meta), path("*.loom"), emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "seurat_merged"
        
        """
        $moduleDir/bin/seurat_intersect_loom.py --loomInput input/${loom} --seuratInput input/${seurat} --annotations input/${annotations} --output ${prefix}.loom
        """
}