#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEURAT_H5AD {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    container "alexthiery/10x-npb-schelper:latest"

    input:
        tuple val(meta), path('input/*')

    output:
        tuple val(meta), file('*h5ad')

    script:

        """
        Rscript $moduleDir/bin/seurat_h5ad.R ${options.args}
        """
}
