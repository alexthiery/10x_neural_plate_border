#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEURAT_H5AD {

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "alexthiery/10x-npb-scvelo:base-1.5"

    input:
        tuple val(meta), path('input/*')

    output:
        tuple val(meta), file('*h5ad')

    script:

        """
        Rscript $moduleDir/bin/seurat_h5ad.R ${options.args}
        """
}