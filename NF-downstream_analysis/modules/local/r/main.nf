#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process R {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    input:
        tuple val(meta), path('input/*')

    output:
        tuple val(meta), file('*')

    script:

        """
        Rscript ${params.script} --cores ${task.cpus} --runtype nextflow ${options.args}
        rm -r input
        """
}