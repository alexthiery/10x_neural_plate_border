#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GTF_RENAME_GENES {
    
    label 'process_low'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/bioframe:0.0.12--pyh3252c3a_0"

    input:
        path gtf

    output:
        path "${prefix}.gtf", emit: gtf

    script:
        prefix = options.suffix ? "${options.suffix}" : "rename_genes"

    """
    gtf_rename_genes.py --input ${gtf} --output ${prefix}.gtf
    """
}


