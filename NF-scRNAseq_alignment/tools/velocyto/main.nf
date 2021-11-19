#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process velocyto_run_10x {

    label "process_high"

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    
    container "quay.io/biocontainers/velocyto.py:0.17.17--py37h97743b1_2"

    input:
        tuple val(meta), path(reads)
        path gtf

    output:
        tuple val(meta), path("*.loom"), emit: velocyto_counts

    script:
        prefix = meta.run ? "${meta.sample_name}_${meta.run}" : "${meta.sample_name}"

        velocyto_command = "velocyto run10x ${options.args} ${reads} ${gtf}"

        if (params.verbose){
            println ("[MODULE] velocyto/run_10x command: " + velocyto_command)
        }

        """
        ${velocyto_command}
        mv ${reads}/velocyto/${reads}.loom ./${prefix}.loom
        """
}

process velocyto_samtools {

    label "process_medium"

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    
    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(reads), emit: sorted_cellranger_out

    script:
        velocyto_samtools_command = "samtools sort -t CB -O BAM -@ ${task.cpus} -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam"

        if (params.verbose){
            println ("[MODULE] velocyto/samtools command: " + velocyto_samtools_command)
        }

        """
        cd ${reads}/outs/
        ${velocyto_samtools_command}
        """
}
