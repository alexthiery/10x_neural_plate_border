#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process cellranger_count {

    label 'process_high'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
        tuple val(meta), path('fastqs/*')
        path reference_genome

    output:
        tuple val(meta), path("${prefix}_cellranger"), emit: cellranger_out
        tuple val(meta), path("${prefix}/*"), emit: read_counts

    script:
        prefix = meta.run ? "${meta.sample_name}_${meta.run}" : "${meta.sample_name}"

        cellranger_count_command = "cellranger count --id='${prefix}_cellranger' --fastqs='fastqs' --sample=${meta.sample_id} --transcriptome=${reference_genome} ${options.args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] cellranger count command: " + cellranger_count_command)
        }

       //SHELL
        """
        ${cellranger_count_command}
        mkdir ${prefix}
        cp ${prefix}_cellranger/outs/filtered_feature_bc_matrix/* ${prefix}
        """
}


process cellranger_mkgtf {

    label 'process_low'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
        path(gtf)

    output:
        path("*.gtf")

    script:
        
        mkgtf_command = "cellranger mkgtf ${gtf} ${gtf.baseName}_mkgtf.gtf ${options.args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] filter gtf command: " + mkgtf_command)
        }

       //SHELL
        """
        ${mkgtf_command}
        """
}

process cellranger_mkref {

    label 'process_medium'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "streitlab/custom-nf-modules-cellranger:latest"

    input:
        path(gtf)
        path(fasta)

    output:
        path("reference_genome")

    script:
        mkref_command = "cellranger mkref --genome=reference_genome --genes=${gtf} --fasta=${fasta} --nthreads=${task.cpus}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] mkref command: " + mkref_command)
        }

        //SHELL
        """
        ${mkref_command}
        """
}
