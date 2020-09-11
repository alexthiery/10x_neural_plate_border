#!/usr/bin/env nextflow

nextflow.preview.dsl=2


process cellranger_count {

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    container "alexthiery/10x-modules-cellranger:latest"

    input:
        val opts
        tuple val(meta), path(reads)
        path reference_genome

    output:
        val "${meta.sample_name}", emit: sampleName
        path "${meta.sample_name}", emit: countFiles

    script:
    args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

    """
    #!/bin/bash

    cellranger count --id="cellrangerOut_${sample_name}" \
    --fastqs="dir1/${sample_id}"\
    --sample=${sample_id} \
    --transcriptome=${reference_genome} \
     ${args}
    

    mkdir ${meta.sample_name}

    mv cellrangerOut_${meta.sample_name}/outs/filtered_feature_bc_matrix/*.gz ${meta.sample_name}
    """
}


process cellranger_filter_gtf {

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }


    container "alexthiery/10x-modules-cellranger:latest"

    input:
        val(opts)
        path(gtf)

    output:
        path("${opts.outfile}")

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        filter_gtf_command = "cellranger mkgtf ${gtf} ${opts.outfile} ${args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] filter_gtf command: " + filter_gtf_command)
        }

    //SHELL
        """
        ${filter_gtf_command}
        """

   
}

process cellranger_mkref {

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }


    container "alexthiery/10x-modules-cellranger:latest"

    input:
        val(opts)
        path(filt_genome)
        path(fasta)

    output:
        path("reference_genome")

    script:
        mkref_command = "cellranger mkref --genome=reference_genome --genes=${filt_genome} --fasta=${fasta} --nthreads=${task.cpus}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] mkref command: " + mkref_command)
        }

    //SHELL
        """
        ${mkref_command}
        """

}
