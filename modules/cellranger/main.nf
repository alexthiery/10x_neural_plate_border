#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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
        tuple val(meta), path("${meta.sample_name}_filtered_feature_bc_matrix"), emit: readCounts
        tuple val(meta), path("cellrangerOut_${meta.sample_name}"), emit: cellrangerOut

    script:
        args = ""
            if(opts.args && opts.args != '') {
                ext_args = opts.args
                args += ext_args.trim()
            }

        cellranger_count_command = "cellranger count --id='cellrangerOut_${meta.sample_name}' --fastqs='./' --sample=${meta.sample_id} --transcriptome=${reference_genome} ${args}"
        
        // Log
        if (params.verbose){
            println ("[MODULE] cellranger count command: " + cellranger_count_command)
        }

       //SHELL
        """
        ${cellranger_count_command}
        mkdir ${meta.sample_name}_filtered_feature_bc_matrix
        cp cellrangerOut_${meta.sample_name}/outs/filtered_feature_bc_matrix/*.gz ${meta.sample_name}_filtered_feature_bc_matrix
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
            println ("[MODULE] filter gtf command: " + filter_gtf_command)
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
