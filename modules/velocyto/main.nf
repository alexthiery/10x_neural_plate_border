#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process velocyto_run_10x {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container "alexthiery/10x-modules-velocyto:latest"

    input:
        val opts
        tuple val(meta), path(reads)
        path gtf

    output:
        // tuple val(meta), path("*[loomhdf5]"), emit: velocyto

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        velocyto_command = "velocyto run10x ${args} --samtools-threads ${task.cpus} --samtools-threads ${task.memory} ${reads} ${gtf}"
        if (params.verbose){
            println ("[MODULE] velocyto/run_10x command: " + velocyto_command)
        }

        """
        ${velocyto_command}
        """
}