#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process r_analysis {
    // publishDir "${params.outdir}/${opts.publish_dir}",
    // mode: "copy",
    // overwrite: true,
    // saveAs: { filename ->
    //                 if (opts.publish_results == "none") null
    //                 else filename }

    container "alexthiery/10x-modules-r_analysis:latest"

    input:
        val opts
        path input

    output:
        // path("plots")
        // path("RDS.files")

    script:

    """
    Rscript ${opts.script} --cores ${task.cpus} --custom_functions ${opts.custom_functions} --runtype nextflow 
    """
}