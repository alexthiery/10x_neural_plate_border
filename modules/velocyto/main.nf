#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process velocyto_run_10x {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container "quay.io/biocontainers/velocyto.py:0.17.17--py37h97743b1_2"

    input:
        val opts
        tuple val(meta), path(reads)
        path gtf

    output:
        tuple val(meta), path("cellrangerOut_${meta.sample_name}/velocyto/cellrangerOut_${meta.sample_name}.loom"), emit: velocytoCounts

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        velocyto_command = "velocyto run10x ${args} ${reads} ${gtf}"
        if (params.verbose){
            println ("[MODULE] velocyto/run_10x command: " + velocyto_command)
        }

        """
        ${velocyto_command}
        """
}

process velocyto_samtools {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }
    
    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(reads), emit: sortedCellrangerOut

    script:
        velocyto_samtools_command = "samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam"
        if (params.verbose){
            println ("[MODULE] velocyto/samtools command: " + velocyto_samtools_command)
        }

        """
        cd cellrangerOut_${meta.sample_name}/outs/
        ${velocyto_samtools_command}
        """
}