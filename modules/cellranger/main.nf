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
        println "${meta.sample_name}"

    """
    #!/bin/bash
    
    cellranger count --id="cellrangerOut_${meta.sample_name}" \
    --fastqs=./\
    --sample=${meta.sample_id} \
    --transcriptome=${reference_genome} \
    --chemistry=SC3Pv3

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

    """
    #!/bin/bash

    # this step filters out genes based on the gene biotypes listed in attributes.
    cellranger mkgtf ${gtf} ${opts.outfile} \
    --attribute=gene_biotype:protein_coding,  \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:pseudogene
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

    """
    #!/bin/bash

    # make reference
    cellranger mkref --genome=reference_genome \
    --genes=${filt_genome} \
    --fasta=${fasta} \
    --nthreads=${task.cpus}

    """
}