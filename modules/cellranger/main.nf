#!/usr/bin/env nextflow

nextflow.preview.dsl=2


process cellranger_count {

    publishDir "${params.alignment_outDir}/cellrangerCounts",
        mode: "copy", overwrite: true

    container "alexthiery/10x_mod_cellranger:latest"

    input:
        tuple val(sample_id), val(sample_name), path('dir1/*'), path(reference_genome)

    output:
        val sample_name, emit: sampleName
        path sample_name, emit: countFiles

    """
    #!/bin/bash
    
    cellranger count --id="cellrangerOut_${sample_name}" \
    --fastqs="dir1/${sample_id}"\
    --sample=${sample_id} \
    --transcriptome=${reference_genome}

    mkdir ${sample_name}

    mv cellrangerOut_${sample_name}/outs/filtered_feature_bc_matrix/*.gz ${sample_name}
    """
}


process cellranger_filter_gtf {

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }


    container "alexthiery/10x_neural_tube:v1.0"

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


    container "alexthiery/10x_neural_tube:v1.0"

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