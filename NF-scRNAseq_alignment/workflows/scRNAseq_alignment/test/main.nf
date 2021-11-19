#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//  include whole alignment workflow
include { scRNAseq_alignment } from '../main.nf' addParams(  gtf_tag_chroms_options: modules['gtf_tag_chroms'],
                                                                                cellranger_mkgtf_options: modules['cellranger_mkgtf'],
                                                                                cellranger_mkref_options: modules['cellranger_mkref'],
                                                                                cellranger_count_options: modules['cellranger_count'],
                                                                                velocyto_samtools_options: modules['velocyto_samtools'],
                                                                                velocyto_run_10x_options: modules['velocyto_run_10x'] )

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

workflow {
    scRNAseq_alignment( ch_fasta, ch_gtf, params.samplesheet )
}