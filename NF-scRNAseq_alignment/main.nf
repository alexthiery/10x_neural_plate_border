#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { gtf_tag_chroms } from '../modules/tools/genome_tools/main.nf' addParams(  options: modules['gtf_tag_chroms'] )

//  include whole alignment workflow
include { scRNAseq_alignment } from '../modules/workflows/scRNAseq_alignment/main.nf' addParams(cellranger_mkgtf_options: modules['cellranger_mkgtf'],
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
    // Add prefix to chromosomes of interest
    gtf_tag_chroms( ch_gtf )

    scRNAseq_alignment( ch_fasta, gtf_tag_chroms.out, params.input )
}