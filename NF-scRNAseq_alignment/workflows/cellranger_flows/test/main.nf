#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//  include whole alignment workflow
include { scRNAseq_alignment_cellranger } from '../main.nf' addParams(  cellranger_mkgtf_options: modules['cellranger_mkgtf'],
                                                                        cellranger_mkref_options: modules['cellranger_mkref'],
                                                                        cellranger_count_options: modules['cellranger_count'])

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

workflow {
    scRNAseq_alignment_cellranger( ch_fasta, ch_gtf, params.samplesheet )
}