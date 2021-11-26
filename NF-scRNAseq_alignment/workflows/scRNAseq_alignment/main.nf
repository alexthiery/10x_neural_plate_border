#!/usr/bin/env nextflow

params.cellranger_mkgtf_options    = [:]
params.cellranger_mkref_options    = [:]
params.cellranger_count_options    = [:]
params.velocyto_samtools_options   = [:]
params.velocyto_run_10x_options    = [:]

//  include cellranger and velocyto subworkflows
include { SCRNASEQ_ALIGNMENT_CELLRANGER } from '../cellranger_flows/main.nf' addParams(  cellranger_mkgtf_options: params.cellranger_mkgtf_options,
                                                                                            cellranger_mkref_options: params.cellranger_mkref_options,
                                                                                            cellranger_count_options: params.cellranger_count_options )

include { VELOCYTO_CELLRANGER } from '../velocyto_flows/main.nf'             addParams(  velocyto_samtools_options: params.velocyto_samtools_options,
                                                                                            velocyto_run_10x_options: params.velocyto_run_10x_options )


workflow SCRNASEQ_ALIGNMENT {
    take:
        fasta
        gtf
        samplesheet

    main:

        // Run cellranger alignment
        SCRNASEQ_ALIGNMENT_CELLRANGER( fasta, gtf, samplesheet )

        // Run RNA velocity on cellranger output
        VELOCYTO_CELLRANGER( gtf, SCRNASEQ_ALIGNMENT_CELLRANGER.out.cellranger_out )
}