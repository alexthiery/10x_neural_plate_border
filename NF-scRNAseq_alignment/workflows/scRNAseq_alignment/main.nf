#!/usr/bin/env nextflow

params.cellranger_mkgtf_options    = [:]
params.cellranger_mkref_options    = [:]
params.cellranger_count_options    = [:]
params.velocyto_samtools_options   = [:]
params.velocyto_run_10x_options    = [:]

//  include cellranger and velocyto subworkflows
include { scRNAseq_alignment_cellranger } from '../cellranger_flows/main.nf' addParams(  cellranger_mkgtf_options: params.cellranger_mkgtf_options,
                                                                                            cellranger_mkref_options: params.cellranger_mkref_options,
                                                                                            cellranger_count_options: params.cellranger_count_options )

include { velocyto_cellranger } from '../velocyto_flows/main.nf'             addParams(  velocyto_samtools_options: params.velocyto_samtools_options,
                                                                                            velocyto_run_10x_options: params.velocyto_run_10x_options )


workflow scRNAseq_alignment {
    take:
        fasta
        gtf
        samplesheet

    main:

        // Run cellranger alignment
        scRNAseq_alignment_cellranger( fasta, gtf, samplesheet )

        // Run RNA velocity on cellranger output
        velocyto_cellranger( gtf, scRNAseq_alignment_cellranger.out.cellranger_out )
}