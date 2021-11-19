#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

params.velocyto_samtools_options    = [:]
params.velocyto_run_10x_options     = [:]

include {VELOCYTO_SAMTOOLS} from "../../tools/velocyto/main.nf" addParams(options: params.velocyto_samtools_options)
include {VELOCYTO_RUN_10X} from "../../tools/velocyto/main.nf" addParams(options: params.velocyto_run_10x_options)

// Define workflow to subset and index a genome region fasta file
workflow velocyto_cellranger {
    take:
        gtf
        cellranger_out
        
    main:
        // Sort bam
        VELOCYTO_SAMTOOLS( cellranger_out )

        // Calculate velocity
        VELOCYTO_RUN_10X( VELOCYTO_SAMTOOLS.out.sorted_cellranger_out, gtf )

    emit:
        velocyto_counts = VELOCYTO_RUN_10X.out.velocyto_counts
}