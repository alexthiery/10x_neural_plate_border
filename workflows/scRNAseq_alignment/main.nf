#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {tenx_fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
include {cellranger_filter_gtf; cellranger_mkref; cellranger_count} from "$baseDir/modules/cellranger/main.nf"

workflow cellranger_alignment {
    take:
        gtf
        genome
        sample_csv

    main:
        tenx_fastq_metadata(sample_csv)
        cellranger_filter_gtf(params.modules['cellranger_filter_gtf'], gtf)
        cellranger_mkref(params.modules['cellranger_mkref'], cellranger_filter_gtf.out, genome)
        cellranger_count(params.modules['cellranger_count'], tenx_fastq_metadata.out, cellranger_mkref.out.collect() )
    
    emit:
        cellranger_count.out.readCounts
        cellranger_count.out.cellrangerOut
}
