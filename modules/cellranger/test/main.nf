#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {tenx_fastq_metadata} from "$baseDir/../../../luslab-nf-modules/tools/metadata/main.nf"
include {cellranger_filter_gtf; cellranger_mkref; cellranger_count} from "$baseDir/../main.nf"

Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.genome)
    .set {ch_genome}
    
workflow {
    tenx_fastq_metadata("$baseDir/../../../test_data/cellranger/sample_info.csv")
    cellranger_filter_gtf(params.modules['cellranger_filter_gtf'], ch_gtf)
    cellranger_mkref(params.modules['cellranger_mkref'], cellranger_filter_gtf.out, ch_genome)
    cellranger_count(params.modules['cellranger_count'], tenx_fastq_metadata.out, cellranger_mkref.out.collect() )
}

