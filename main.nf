#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {cellranger_filter_gtf} from "$baseDir/modules/cellranger/main.nf"

params.gtf = "$baseDir/test_data/genome/chr1.gtf"


Channel
    .from(params.gtf)
    .set {ch_gtf}


workflow {
    cellranger_filter_gtf(params.modules['cellranger_filter_gtf'], ch_gtf)
    cellranger_filter_gtf.out.view()
}

