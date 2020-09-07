#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {cellranger_filter_gtf; cellranger_mkref} from "$baseDir/modules/cellranger/main.nf"

params.gtf = "$baseDir/test_data/genome/chr1.gtf"
params.fasta = "$baseDir/test_data/genome/Gallus_gallus.sub.fa"


Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.fasta)
    .set {ch_fasta}

workflow {
    cellranger_filter_gtf(params.modules['cellranger_filter_gtf'], ch_gtf)
    cellranger_filter_gtf.out.view()

    cellranger_mkref(params.modules['cellranger_mkref'], cellranger_filter_gtf.out, ch_fasta)
}

