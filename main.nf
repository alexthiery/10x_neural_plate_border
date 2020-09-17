#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {cellranger_alignment} from "$baseDir/workflows/scRNAseq_alignment/main.nf"

Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.genome)
    .set {ch_genome}
    
workflow {
    cellranger_alignment( ch_gtf, ch_genome, params.sample_csv )
}