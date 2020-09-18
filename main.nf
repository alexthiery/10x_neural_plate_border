#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include {modify_gtf} from "$baseDir/modules/modify_gtf/main.nf"
include {cellranger_alignment} from "$baseDir/workflows/scRNAseq_alignment/main.nf"
include {velocyto_run_10x} from "$baseDir/modules/velocyto/main.nf"

Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.genome)
    .set {ch_genome}
    
workflow {
    // tag MT, W and Z genes in GTF
    // modify_gtf( ch_gtf )

    // run cellranger
    cellranger_alignment( ch_gtf, ch_genome, params.sample_csv )

    // run velocyto on cellranger output
    velocyto_run_10x( params.modules['velocyto_run_10x'], cellranger_alignment.out.cellranger_out, ch_gtf )
}