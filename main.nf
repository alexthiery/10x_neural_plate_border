#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {tenx_fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
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
    tenx_fastq_metadata( params.sample_csv )
    cellranger_alignment( ch_gtf, ch_genome, tenx_fastq_metadata.out.metadata )

    // run velocyto on cellranger output
    velocyto_run_10x( params.modules['velocyto_run_10x'], cellranger_alignment.out.cellranger_out, ch_gtf )
}