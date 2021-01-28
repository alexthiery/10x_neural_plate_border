#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {modify_gtf} from "$baseDir/../custom-nf-modules/tools/modify_gtf/main.nf"
include {cellranger_alignment} from "$baseDir/../custom-nf-modules/workflows/cellranger_flows/main.nf"
include {velocyto_cellranger} from "$baseDir/../custom-nf-modules/workflows/velocyto_flows/main.nf"


Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.genome)
    .set {ch_genome}
    
workflow {
    // tag MT, W and Z genes in GTF
    modify_gtf( params.modules['modify_gtf'], ch_gtf )

    // run cellranger
    cellranger_alignment( ch_genome, modify_gtf.out.gtf, params.sample_csv )

    // run velocyto on cellranger output
    velocyto_cellranger( cellranger_alignment.out.cellranger_out, modify_gtf.out.gtf )
}