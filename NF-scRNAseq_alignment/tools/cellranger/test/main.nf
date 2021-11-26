#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {tenx_fastq_metadata} from "$baseDir/../../../luslab-nf-modules/tools/metadata/main.nf" 
include {cellranger_mkgtf} from "$baseDir/../main.nf" addParams(options: modules['cellranger_mkgtf'])
include {cellranger_mkref} from "$baseDir/../main.nf" addParams(options: modules['cellranger_mkref'])
include {cellranger_count} from "$baseDir/../main.nf" addParams(options: modules['cellranger_count'])

Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.fasta)
    .set {ch_fasta}
    
workflow {
    tenx_fastq_metadata(params.samplesheet)
    cellranger_mkgtf(ch_gtf)
    cellranger_mkref(cellranger_mkgtf.out, ch_fasta)
    cellranger_count(tenx_fastq_metadata.out, cellranger_mkref.out.collect() )
}