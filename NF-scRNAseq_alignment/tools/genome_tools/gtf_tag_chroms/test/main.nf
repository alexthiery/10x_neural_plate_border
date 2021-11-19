#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { gtf_tag_chroms } from '../main.nf' addParams(options: modules['gtf_tag_chroms'] )

params.gtf = "$baseDir/../../../test_data/genome_tools/test.gtf"

Channel
    .from(params.gtf)
    .set {ch_gtf}
    
workflow {
    gtf_tag_chroms(ch_gtf)
}