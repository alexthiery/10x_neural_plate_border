#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {modifyGTF} from "$baseDir/../main.nf"
params.gtf = "$baseDir/../../../test_data/modify_gtf/test.gtf"

Channel
    .from(params.gtf)
    .set {ch_gtf}
    
workflow {
    modifyGTF(params.modules['modify_GTF'], ch_gtf)
}