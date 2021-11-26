#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {metadata} from "../../metadata/main.nf"

include {r as r_1} from '../main.nf' addParams(options: modules['r_1'], script: modules['r_1'].script)
include {r as r_2} from '../main.nf' addParams(options: modules['r_2'], script: modules['r_2'].script)


params.input = "$baseDir/samplesheet.csv"

workflow {
    metadata(params.input)
    r_1 (metadata.out.filter{ it[0].sample_id == 'test' } )
    r_1.out.view()
    r_2 (r_1.out)
}