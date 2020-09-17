#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process modifyGTF {

    //publishDir "${params.alignment_outDir}/modifiedGTF",
    //mode: "copy", overwrite: true

    //publishDir "test/modifiedGTF",
    //mode: "copy", overwrite: true

    input:
        path(pythonFile)
        path(gtf)

    output:
        path('./modified.gtf')

    """
    #!/bin/bash
    python ${pythonFile} ${gtf} ./modified.gtf

    """
}

params.python = "$baseDir/modifyGTF.py"
params.gtf = "$baseDir/test.gtf"


Channel
    .from(params.python)
    .set {ch_python}

Channel
    .from(params.gtf)
    .set {ch_gtf}
    
workflow {
    modifyGTF(ch_python, ch_gtf)
}