#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process gtf_tag_chroms {
    
    label 'process_low'

    publishDir "${params.outdir}",
        mode: 'copy',
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/bioframe:0.0.12--pyh3252c3a_0"

    input:
        path gtf

    output:
        path "${prefix}.gtf", emit: gtf

    script:
        prefix = options.suffix ? "${options.suffix}" : "tag_chroms"

    """
    #!/usr/local/bin/python
    
    import pandas as pd
    import re
  
    # make array of chromosomes to tag
    lab = [${options.args}]

    # create output file to write to
    outfile = open("${prefix}.gtf", 'a')

    with open("${gtf}", 'rt') as gtf:
        for line in gtf:
            # only search lines with gene_id in order to skip header lines
            if 'gene_id' in line:

                # change names for genes from chroms of interest (lab)
                chr = line.split()[0]
                if chr in lab:
                    match = [s for s in lab if chr in s]
                    line = re.sub('gene_id "', 'gene_id "'+''.join(match)+'-', line)
                    if 'gene_name' in line:
                        line = re.sub('gene_name "', 'gene_name "'+''.join(match)+'-', line)

            outfile.write(line)
    """
}


