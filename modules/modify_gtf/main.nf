#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process modify_gtf {

    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    container "alexthiery/10x-modules-modify_gtf:latest"

    input:
        val opts
        path gtf

    output:
        path "${opts.filename}.gtf", emit: GTF

    """
    #!/opt/conda/bin/python
    
    import pandas as pd
    import re
  
    # make list of chromosomes to edit
    lab = [${opts.edit_chroms}]

    # create output file to write to
    outfile = open("./modified.gtf", 'a')

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


