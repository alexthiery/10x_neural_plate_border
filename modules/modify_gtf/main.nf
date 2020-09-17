#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process modifyGTF {

    //publishDir "${params.alignment_outDir}/modifiedGTF",
    //mode: "copy", overwrite: true

    //publishDir "test/modifiedGTF",
    //mode: "copy", overwrite: true

    input:
        path(gtf)

    output:
        path('./modified.gtf')

    """

    #!/Users/evahamrud/opt/anaconda3/bin python
    import pandas as pd
    import re
    import sys
    import os
  
    # make list of chromosomes to edit
    lab = ['MT', 'Z', 'W']

    # create output file to write to
    outfile = open("./modified.gtf", 'a')

    fn = ${params.gtf}
    if os.path.isfile(fn) == False:
    sys.exit("GTF file missing")

    with open(fn, 'rt') as gtf:
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

//params.python = "$baseDir/modifyGTF.py"
params.gtf = "$baseDir/test.gtf"


//Channel
//    .from(params.python)
//    .set {ch_python}

Channel
    .from(params.gtf)
    .set {ch_gtf}
    
workflow {
    modifyGTF(ch_gtf)
}