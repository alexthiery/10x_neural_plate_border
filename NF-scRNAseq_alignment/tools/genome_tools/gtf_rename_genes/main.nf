#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GTF_RENAME_GENES {
    
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
        prefix = options.suffix ? "${options.suffix}" : "rename_genes"

    """
    #!/usr/local/bin/python
    
    import re
    
    def rename_genes(gtf, outfile):
        gtf =  open(gtf, 'rt')
        
        # create output file to write to
        outfile = open(outfile, 'a')

        new_gtf = []

        for line in gtf:
            if 'gene_id' in line:
        #       Check if any gene_ids in line
                if any(gene_id in line for gene_id in goi.keys()):
        #           Get matching gene name from dict
                    gene_name = [gene_name for gene_id, gene_name in goi.items() if gene_id in line]
        #           If gene_name is already in line - remove and replace
                    if 'gene_name' in line:
                        line = re.sub('gene_name.*?;', '', line, flags=re.DOTALL)
                    line = line.rstrip() + ' gene_name "' + gene_name[0] + '";\n'
            outfile.write(line)

    goi = {'ENSGALG00000030902': 'SNAI2'}
    rename_genes(gtf = ${gtf}, 'rt', outfile = ${prefix}.gtf)
    """
}


