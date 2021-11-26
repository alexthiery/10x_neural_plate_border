#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {GTF_TAG_CHROMS} from './tools/genome_tools/gtf_tag_chroms/main.nf'     addParams(  options: modules['gtf_tag_chroms'] )

include {GTF_RENAME_GENES} from './tools/genome_tools/gtf_rename_genes/main.nf' addParams(  options: modules['gtf_rename_genes'] )

//  include whole alignment workflow
include {SCRNASEQ_ALIGNMENT} from './workflows/scRNAseq_alignment/main.nf'      addParams(  cellranger_mkgtf_options: modules['cellranger_mkgtf'],
                                                                                            cellranger_mkref_options: modules['cellranger_mkref'],
                                                                                            cellranger_count_options: modules['cellranger_count'],
                                                                                            velocyto_samtools_options: modules['velocyto_samtools'],
                                                                                            velocyto_run_10x_options: modules['velocyto_run_10x'] )

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}


workflow {
    // Add prefix to chromosomes of interest
    GTF_TAG_CHROMS( ch_gtf )

    GTF_RENAME_GENES(GTF_TAG_CHROMS.out)

    SCRNASEQ_ALIGNMENT( ch_fasta, GTF_RENAME_GENES.out, params.input )
}