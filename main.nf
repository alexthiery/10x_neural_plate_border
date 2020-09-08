#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {cellranger_filter_gtf; cellranger_mkref; cellranger_count} from "$baseDir/modules/cellranger/main.nf"

params.gtf = "$baseDir/test_data/genome/chr1.gtf"
params.fasta = "$baseDir/test_data/genome/chr1.fa"
params.metadata = "$baseDir/sampleInfo.csv"

Channel
    .from(params.gtf)
    .set {ch_gtf}

Channel
    .from(params.fasta)
    .set {ch_fasta}

Channel
    .fromPath( params.metadata )
    .splitCsv(header: ['sample_id', 'sample_name', 'dir1'], skip: 1 )
    .map { row -> [row.sample_id, row.sample_name, file(row.dir1)] }
    .set { ch_fastq }

workflow {
    cellranger_filter_gtf(params.modules['cellranger_filter_gtf'], ch_gtf)
    cellranger_filter_gtf.out.view()

    cellranger_mkref(params.modules['cellranger_mkref'], cellranger_filter_gtf.out, ch_fasta)
    cellranger_count( ch_fastq.combine(cellranger_mkref.out) )
}
