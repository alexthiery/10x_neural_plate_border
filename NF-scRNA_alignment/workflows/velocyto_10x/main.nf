#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {velocyto_run_10x; velocyto_samtools} from "$baseDir/modules/velocyto/main.nf"

workflow velocyto_10x {
    take:
        cellranger_out
        gtf

    main:
        velocyto_samtools( params.modules['velocyto_samtools'], cellranger_out )
        velocyto_run_10x( params.modules['velocyto_run_10x'], velocyto_samtools.out.sortedCellrangerOut, gtf )
    
    emit:
        velocyto_counts = velocyto_run_10x.out.velocytoCounts
}
