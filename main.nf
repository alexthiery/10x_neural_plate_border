#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include {cellranger_alignment} from "$baseDir/workflows/scRNAseq_alignment/main.nf"
include {velocyto_run_10x; velocyto_samtools} from "$baseDir/modules/velocyto/main.nf"
include {modify_gtf} from "$baseDir/modules/modify_gtf/main.nf"

Channel
    .from(params.gtf)
    .set {ch_gtf}

// Channel
//     .from(params.genome)
//     .set {ch_genome}
    
// workflow {
//     // tag MT, W and Z genes in GTF
//     modify_gtf( params.modules['modify_GTF'], ch_gtf )

//     // run cellranger
//     cellranger_alignment( modify_gtf.out.GTF, ch_genome, params.sample_csv )

//     // run velocyto on cellranger output
//     velocyto_samtools( params.modules['velocyto_samtools'], cellranger_alignment.out.cellranger_out )
//     velocyto_run_10x( params.modules['velocyto_run_10x'], velocyto_samtools.out.sorted_cellrangerOut, modify_gtf.out.GTF )
// }

homerData = [
    [[sample_id:'S1'], "/camp/home/thierya/scratch/10x_neural_plate_border/work/f0/0291b794d2877e1a550a016d19c3ae/cellrangerOut_hh4"]
]

Channel
    .from(homerData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_homerData}


workflow {
    // tag MT, W and Z genes in GTF
    modify_gtf( params.modules['modify_GTF'], ch_gtf )

    // run cellranger
    // cellranger_alignment( modify_gtf.out.GTF, ch_genome, params.sample_csv )

    // run velocyto on cellranger output
    // velocyto_samtools( params.modules['velocyto_samtools'], cellranger_alignment.out.cellranger_out )
    velocyto_run_10x( params.modules['velocyto_run_10x'], ch_homerData, modify_gtf.out.GTF )
}