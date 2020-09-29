#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include {cellranger_alignment} from "$baseDir/workflows/scRNAseq_alignment/main.nf"
// include {velocyto_10x} from "$baseDir/workflows/velocyto_10x/main.nf"
// include {modify_gtf} from "$baseDir/modules/modify_gtf/main.nf"

// Channel
//     .from(params.gtf)
//     .set {ch_gtf}

// Channel
//     .from(params.genome)
//     .set {ch_genome}
    
// workflow {
//     // tag MT, W and Z genes in GTF
//     modify_gtf( params.modules['modify_GTF'], ch_gtf )

//     // run cellranger
//     cellranger_alignment( modify_gtf.out.GTF, ch_genome, params.sample_csv )

//     // run velocyto on cellranger output
//     velocyto_10x( cellranger_alignment.out.cellranger_out, modify_gtf.out.GTF )
// }

include {r_analysis as r_analysis_seurat} from "$baseDir/modules/r_analysis/main.nf"

Channel
    .fromPath("$baseDir/alignment_out/*matrix", type: 'dir')
    .collect()
    .set { ch_read_counts }


// Channel
//     .fromPath(params.seurat_rds, type: 'dir')
//     .collect()
//     .set { ch_read_counts }


workflow {
    r_analysis_seurat( params.modules['seurat_full'], ch_read_counts )
    // r_analysis_antler( params.modules['antler'], ch_read_counts )
}