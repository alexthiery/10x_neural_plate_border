/*
 * Subset clusters from filtered Seurat object and run downstream analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.subset                             = file("$baseDir/bin/seurat/subset_cells.R", checkIfExists: true)
analysis_scripts.cluster                            = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.gene_modules                       = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification               = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)
analysis_scripts.phate                              = file("$baseDir/bin/other/phateR.R", checkIfExists: true)

// Include Seurat R processes
include {R as SUBSET} from "$baseDir/modules/local/r/main"                          addParams( script:                         analysis_scripts.subset )
include {R as CLUSTER} from "$baseDir/modules/local/r/main"                         addParams( script:                         analysis_scripts.cluster )
include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                    addParams( script:                         analysis_scripts.gene_modules )

include {SEURAT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"
include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"


/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_SUBSET_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]
    loom            //Channel: merged.loom
    annotations     //Channel: seurat_annotations.csv

    main:
    // Run Seurat pipeline
    SUBSET( seurat_out )

    SUBSET.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    CLUSTER( ch_split_run )
    GENE_MODULES( CLUSTER.out )

    // Run scVelo
    SEURAT_H5AD( CLUSTER.out )
    SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv

    emit:
    cluster_out                     = CLUSTER.out                               //Channel: [[meta], [output]]
    gene_modules_out                = GENE_MODULES.out                          //Channel: [[meta], [output]]
}

