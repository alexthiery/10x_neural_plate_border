/*
 * Split Seurat object by batch and run downstream analysis
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.split                              = file("$baseDir/bin/seurat/split_seurat.R", checkIfExists: true)
analysis_scripts.cluster                            = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)
analysis_scripts.gene_modules                       = file("$baseDir/bin/other/subset_gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification               = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)
analysis_scripts.gene_modules_latent_time           = file("$baseDir/bin/other/gene_modules_latent_time.R", checkIfExists: true)

params.split_options                                = [:]
params.cluster_options                              = [:]
params.gene_modules_options                         = [:]
params.state_classification_options                 = [:]
params.seurat_h5ad_options                          = [:]
params.seurat_intersect_loom_options                = [:]
params.scvelo_run_options                           = [:]
params.gene_modules_latent_time_options             = [:]



// Include Seurat R processes
include {R as SPLIT} from "$baseDir/modules/local/r/main"                           addParams(  options:                        params.split_options,
                                                                                                script:                         analysis_scripts.split )

include {R as CLUSTER} from "$baseDir/modules/local/r/main"                         addParams(  options:                        params.cluster_options,
                                                                                                script:                         analysis_scripts.cluster )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                    addParams(  options:                        params.gene_modules_options,
                                                                                                script:                         analysis_scripts.gene_modules )

include {R as STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"            addParams(  options:                        params.state_classification_options,
                                                                                                script:                         analysis_scripts.state_classification )

// include {SEURAT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"                addParams(  options:                        params.seurat_h5ad_options )

// include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"             addParams(  seurat_intersect_loom_options:  params.seurat_intersect_loom_options,
//                                                                                                 scvelo_run_options:             params.scvelo_run_options)

// include {R as GENE_MODULES_LATENT_TIME} from "$baseDir/modules/local/r/main"        addParams(  options:                        params.gene_modules_latent_time_options,
//                                                                                                 script:                         analysis_scripts.gene_modules_latent_time )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_RUN_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]
    // loom            //Channel: merged.loom
    // annotations     //Channel: seurat_annotations.csv

    main:
    // Run Seurat pipeline
    SPLIT( seurat_out )

    SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    CLUSTER( ch_split_run )
    GENE_MODULES( CLUSTER.out )
    // STATE_CLASSIFICATION( CLUSTER.out )

    // SEURAT_H5AD( CLUSTER.out )

    // // Run scVelo
    // SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv

    // // Run gene module analysis across latent time
    // ch_cluster_rds              = CLUSTER.out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    // ch_gene_modules_rds         = GENE_MODULES.out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    // ch_gene_module_latent_time  = ch_cluster_rds.combine(ch_gene_modules_rds, by: 0).combine(SEURAT_SCVELO.out.scvelo_run_out_metadata, by: 0)
    // ch_gene_module_latent_time  = ch_gene_module_latent_time.map{[it[0], [it[1], it[2], it[3]]]}
    
    // GENE_MODULES_LATENT_TIME(ch_gene_module_latent_time)

    emit:
    cluster_out                     = CLUSTER.out                               //Channel: [[meta], [output]]
    gene_modules_out                = GENE_MODULES.out                          //Channel: [[meta], [output]]
    // state_classification_out        = STATE_CLASSIFICATION.out                  //Channel: [[meta], [output]]   

    // scvelo_run_out_metadata         = SEURAT_SCVELO.out.scvelo_run_out_metadata     //Channel: [[meta], csv]
    // scvelo_run_out_5had             = SEURAT_SCVELO.out .scvelo_run_out_h5ad        //Channel: [[meta], h5ad]
}
