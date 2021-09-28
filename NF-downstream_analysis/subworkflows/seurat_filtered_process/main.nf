/*
 * Downstream exploratory analysis on full Seurat filtered dataset
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream exploratory analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.scatterplot3d                      = file("$baseDir/bin/other/scatterplot3d.R", checkIfExists: true)
analysis_scripts.gene_modules                       = file("$baseDir/bin/other/gene_modules.R", checkIfExists: true)
analysis_scripts.state_classification               = file("$baseDir/bin/seurat/state_classification.R", checkIfExists: true)
analysis_scripts.gene_modules_latent_time           = file("$baseDir/bin/other/gene_modules_latent_time.R", checkIfExists: true)
analysis_scripts.phate                              = file("$baseDir/bin/other/phateR.R", checkIfExists: true)

params.scatterplot3d_options                        = [:]
params.gene_module_options                          = [:]
params.state_classification_options                 = [:]
params.phate_options                                = [:]
params.seurat_h5ad_options                          = [:]
params.seurat_intersect_loom_options                = [:]
params.scvelo_run_options                           = [:]
params.cellrank_run_options                         = [:]

// Include R processes

include {R as SCATTERPLOT3D} from "$baseDir/modules/local/r/main"                   addParams(  options:                        params.scatterplot3d_options,
                                                                                                script:                         analysis_scripts.scatterplot3d )

include {R as GENE_MODULES} from "$baseDir/modules/local/r/main"                    addParams(  options:                        params.gene_module_options,
                                                                                                script:                         analysis_scripts.gene_modules )

include {R as STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"            addParams(  options:                        params.state_classification_options,
                                                                                                script:                         analysis_scripts.state_classification )

include {R as PHATE} from "$baseDir/modules/local/r/main"                           addParams(  options:                        params.phate_options,
                                                                                                script:                         analysis_scripts.phate )

include {SEURAT_H5AD} from "$baseDir/modules/local/seurat_h5ad/main"                addParams(  options:                        params.seurat_h5ad_options)

include {SEURAT_SCVELO} from "$baseDir/subworkflows/seurat_scvelo/main"             addParams(  seurat_intersect_loom_options:  params.seurat_intersect_loom_options,
                                                                                                scvelo_run_options:             params.scvelo_run_options )

include {CELLRANK_RUN} from "$baseDir/modules/local/scvelo/cellrank_run/main"       addParams(  options:                        params.cellrank_run_options )

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow SEURAT_FILTERED_PROCESS {
    take:
    seurat_out      //Channel: [[meta], [plot_dir, rds_dir]]
    loom            //Channel: merged.loom
    annotations     //Channel: seurat_annotations.csv

    main:
    // Run processes on full filtered dataset
    SCATTERPLOT3D(seurat_out)
    STATE_CLASSIFICATION( seurat_out )
    GENE_MODULES( STATE_CLASSIFICATION.out )
    PHATE( STATE_CLASSIFICATION.out )

    // // Run scVelo
    SEURAT_H5AD( STATE_CLASSIFICATION.out )
    SEURAT_SCVELO( SEURAT_H5AD.out, loom, annotations ) // Channel: [[meta], seurat.h5ad], Channel: merged.loom, Channel: seurat_annotations.csv
    CELLRANK_RUN( SEURAT_SCVELO.out.scvelo_run_out_h5ad )

    // // Run gene module analysis across latent time
    // ch_cluster_rds              = CLUSTER.out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    // ch_gene_modules_rds         = GENE_MODULES.out.map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //Channel: [[meta], *.rds_file]
    // ch_gene_module_latent_time  = ch_cluster_rds.combine(ch_gene_modules_rds, by: 0).combine(SEURAT_SCVELO.out.scvelo_run_out_metadata, by: 0)
    // ch_gene_module_latent_time  = ch_gene_module_latent_time.map{[it[0], [it[1], it[2], it[3]]]}
    
    // GENE_MODULES_LATENT_TIME(ch_gene_module_latent_time)

    emit:
    state_classification_out        = STATE_CLASSIFICATION.out                  //Channel: [[meta], [output]]
    gene_modules_out                = GENE_MODULES.out                          //Channel: [[meta], [output]]

    cellrank_run_out_metadata       = CELLRANK_RUN.out.scvelo_run_out_metadata     //Channel: [[meta], csv]
}

