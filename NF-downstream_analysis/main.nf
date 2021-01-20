#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis} from "$baseDir/../modules/r_analysis/main.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
// log.info projectHeader()

// Header log info
// def summary = [:]
// summary['Run Name']               = workflow.runName
// summary['Input File']          = params.input
// summary['Fasta File']             = params.fasta
// summary['GTF File']               = params.gtf
// summary['Max Resources']          = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
// if (workflow.containerEngine)     summary['Container'] = "$workflow.containerEngine - $workflow.container"
// summary['Output Dir']             = params.outdir
// summary['Launch Dir']             = workflow.launchDir
// summary['Working Dir']            = workflow.workDir
// summary['Script Dir']             = workflow.projectDir
// summary['User']                   = workflow.userName
// log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join('\n')
// log.info "-\033[2m--------------------------------------------------\033[0m-"


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists:true))
    .set {ch_gtf}

Channel
    .fromPath( params.input )
    .splitCsv(header:true)
    .map { row -> processRow(row) }
    .set { metadata }

metadata
    .filter{ it[0].sample_id == '10x_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/merged_counts/output") }] }
    .map { row -> listFiles(row, '.*.csv') }
    .flatMap { row -> row[1] }
    .set { ch_10x_counts }

metadata
    .filter{ it[0].sample_id == '10x_alignment_out' }
    .map { row -> row[1].collect{ file(it+"/velocyto/onefilepercell_ss8-TSS_P2_C10_75_and_others_TVIUH.loom", checkIfExists: true) } }
    .set { ch_10x_velocyto }

    // need to sort out the reading in from the metadata, reading it both counts and velocity from the same sample??



// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    //  test - print out the sample name

}
