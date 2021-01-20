#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis as test} from "$baseDir/../modules/r_analysis/main.nf"


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/


Channel
    .fromPath( params.input ) //passed in as a param in the command line - path to local_samplesheet.csv
    .splitCsv(header:true)
    .map { row -> processRow(row) } //passing each row of the csv at a time into processRow, which makes dictionary linking sample id and path
    .set { metadata } //name of channel is metadata

metadata
    .filter{ it[0].sample_id == 'scRNA_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it) }] }
    .set { ch_scRNA }

// metadata
//     .filter{ it[0].sample_id == '10x_alignment_out' }
//     .map { row -> row[1].collect{ file(it+"/velocyto/onefilepercell_ss8-TSS_P2_C10_75_and_others_TVIUH.loom", checkIfExists: true) } }
//     .set { ch_10x_velocyto }

    // need to sort out the reading in from the metadata, reading it both counts and velocity from the same sample??



// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    //  Run differential expression analysis for lmx1a vs sox3U3
    test( params.modules['test'], ch_scRNA )
}



def processRow(LinkedHashMap row) { //defining a function in Java, this is read in before anything else
    def meta = [:] //meta is an empty dictionary
    meta.sample_id = row.sample_id //meta sample id will equal the sample id of the row we are looking at

    for (Map.Entry<String, ArrayList<String>> entry : row.entrySet()) {
        String key = entry.getKey();
        String value = entry.getValue();
        
        if(key != "sample_id" && key != "data") {
            meta.put(key, value)
        }
        
    }

    def array = [ meta, [ row.data ] ] //it can put in a whole list of paths which is useful if you want to read lots of files into one channel
    return array
}