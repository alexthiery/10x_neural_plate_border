#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis as test_1} from "$baseDir/../modules/tools/r_analysis/main.nf" addParams(options: modules['test_1'])
include {r_analysis as test_2} from "$baseDir/../modules/tools/r_analysis/main.nf" addParams(options: modules['test_2'])


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .fromPath( params.input ) //passed in as a param in the command line - path to local_samplesheet.csv
    .splitCsv(header:true)
    .map { row -> processRow(row) } //passing each row of the csv at a time into processRow, which makes dictionary linking sample id and path
    .set { metadata } //name of channel is metadata


metadata
    .filter{ it[0].sample_id == 'test_data_1' }
    .map { row -> row[1].collect{file(it)} }
    .set { ch_test_data_1 }

// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    test_1 (ch_test_data_1)
    test_2 (test_1.out)
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