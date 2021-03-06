#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

//Input parameters
/// Reference data
reference_genome_path = wfi.parameters.input.genome_path
fastq_path = wfi.parameters.input.fastq_path

species = wfi.parameters.input.species
dataset_name = wfi.parameters.input.name //change this?
target_path = "$params.s3target"
target_path_count = target_path + 'count/'
target_path_agg = target_path + 'agg/'
target_path_hdf5 = target_path + 'pubweb/'
source_path = "$params.s3source"
scratch_path = '/opt/work'
agg_path = 'input/count'


Channel.fromList(wfi.parameters.input.samples)
   .into { sample_list; sample_list_2}
//sample_list.set { read_list; agg_list }

sample_list
  .map { [ it, file(fastq_path + '/' + it) ] }
  .into { read_folder_ch }




sample_list_2.map {
  "${it},$agg_path/$it/molecule_info.h5"
}.collectFile(
    name: 'molecule_info.csv',
    newLine: true,
    seed: "library_id,molecule_h5"
).into { agg_info_csv }

agg_info_csv.subscribe {  println "Got: $it"  }

