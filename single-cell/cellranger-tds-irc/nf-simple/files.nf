#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)



//Input parameters
/// Reference data
fastq_path = wfi.parameters.input.fastq_path

target_path = "$params.source"
source_path = "$params.target"

sample_folders = Channel.fromList(wfi.parameters.input.samples)

Channel
    .from(sample_folders)
    .map { fold -> ["sample$fold", file(fastq_path + '/' + fold)] }
    .set ( folder_ch )
//folder_ch.into {view_folder_ch; debug_folder_ch; read_folder_ch}
//read_file_ch = read_folder_ch

folder_ch.view()


process DEBUG_VAL {
  echo true
  input:
    val fastq_path
    val debug_folder_ch

  script:
    """
    echo "File channel is $debug_folder_ch"
    """
}


process SINGLE_FILE {
  echo true
  input:
    path textx from read_file_ch

  script:
    """
    echo "Test file val is $textx"
    """
}
