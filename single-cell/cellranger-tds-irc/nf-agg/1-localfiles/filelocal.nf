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


sample_folders
  .map { fastq_path + '/' + it }
  .into {view_folder_ch; read_folder_ch}
//folder_ch.into {view_folder_ch; read_folder_ch}

view_folder_ch.view()



process SINGLE_FILE {
  echo true
  input:
    path textx from read_folder_ch

  output:
    file *.txt into count_ch

  script:
    """
    echo "Looking at \$PWD $textx"
    ls $textx
    echo "fakeentry" > $textx\\testoutput.txt
    """
}


process AGG_FILE {
  echo true
  input:
    path textx from count_ch

  script:
    """
    echo "Looking at \$PWD $textx"
    ls $textx
    echo "fakeentry" > $textx\\testoutput.txt
    """
}

