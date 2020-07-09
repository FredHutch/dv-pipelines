#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

//Input parameters
/// Reference data
reference_genome_path = wfi.parameters.input.genome_path
fastq_path = wfi.parameters.input.fastq_path

target_path = "$params.target"
source_path = "$params.source"
scratch_path = '/opt/work'

sample_list = Channel.fromList(wfi.parameters.input.samples)

sample_list
  .map { [ it, file(fastq_path + '/' + it) ] }
  .into {view_folder_ch; read_folder_ch}


process SINGLE_FILE {
  echo true
  input:
    set val(x), file('sample/*') from read_folder_ch

  script:
    """
    echo "Looking at \$PWD"
    echo "It value is $x"
    echo "Folder contents:"
    ls sample
    echo "fakeentry" > sample/testoutput.txt
    """
}
