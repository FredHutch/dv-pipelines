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


process CREATE_OUTPUTS {
  echo true
  input:
    path test from read_folder_ch

  output:
    file 'placeholder/*' into testout_ch

  script:
    """
    SOURCEPATH="$test"
    echo "Path is \$SOURCEPATH"
    mkdir -p placeholder
    mkdir -p placeholder/\$SOURCEPATH
    echo "Now running cat"
    cat \$SOURCEPATH/* >> placeholder/\$SOURCEPATH/test.txt

    """
}

process READ_OUTPUTS {
  echo true
  input:
    path test from testout_ch.flatMap()

  script:
    """
    if [ -f "$test" ]
    then
      echo "File $test"
      cat $test
    else
      echo "Directory $test"
      cd $test
      cat *
    fi
    """
}