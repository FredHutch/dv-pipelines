#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

//Input parameters
/// Reference data
species = wfi.parameters.input.species
dataset_name = wfi.parameters.input.name
dataset_type = wfi.parameters.input.type //new
target_path = "$params.s3target"
source_path = "$params.s3source" + '/*.tar.gz'
scratch_path = '/opt/work'
s3_pubweb_source = 's3://dv-code-dev/pubweb/'


input_files = Channel.fromPath( source_path )


process PROCESS_PUBWEB {
  echo true
  publishDir target_path, mode: 'copy'

  input:
    path x from input_files
    val species
    val dataset_name
    val dataset_type
    val s3_pubweb_source

  output:
    path "pubweb/*" into pub_ch

  script:
    """
    mkdir -p input
    mkdir -p pubweb
    INPUTAR="\$(ls | grep .tar.gz)"
    echo "Now untar'ing \$INPUTAR"
    tar -xzf \$INPUTAR -C input
    echo "List of untar'd files"
    ls input/*

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    OLDDIR=\$PWD
    rm -rf \$LIBRARYDIR
    mkdir -p \$LIBRARYDIR
    aws s3 cp $s3_pubweb_source \$LIBRARYDIR --recursive
    python -m pip install \$LIBRARYDIR

    python \$LIBRARYDIR/pubweb/invoke-pubweb.py \
      --input 'input' \
      --output 'pubweb' \
      --name $dataset_name \
      --type $dataset_type \
      --species $species
    """
}