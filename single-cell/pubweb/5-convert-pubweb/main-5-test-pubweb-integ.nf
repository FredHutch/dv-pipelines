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
pubweb_path = target_path
scratch_path = '/opt/work'
s3_pubweb_source = 's3://dv-code-dev/pubweb/'


pub_ch = Channel.fromPath( "$params.s3source" )



process PROCESS_PUBWEB {
  echo true
  publishDir pubweb_path, mode: 'copy'

  input:
    path x from pub_ch
    val species
    val dataset_name
    val dataset_type
    val s3_pubweb_source

  output:
    path "pubweb/*" into hdf5_ch

  script:
    """
    echo "List of input files"
    ls *
    mkdir -p pubweb

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    OLDDIR=\$PWD
    rm -rf \$LIBRARYDIR
    mkdir -p \$LIBRARYDIR
    aws s3 cp $s3_pubweb_source \$LIBRARYDIR --recursive
    python -m pip install \$LIBRARYDIR

    python \$LIBRARYDIR/pubweb/invoke-pubweb.py \
      --input '$x' \
      --output 'pubweb' \
      --name $dataset_name \
      --type $dataset_type \
      --species $species
    """
}