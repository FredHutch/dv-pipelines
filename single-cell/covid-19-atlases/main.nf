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
source_url = wfi.parameters.input.url
target_path = "$params.s3target"

s3_pubweb_source = 's3://dvc-wf-metadata/code/pubweb/'

input_json_loc = Channel.fromPath ("${params.wfconfig}")

process GET_DATA_AND_PUBWEB {
  echo true
  publishDir target_path, mode: 'copy'

  input:
    val species
    val dataset_name
    val dataset_type
    val source_url

  output:
    path "*" into pub_ch

  script:
    """
    apt-get -y update
    apt-get -y install curl
    # now get the data
    echo "source_url is $source_url"
    curl $source_url -o input.h5ad

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    rm -rf \$LIBRARYDIR
    mkdir -p \$LIBRARYDIR
    aws s3 cp $s3_pubweb_source \$LIBRARYDIR --recursive
    python -m pip install \$LIBRARYDIR
    
    # convert the source data to an hdf5 file
    echo "Converting source file to hdf5"
    python -m pip install anndata
    python \$LIBRARYDIR/pubweb/convert-covid-h5ad.py \
      --input input.h5ad \
      --output convert.hdf5

    # now run pubweb
    mkdir -p pubweb
    echo "Converting hdf5 to pubweb"

    python \$LIBRARYDIR/pubweb/invoke-pubweb.py \
      --input 'convert.hdf5' \
      --output 'pubweb' \
      --name $dataset_name \
      --type $dataset_type \
      --species $species

    rm input.h5ad convert.hdf5
    mv pubweb/* .
    rm -rf pubweb/*
    """
}

