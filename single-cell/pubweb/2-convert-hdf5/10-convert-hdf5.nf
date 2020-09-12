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
pubweb_path = target_path + '/pubweb/'
hdf5_path = target_path + '/hdf5/'
src_path = target_path + '/src/'

source_path = "$params.s3source" + '/*.tar.gz'
scratch_path = '/opt/work'
s3_pubweb_source = 's3://dv-code-dev/pubweb/'


input_files = Channel.fromPath( source_path )
input_json_loc = Channel.fromPath ("${params.wfconfig}")

process CONVERT_MATRIXMARKET_TO_HDF5 {
  echo true
  publishDir hdf5_path, mode: 'copy'

  input:
    path x from input_files
    val species
    val dataset_name
    val dataset_type
    path input_json from input_json_loc

  output:
    file "output/*" into pub_ch

  script:
    """
    mkdir -p output
    mkdir -p input
    echo "input_json is $input_json . Contents:"
    cat $input_json
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

    python \$LIBRARYDIR/pubweb/convert-to-hdf5.py --params $input_json \
      --matrix 'input/matrix' --var 'input/var' --obs 'input/obs' \
      --output 'output/output.hdf5'

    echo "List of output files"
    ls output/*
    """
}

