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
hdf5_path = target_path
src_path = target_path + '/src/*.tar.gz'

scratch_path = '/opt/work'
s3_pubweb_source = 's3://dv-code-dev/pubweb/'


input_files = Channel.fromPath( src_path )
input_json_loc = Channel.fromPath ("${params.wfconfig}")

process CONVERT_MATRIXMARKET_TO_HDF5 {
  echo true
  publishDir hdf5_path, mode: 'copy'

  input:
    path x from input_files
    val species
    val dataset_name
    val dataset_type
    val s3_pubweb_source
    path input_json from input_json_loc

  output:
    file "output.hdf5" into hdf5_ch

  script:
    """
    mkdir -p input
    echo "input_json is $input_json"
    INPUTAR="\$(ls | grep .tar.gz)"
    echo "Local files are:"
    ls *
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
      --output 'output.hdf5'

    echo "List of output files"
    ls *
    """
}


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