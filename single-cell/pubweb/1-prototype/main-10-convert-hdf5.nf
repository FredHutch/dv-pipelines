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
dataset_name = wfi.parameters.input.name
dataset_type = wfi.parameters.input.type //new
pubweb_path = target_path + '/pubweb/'
hdf5_path = target_path + '/hdf5/'
src_path = target_path + '/src/'
source_path = "$params.s3source" + '/*.tar.gz'
scratch_path = '/opt/work'
s3_pubweb_source = 's3://dv-code-dev/pubweb/'


input_files = Channel.fromPath( source_path )

process CONVERT_HDF5 {
  echo true
  publishDir hdf5_path, mode: 'copy'

  input:
    path x from input_files
    val species
    val dataset_name
    val dataset_type

  output:
    file "output/*.tar.gz" into hdf5_ch

  script:
    """
    
    """
}


process PROCESS_PUBWEB {
  echo true
  publishDir pubweb_path, mode: 'copy'

  input:
    path x from hdf5_ch
    val species
    val dataset_name
    val dataset_type
    val s3_pubweb_source

  output:
    file "output/*" into pub_ch

  script:
    """
    mkdir -p input
    mkdir -p output
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
      --output 'output' \
      --name $dataset_name \
      --type $dataset_type \
      --species $species
    """
}