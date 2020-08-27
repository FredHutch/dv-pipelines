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
target_path = "$params.s3target" + '/pubweb/'
source_path = "$params.s3source"
scratch_path = '/opt/work'
agg_path = 'input/count'
s3_pubweb_source = 's3://dvc-wf-metadata/code/pubweb/'




process PROCESS_PUBWEB {
  echo true
  publishDir target_path, mode: 'copy'

  input:
    path('source.tar.gz') from agg_ch
    val species
    val dataset_name
    val dataset_type
    val s3_pubweb_source

  output:
    file "*" into pub_ch

  script:
    """
    mkdir -p input
    mkdir -p output
    tar -xzf source.tar.gz -C input

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    OLDDIR=\$PWD
    rm -rf /opt/pubweb
    mkdir -p \$LIBRARYDIR
    aws s3 cp $s3_pubweb_source \$LIBRARYDIR --recursive
    python -m pip install /opt/pubweb

    if [ $dataset_type == "cellranger" ]
    then
      python /opt/pubweb/pubweb/invoke-cellranger.py \
        --input 'input/outs' \
        --output 'output' \
        --name $dataset_name \
        --species $species
    fi
    elif [ $dataset_type == "anndata ]
    then
      python /opt/pubweb/pubweb/invoke-anndata.py \
        --input 'input' \
        --output 'output' \
        --name $dataset_name \
        --species $species
    fi

    echo "List of folders"
    ls -d .
    echo "Listing output"
    ls output/*
    """
}