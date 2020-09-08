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

process GET_DATA {
  echo true
  publishDir src_path, mode: 'copy'

  input:
    val species
    val dataset_name
    val dataset_type
    path input_json from input_json_loc

  output:
    path "input.tar.gz" into pub_ch

  script:
    """
    mkdir -p output
    echo "input_json is $input_json . Contents:"
    cat $input_json \
      | jq -r '.parameters.input.source | map( [.name, .file] | join(", ")) | join("\n")' > input.csv
    # now get the data
    while IFS="," read -r name url
    do
      echo "Now getting \$name from \$url"
      curl \$url -o \$name
    done < input.csv
    rm input.csv
    # tar the results
    find . -type f -print | tar -czf output/input.tar.gz --files-from -
    # move the tarball to the publish location
    mv output/input.tar.gz .
    """
}

