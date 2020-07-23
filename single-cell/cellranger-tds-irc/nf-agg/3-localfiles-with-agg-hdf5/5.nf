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
dataset_name = wfi.parameters.input.name //change this?
target_path = "$params.s3target"
target_path_count = target_path + 'count/'
target_path_agg = target_path + 'agg/'
target_path_hdf5 = target_path + 'pubweb/'
source_path = "$params.s3source"
scratch_path = '/opt/work'
agg_path = 'input/count'


Channel.fromList(wfi.parameters.input.samples)
   .into { sample_list; sample_list_2}

// generate the mapping for cellranger counts
sample_list
  .map { [ it, file(fastq_path + '/' + it) ] }
  .set { read_folder_ch }


process SINGLE_FILE {
  echo true
  publishDir target_path_count, mode: 'copy'

  input:
    tuple val(x), file('sample/*') from read_folder_ch

  output:
    path("*.tar.gz") into count_ch

  script:
    """
    ID="$x"
    echo "x is $x"
    mkdir -p $x
    mv sample/* $x

    cd $x

    tar -czvf "\$ID.tar.gz" *
    mv *.tar.gz ../
    """
}


// generate the agg csv for cellranger aggr
sample_list_2.map {
  "${it},PLACEHOLDERDIR/$it/molecule_info.h5"
}.collectFile(
    name: 'molecule_info.csv',
    newLine: true,
    seed: "library_id,molecule_h5"
).set { agg_info_csv }


process AGG_FILE {
  echo true
  publishDir target_path_agg, mode: 'copy'

  input: 
    path('*.tar.gz') from count_ch.collect()
    path "sourcemap.csv" from agg_info_csv

  output:
    path("agg/*") into agg_ch

  script:
    """
    ls *.tar.gz | xargs -I % sh -c 'tar -xzf %'
    rm *.tar.gz

    echo "Finding text files"
    find -L . -name *.txt
    
    echo "Current working directory is \$PWD"
    echo "Looking for the CSV"
    cat sourcemap.csv
    echo "Changing the source"
    sed 's@PLACEHOLDERDIR@'"\$PWD"'@' sourcemap.csv > mapping.csv
    cat mapping.csv

    tar -czvf "agg.tar.gz" *
    mkdir -p agg
    mv agg.tar.gz agg
    """
}


process HDF5_FILE {
  echo true
  publishDir target_path_hdf5, mode: 'copy'

  input:
    path('agg.tar.gz') from agg_ch
    val species
    val dataset_name

  output:
    file "output/*" into pub_ch

  script:
    """
    mkdir -p input
    tar -xzf agg.tar.gz -C input
    echo "Contents of local, \$(pwd)"
    echo "Finding the file location"
    find -L . -name *.txt
    
    mkdir -p output
    echo 'Output dir' > output/note.txt
    cp -r input output/
    """
}
