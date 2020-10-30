#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

//Input parameters
/// Reference data
reference_genome_path = wfi.parameters.input.genome_path

species = wfi.parameters.input.species
dataset_name = wfi.parameters.input.name //change this?
target_path = "$params.s3target"
target_path_count = target_path + '/count/'
target_path_agg = target_path + '/agg/'
target_path_hdf5 = target_path + '/pubweb/'
source_path = "$params.s3source"
scratch_path = '/opt/work'
agg_path = 'input/count'


process GET_SAMPLE_LIST {
  echo true

  input: 
    val source_path

  output:
    file "sample*" into sample_list, sample_list_2

  script:
    """
    apt-get -y install perl
    # get the list of prefixes
    # parse to just the prefix name, strip out trailing /
    # split into one prefix per file
    aws s3 ls $source_path | grep 'PRE ' | awk '{print \$2}' \
     | rev | cut -c 2- | rev | split -l 1 - sample
    # now remove the trailing newLine
    perl -p -i -e 's/\\R//g;' sample*

    # print the list of samples for debugging purposes
    echo "Sample list:"
    cat sample*
    """
}


// generate the mapping for cellranger counts
sample_list
  .flatMap()
  .map { [ it.text, file(source_path + it.text) ] }
  .set { read_folder_ch }


process CELLRANGER_COUNT {
  echo true
  publishDir target_path_count, mode: 'copy'
  errorStrategy 'finish'

  input: 
    tuple val(x), file('sample/*') from read_folder_ch
    path genome, stageAs: 'genome/*' from reference_genome_path
    val target_path

  output:
    path("*.tar.gz") into count_ch
    path("${x}/*") into debugging_ch

  // 'sample' has the input fastq files
  // $x (a.k.a $ID) is where the outputs go
  script:
    """
    ID="$x"
    mkdir -p $x
    echo "ID is \$ID"
    echo "Sample folder contents"
    ls sample/*

    COMMAND="cellranger count --id=\$ID --transcriptome=$genome"
    COMMAND="\$COMMAND --fastqs=sample"

    echo "Command: \$COMMAND"
    eval \$COMMAND || true

    tar -czf "\$ID.tar.gz" \$ID/*
    """
}



// generate the agg csv for cellranger aggr
sample_list_2
  .flatMap()
  .map {
  "${it.text},PLACEHOLDERDIR/${it.text}/outs/molecule_info.h5"
}.collectFile(
    name: 'molecule_info.csv',
    newLine: true,
    seed: "library_id,molecule_h5"
).set { agg_info_csv }




process CELLRANGER_AGG {
  echo true
  publishDir target_path_agg, mode: 'copy'

  input: 
    path tarball from count_ch.collect()
    path "sourcemap.csv" from agg_info_csv

  output:
    path("agg.tar.gz") into agg_ch

  script:
    """
    echo "List of tarballs"
    ls *
    echo "List of tarball contents"
    # ls *.tar.gz | xargs -I % sh -c 'tar -tzf %'
    ls *.tar.gz | xargs -I % sh -c 'tar -xzf %'
    rm *.tar.gz

    echo "List of un-tar'd folders"
    ls -d *
    
    echo "CSV mapping:"
    sed 's@PLACEHOLDERDIR@'"\$PWD"'@' sourcemap.csv > mapping.csv
    cat mapping.csv

    ID="aggregated"
    COMMAND="cellranger aggr --id=\$ID --csv=mapping.csv"

    echo "Command: \$COMMAND"
    eval \$COMMAND
    
    echo "Cellranger agg output" > aggregated/semaphore
    cd \$ID
    tar -czf agg.tar.gz *
    mv agg.tar.gz ../
    """
}




process CELLRANGER_HDF5 {
  echo true
  publishDir target_path_hdf5, mode: 'copy'

  input:
    path('agg.tar.gz') from agg_ch
    val species
    val dataset_name

  output:
    file "*" into pub_ch

  script:
    """
    mkdir -p input
    tar -xzf agg.tar.gz -C input
    echo "Contents of local, \$(pwd)"
    echo "Finding the file location"
    find -L . -name "filtered_feature_bc_matrix.h5"
    echo "Finding the PCA location"
    find -L . -type d -name "pca"
    
    mkdir -p output

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    OLDDIR=\$PWD
    rm -rf /opt/pubweb
    mkdir -p \$LIBRARYDIR
    aws s3 cp s3://dvc-wf-metadata/code/pubweb/ \$LIBRARYDIR --recursive
    python -m pip install /opt/pubweb
    
    python /opt/pubweb/pubweb/invoke-cellranger.py \
      --input 'input/outs' \
      --output 'output' \
      --name $dataset_name \
      --species $species

    echo "List of folders"
    ls -d .
    echo "Listing output"
    ls output/*
    """
}