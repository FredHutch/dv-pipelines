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
target_path_hdf5 = target_path + 'pubweb/'
source_path = "$params.s3source"
scratch_path = '/opt/work'

sample_list = Channel.fromList(wfi.parameters.input.samples)
sample_list
  .map { [ it, file(fastq_path + '/' + it) ] }
  .into {view_folder_ch; read_folder_ch }


process CELLRANGER_COUNT {
  echo true
  publishDir target_path_count, mode: 'copy'

  input: 
    tuple val(x), file('sample/*') from read_folder_ch
    path genome, stageAs: 'genome/*' from reference_genome_path
    val target_path

  output:
    tuple val(x), path("count.tar.gz") into count_ch

  script:
    """
    echo "Sample folder contents"
    ls sample
    IGNORE="$x"
    ID="aligned"
    echo "ID is \$ID"s
    COMMAND="cellranger count --id=\$ID --transcriptome=$genome"
    COMMAND="\$COMMAND --fastqs=sample" 

    echo "Command: \$COMMAND"
    eval \$COMMAND

    tar -czvf count.tar.gz aligned
    """
}


process CELLRANGER_HDF5 {
  echo true
  publishDir target_path_hdf5, mode: 'copy'

  input:
    tuple val(x), path('count.tar.gz') from count_ch
    val species
    val dataset_name

  output:
    file "${x}/*" into pub_ch

  script:
    """
    mkdir -p input
    tar -xzf count.tar.gz -C input
    echo "Contents of local, \$(pwd)"
    echo "Finding the file location"
    find . -name "filtered_feature_bc_matrix.h5"
    echo "Finding the PCA location"
    find . -type d -name "pca"
    
    mkdir -p $x

    # reinstall the library
    export LIBRARYDIR=/opt/pubweb
    OLDDIR=\$PWD
    rm -rf /opt/pubweb
    mkdir -p \$LIBRARYDIR
    aws s3 cp s3://dv-code-dev/pubweb/ \$LIBRARYDIR --recursive
    python -m pip install /opt/pubweb
    
    python /opt/pubweb/pubweb/invoke-cellranger.py \
      --input 'input/aligned/outs' \
      --output '$x' \
      --name $dataset_name \
      --species $species
    """
}
