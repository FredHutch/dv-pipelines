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

sample_list = Channel.fromList(wfi.parameters.input.samples)
sample_list.into { read_list, agg_list }

read_list
  .map { [ it, file(fastq_path + '/' + it) ] }
  .into { read_folder_ch }


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



// use the sample sheet and the output of the count process to make
// a new sample sheet for the aggregate process
molecule_info.join(sampleSheetRows).map {
    it[2].remove('library_id')
    values = it[2].values().join(',')
    return [it[0], it[1], values].join(',')
}.collectFile(
    name: 'molecule_info.csv',
    newLine: true,
    seed: "library_id,molecule_h5," + keys.drop(1).join(',')
).set { agg_info_csv }



process MAKE_AGG_CSV {
  echo true
  input: 
    val(x) from agg_list

  output:
    file agg.csv into aggcsv_ch
    tuple val(x), path("count.tar.gz") into 

  // make the CSV file
  script:
    """
    
    """
}


process CELLRANGER_AGG {
  echo true
  publishDir target_path_agg, mode: 'copy'

  input: 
    path sample, stageAs: 'sample/*' from count_ch
    path "molecule_info.csv" from agg_info_csv

  output:
    tuple val(x), path("agg.tar.gz") into agg_ch

  script:
    """
    ID="aggregated"
    COMMAND="cellranger aggr --id=\$ID --csv=molecule_info.csv" 

    echo "Command: \$COMMAND"
    eval $COMMAND

    tar -czvf agg.tar.gz aggregated
    """
}



process CELLRANGER_HDF5 {
  echo true
  publishDir target_path_hdf5, mode: 'copy'

  input:
    tuple val(x), path('agg.tar.gz') from agg_ch
    val species
    val dataset_name

  output:
    file "${x}/*" into pub_ch

  script:
    """
    mkdir -p input
    tar -xzf agg.tar.gz -C input
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
