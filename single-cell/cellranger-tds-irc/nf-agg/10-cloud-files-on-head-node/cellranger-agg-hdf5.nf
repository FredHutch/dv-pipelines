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




process CELLRANGER_COUNT {
  echo true
  publishDir target_path_count, mode: 'copy'

  input: 
    tuple val(x), file('sample/*') from read_folder_ch
    path genome, stageAs: 'genome/*' from reference_genome_path
    val target_path

  output:
    path("*.tar.gz") into count_ch

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
    eval \$COMMAND

    cd $x

    tar -czf "\$ID.tar.gz" *
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




process CELLRANGER_AGG {
  echo true
  publishDir target_path_agg, mode: 'copy'

  input: 
    path('*.tar.gz') from count_ch.collect()
    path "sourcemap.csv" from agg_info_csv

  output:
    path("agg.tar.gz") into agg_ch

  script:
    """
    ls *.tar.gz | xargs -I % sh -c 'tar -xzf %'
    rm *.tar.gz
    
    echo "Changing the CSV mapping"
    sed 's@PLACEHOLDERDIR@'"\$PWD"'@' sourcemap.csv > mapping.csv

    ID="aggregated"
    mkdir -p \$ID
    COMMAND="cellranger aggr --id=\$ID --csv=mapping.csv" 

    echo "Command: \$COMMAND"
    eval $COMMAND

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
    file "output/*" into pub_ch

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
    aws s3 cp s3://dv-code-dev/pubweb/ \$LIBRARYDIR --recursive
    python -m pip install /opt/pubweb
    
    python /opt/pubweb/pubweb/invoke-cellranger.py \
      --input 'input/aligned/outs' \
      --output 'output' \
      --name $dataset_name \
      --species $species
    """
}