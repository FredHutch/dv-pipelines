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
    set val(x), file('sample/*') from read_folder_ch
    path genome, stageAs: 'genome/*' from reference_genome_path
    val target_path

  output:
    file "${x}/*" into count_ch
    file 'semaphore.txt' into count_hdf5_ch

  script:
    """
    echo "Reference genome contents:"
    ls $genome
    echo "Sample folder contents"
    ls sample
    ID="$x"
    echo "ID is \$ID"
    COMMAND="cellranger count --id=\$ID --transcriptome=$genome"
    COMMAND="\$COMMAND --fastqs=sample" 

    echo "Command: \$COMMAND"
    eval \$COMMAND

    echo "$x" > semaphore.txt
    """
}


process CELLRANGER_HDF5 {
  echo true
  publishDir target_path_hdf5, mode: 'copy'

  input:
    path counts, stageAs: 'input/*' from target_path_count
    path semaphore, stageAs: 'wait/*' from count_hdf5_ch
    val species
    val dataset_name

  output:
    file "output/*" into pub_ch

  script:
    """
    mkdir -p output
    echo "Contents of local, \$(pwd)"
    ls -lhR $counts
    
    python /opt/pubweb/invoke-cellranger.py \
      --input '$counts' \
      --output 'output' \
      --name $dataset_name \
      --species $species
    """
}
