#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)

//Input parameters
/// Reference data
reference_genome_path = wfi.parameters.input.genome_path
fastq_path = wfi.parameters.input.fastq_path

target_path = "$params.source"
source_path = "$params.target"
scratch_path = '/opt/work'

sample_list = Channel.fromList(wfi.parameters.input.samples)
sample_list
  .map { fastq_path + '/' + it }
  .into {view_folder_ch; read_folder_ch}


process CELLRANGER_COUNT {
  echo true
  input: 
    path sample, stageAs: 'sample/*' from read_folder_ch
    path genome, stageAs: 'genome/*' from reference_genome_path
    val target_path
    val source_path

  output:
    file * into count_ch, count_ch_2

  script:
    """
    echo "Reference genome contents:"
    ls $genome
    echo "Sample folder contents"
    ls $sample
    ID="placeholder"
    echo "ID is \$ID"
    COMMAND="cellranger count --id=\$ID --transcriptome=$genome"
    COMMAND="\$COMMAND --fastqs=$sample" 

    echo "Command: \$COMMAND"
    eval \$COMMAND
    """
}

process MAKE_AGG_CSV {
  echo true
  input: 
    path sample, stageAs: 'sample/*' from count_ch_2

  output:
    file agg.csv into aggcsv_ch

  // make the CSV file
  script:
    """
    
    """
}


process CELLRANGER_AGG {
  echo true
  input: 
    path sample, stageAs: 'sample/*' from count_ch
    path aggcsv from aggcsv_ch
    val target_path
    val source_path

  script:
    """
    $\ID="placeholder"
    COMMAND="cellranger aggr --id=\$ID --csv=$aggcsv" 

    echo "Command: \$COMMAND"
    """
}
