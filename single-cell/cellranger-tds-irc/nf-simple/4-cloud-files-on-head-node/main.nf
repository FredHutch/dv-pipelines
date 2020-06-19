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

sample_list = Channel.from(wfi.parameters.input.samples)



//cellranger count
// --id=GS_PBMC_8_2013_MCCB1050
//    --transcriptome=/scratch/data/genome/hg38_merkerCellPolyomavirus_CellRanger_v2
//    --sample=GS_PBMC_8_2013_MCCB1050
//    --fastqs=/scratch/data/patient1/outs/fastq_path/HKMNCBCXY/GS_PBMC_8_2013_MCCB1050

process CELLRANGER_COUNT {
  echo true
  input: 
    val sample from sample_list
    val target_path
    val source_path
    val reference_genome_path

  script:
    """
    echo "Reference genome $reference_genome_path"
    echo "Source is $source_path/$sample"
    echo "Target is $target_path/$sample"
    COMMAND="cellranger count --id=$sample --transcriptome=$reference_genome_path"
    COMMAND="\$COMMAND --fastqs=$fastq_path/$sample" 

    echo "Command: \$COMMAND"
    """
}
