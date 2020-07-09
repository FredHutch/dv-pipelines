#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)
source_path = "$params.s3bucket" + "/"

//Input parameters
/// Reference data
reference_path = wfi.parameters.reference_path
gex_reference = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)

/// Workflow data
metadata = Channel.fromPath(source_path + wfi.parameters.input.metadata)
gex_fastq_paths = Channel.of(source_path + wfi.parameters.input.gex_fastq_path)
gex_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
vdj_fastq_paths = Channel.of(source_path + wfi.parameters.input.vdj_fastq_path)
vdj_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.vdj_fastq_path)
study_id = Channel.from(wfi.parameters.input.study_id)
modes = Channel.from(wfi.parameters.aggr.modes)
run_gex = Channel.from(wfi.parameters.input.gex)
run_vdj = Channel.from(wfi.parameters.input.vdj)
fastq_type = Channel.from(wfi.parameters.count.fastq_type)


//Split input channels as needed
run_gex.into{map_gex ; count_gex ; aggr_gex ; mat_gex}
metadata.into{gex_metadata ; vdj_metadata}
run_vdj.into{map_vdj ; ana_vdj}
study_id.into{gex_study_id; vdj_study_id}



process TENX_GEX_MAP {
  echo true
  publishDir "$params.output.folder/Metadata", mode : 'copy'
  
  input:
    file meta_file from gex_metadata
    val fastq_path from gex_fastq_paths
    val study from gex_study_id
    val map from map_gex

  output:
    path "GEX_h5_samplesheet.csv" into gex_h5sheet
    path "GEX_samplesheet.csv" into gex_samplesheet
  
  when:
    map == 1
        
  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    # Metadata file provide should have the 12 mandatory fields mentioned in the README 
    data <- read_csv("${meta_file}")
    data <- data %>% filter(platform == '10XGenomics' & str_detect(locus, "GEX"))
    data <- data %>% mutate(library= str_extract(locus, "GEX"))
    data <- data %>% mutate(fastqType = "$params.count.fastq_type")
    # Account for folder structure differences between mkfastq and demux
    if ("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq") {
      data <- data %>% mutate(fastqSample = sampleName)
      data <- data %>% mutate(fastqPath = "${fastq_path}")
    } else if ("$params.count.fastq_type" == "demux") {
      data <- data %>% mutate(fastqSample = paste(sampleName, "demux_id", sep="/"))
      data <- data %>% mutate(fastqFolder = "${fastq_path}")
      data <- data %>% mutate(fastqPath =  paste(fastqFolder, fastqSample, sep="/"))
    }
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    # Split VDJ and GEX libraries if present
    data_gex <- data %>% filter(library == "GEX")
    data_gex <- data_gex %>% mutate(molecule_h5 = paste("$params.output.folder", "Counts", 
                                    repoName, "outs/molecule_info.h5", sep="/"))
    # Create samplesheet for mkfastq, and H5 and analysis sheets for VDJ analysis
    data_gex_h5 <- data_gex %>% rename(library_id = repoName)
    data_gex_h5 <- data_gex_h5 %>% select(library_id, molecule_h5, everything())
    data_gex_ss <- data_gex %>% select(repoName, indices) %>% rename(Sample = repoName, Index = indices) 
    data_gex_ss <- data_gex_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    # Write metadata files
    write_csv(data_gex_h5, "GEX_h5_samplesheet.csv")
    write_csv(data_gex_ss, "GEX_samplesheet.csv")
    """
}

//Copy GEX output channels into channels for counting and aggreagtion
gex_h5sheet.into {count_gex_h5sheet; aggr_gex_h5sheet}

process TENX_COUNT {
  echo true 
  publishDir "$params.output.folder/Counts" , mode : 'copy'

  input:
    each sample from count_gex_h5sheet.splitCsv(header: true, quote: '\"')
    file samples from gex_fastq_files.collect()
    val run_count from count_gex
  
  output:
    path "${sample.library_id}" into count_path
    val task.exitStaus into count_status

  when:
    run_count == 1

  script:
    if("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq") {
      samplecmd = "--sample=$sample.fastqSample"
    } else { samplecmd = "" }

    """
    mkdir fastq
    mv $samples fastq

    COMMAND="cellranger count --id=$sample.library_id --transcriptome=$params.input.gex_reference"
    COMMAND="\$COMMAND --fastqs=fastq --expect-cells=$sample.expected_cells"
    COMMAND="\$COMMAND --chemistry=$sample.chemistry"

    echo "Command: \$COMMAND"
    # eval \$COMMAND
    """
}

