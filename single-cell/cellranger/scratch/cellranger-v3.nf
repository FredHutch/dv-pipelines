#!/usr/bin/env nextflow
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

String configJSON = new File("${params.wfconfig}").text
def wfi = jsonSlurper.parseText(configJSON)
target_path = "$params.s3bucket" + "/"
target_path_val = Channel.value("$params.s3bucket" + "/")
source_path = wfi.parents[0].s3path + "/"
source_path_val = source_path

//Input parameters
/// Reference data
reference_path = wfi.parameters.input.reference_path
gex_reference = Channel.fromPath(reference_path + wfi.parameters.input.gex_reference)
vdj_reference = Channel.fromPath(reference_path + wfi.parameters.input.vdj_reference)

/// Workflow data
metadata = Channel.fromPath(target_path + wfi.parameters.input.metadata)
gex_fastq_paths = Channel.of(source_path + wfi.parameters.input.gex_fastq_path)
gex_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.gex_fastq_path)
vdj_fastq_paths = Channel.of(source_path + wfi.parameters.input.vdj_fastq_path)
vdj_fastq_files = Channel.fromPath(source_path + wfi.parameters.input.vdj_fastq_path)
study_id = Channel.from(wfi.parameters.input.study_id)
modes = Channel.from(wfi.parameters.aggr.modes)
run_gex = Channel.from(wfi.parameters.input.gex)
run_vdj = Channel.from(wfi.parameters.input.vdj)
fastq_type = Channel.value(wfi.parameters.count.fastq_type)


//Split input channels as needed
run_gex.into{map_gex ; count_gex ; aggr_gex ; mat_gex}
metadata.into{gex_metadata ; vdj_metadata}
run_vdj.into{map_vdj ; ana_vdj}
study_id.into{gex_study_id; vdj_study_id}



process TENX_GEX_MAP {
  echo true
  publishDir "$target_path_dir/Metadata", mode : 'copy'
  
  input:
    path meta_file from gex_metadata
    val fastq_path from gex_fastq_paths
    val study from gex_study_id
    val map from map_gex
    val fastq_type from fastq_type
    val target_path_dir from target_path_val

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
    cat("Meta file: $meta_file ")
    cat("Fastq type $fastq_type ")
    cat("Fastq path: $fastq_path ")
    cat("Target path: $target_path_val ")
    data <- read_csv("${meta_file}")
    cat('CSV rowcount: ', nrow(data))
    data <- data %>% filter(platform == '10XGenomics' & str_detect(locus, "GEX"))
    data <- data %>% mutate(library= str_extract(locus, "GEX"))
    data <- data %>% mutate(fastqType = "$fastq_type")
    # Account for folder structure differences between mkfastq and demux
    if ("$fastq_type" == "mkfastq" | "$fastq_type" == "bcl2fastq") {
      data <- data %>% mutate(fastqSample = sampleName)
      data <- data %>% mutate(fastqPath = "${fastq_path}")
    } else if ("$fastq_type" == "demux") {
      data <- data %>% mutate(fastqSample = paste(sampleName, "demux_id", sep="/"))
      data <- data %>% mutate(fastqFolder = "${fastq_path}")
      data <- data %>% mutate(fastqPath =  paste(fastqFolder, fastqSample, sep="/"))
    }
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    # Split VDJ and GEX libraries if present
    data_gex <- data %>% filter(library == "GEX")
    data_gex <- data_gex %>% mutate(molecule_h5 = paste("$target_path_val", "Counts", 
                                    repoName, "outs/molecule_info.h5", sep="/"))
    # Create samplesheet for mkfastq, and H5 and analysis sheets for VDJ analysis
    data_gex_h5 <- data_gex %>% rename(library_id = repoName)
    data_gex_h5 <- data_gex_h5 %>% select(library_id, molecule_h5, everything())
    data_gex_ss <- data_gex %>% select(repoName, indices) %>% rename(Sample = repoName, Index = indices) 
    data_gex_ss <- data_gex_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    # Write metadata files
    write_csv(data_gex_h5, "GEX_h5_samplesheet.csv")
    write_csv(data_gex_ss, "GEX_samplesheet.csv")
    cat('Done writing files. Done with R')
    """
}

process TENX_VDJ_MAP {
  echo true
  publishDir "$target_path_dir/Metadata", mode : 'copy'
  
  input:
    path meta_file from vdj_metadata
    val fastq_path from vdj_fastq_paths
    val study from vdj_study_id
    val map from map_vdj
    val fastq_type from fastq_type
    val target_path_dir from target_path_val

  output:
    path "VDJ_samplesheet.csv" into vdj_samplesheet
    path "VDJ_analysis_samplesheet.csv" into vdj_analysissheet
  
  when:
    map == 1
        
  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    # Metadata file provide should have the 12 mandatory fields mentioned in the README 
    data <- read_csv('${meta_file}')
    data <- data %>% filter(platform == '10XGenomics' & str_detect(locus, "VDJ"))
    data <- data %>% mutate(library = str_extract(locus, "VDJ"))
    data <- data %>% mutate(fastqType = "$fastq_type")
    # Account for folder structure differences between mkfastq and demux
    if ("$fastq_type" == "mkfastq" | "$fastq_type" == "bcl2fastq") {
      data <- data %>% mutate(fastqSample = sampleName)
      data <- data %>% mutate(fastqPath = "${fastq_path}")
    } else if ("$fastq_type" == "demux") {
      data <- data %>% mutate(fastqSample = paste(sampleName, "demux_id", sep="/"))
      data <- data %>% mutate(fastqFolder = "${fastq_path}")
      data <- data %>% mutate(fastqPath =  paste(fastqFolder, fastqSample, sep="/"))
    }
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    # Split VDJ and GEX libraries if present
    data_vdj <- data %>% filter(library == "VDJ")
    data_vdj <- data_vdj %>% mutate(vdj_sequences = paste("$target_path_val", "VDJ", 
                                    repoName, "outs/all_contig_annotations.csv", sep="/"))
    # Create samplesheet for mkfastq, and H5 and analysis sheets for VDJ analysis
    data_vdj_ss <- data_vdj %>% select(repoName, indices) %>% rename(Sample = repoName, Index = indices) 
    data_vdj_ss <- data_vdj_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    # Write metadata files
    write_csv(data_vdj_ss, "VDJ_samplesheet.csv")
    write_csv(data_vdj, "VDJ_analysis_samplesheet.csv")
    """
}


//Copy GEX output channels into channels for counting and aggreagtion
gex_h5sheet.into {count_gex_h5sheet; aggr_gex_h5sheet}

process TENX_COUNT {
  echo true 
  publishDir "$target_path_dir/Counts" , mode : 'copy'

  input:
    each sample from count_gex_h5sheet.splitCsv(header: true, quote: '\"')
    file sample_file from gex_fastq_files.collect()
    file gex_ref from gex_reference
    val run_count from count_gex
    val target_path_dir from target_path_val
  
  output:
    path "${sample.library_id}" into count_path
    val task.exitStaus into count_status

  when:
    run_count == 1

  script:
    if("$fastq_type" == "mkfastq" | "$fastq_type" == "bcl2fastq") {
      samplecmd = "--sample=$sample.fastqSample"
    } else { samplecmd = "" }

    """
    mkdir fastq
    mv $sample_file fastq
    du -sh fastq

    mkdir ref
    tar -zxf $gex_ref -C ref
    REFPATH="\$(ls ref/)"
    COMMAND="cellranger count --id=$sample.library_id --transcriptome=ref/\$REFPATH"
    COMMAND="\$COMMAND --fastqs=fastq --expect-cells=$sample.expected_cells"
    COMMAND="\$COMMAND --chemistry=$sample.chemistry"

    echo "Command: \$COMMAND"
    eval \$COMMAND
    """
}


process TENX_AGGR {
  echo true
  publishDir "$target_path_dir/Counts" , mode : 'copy'

  input:
    path samplesheet from aggr_gex_h5sheet
    each mode from modes
    val status from count_status.collect()
    val run_aggr from aggr_gex
    val target_path_dir from target_path_val
  
  output:
    path "Aggregate_${mode}_normalized" into aggr_path

  when:
    run_aggr == 1

  script:
    """
    mkdir counts
    mv $ counts

    cd counts
    cellranger aggr --id=Aggregate_${mode}_normalized --csv=${samplesheet} --normalize=${mode} 
    """
}

process TENX_MATRIX {
  echo true
  publishDir "$params.output.folder/Counts/${aggr_out}/out/filtered_feature_bc_matrix" , mode : 'copy'

  input:
    path aggr_out from aggr_path
    val run_matrix from mat_gex

  output:
    path "Filtered_expression_matrix.csv" into mtx_path

  when:
    run_matrix == 1

  script:
    """
    cellranger mat2csv ${aggr_out}/outs/filtered_feature_bc_matrix Filtered_expression_matrix.csv
    """
}

process TENX_VDJ {
  echo true
  publishDir "$params.output.folder/VDJ" , mode : 'copy'

  input:
    each sample from vdj_analysissheet.splitCsv(header: true, quote: '\"') 
    val analyze from ana_vdj

  output:
    path "${sample.sampleName}"

  when:
    analyze == 1

  script:
    if("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq") {
      samplecmd = "--sample=$sample.fastqSample"
    } else { samplecmd = "" }

    """
    mkdir fastq
    mv $samples fastq

    COMMAND="cellranger vdj --id=$sample.sampleName --reference=$params.input.vdj_reference"
    COMMAND="\$COMMAND --fastqs=$sample.fastqPath $samplecmd" 

    echo "Command: \$COMMAND"
    eval \$COMMAND
    """
}
