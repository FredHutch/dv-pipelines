#!/usr/bin/env nextflow

metadata = Channel.fromPath(params.input.metadata)
gex_reference = Channel.fromPath(params.input.gex_reference)
vdj_reference = Channel.fromPath(params.input.vdj_reference)
fastq_paths = Channel.from(params.input.fastq_paths)
study_id = Channel.from(params.input.study_id)
modes = Channel.from(params.aggr.modes)
run_gex = Channel.from(params.input.gex)
run_vdj = Channel.from(params.input.vdj)
fastq_type = Channel.from(params.count.fastq_type)


run_gex.into{ count_gex ; aggr_gex ; mat_gex}

process TENX_PROCESS {
  echo false
  publishDir "$params.output.folder/Metadata"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    file 'metadata.csv' from metadata
    val fastq_path from fastq_paths
    val study from study_id

  output:
    path "${study}_GEX_h5_samplesheet.csv" into gex_h5sheet
    path "${study}_GEX_samplesheet.csv" into gex_samplesheet
    path "${study}_VDJ_samplesheet.csv" into vdj_samplesheet
    path "${study}_VDJ_analysis_samplesheet.csv" into vdj_analysissheet
        
  script:

  if ("$params.count.fastq_type" == 'mkfastq')
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    data <- read_csv('metadata.csv')
    data <- data %>% filter(platform == '10XGenomics')
    data <- data %>% mutate(library= str_extract(locus, "(GEX|VDJ)"))
    data <- data %>% mutate(fastqSample = paste(library, sampleName, sep="_"))
    data <- data %>% mutate(fastqFolder = "${fastq_path}")
    data <- data %>% mutate(fastqPath =  paste(fastqFolder, fastqSample, sep="/${study}/"))
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    data <- data %>% mutate(molecule_h5 = paste("$params.output.folder", "Counts", sampleName, "outs/molecule_info.h5", sep="/"))
    data_gex <- data %>% filter(library == 'GEX')
    data_vdj <- data %>% filter(library == 'VDJ')
    data_gex_h5 <- data_gex %>% rename(library_id = sampleName)
    data_gex_h5 <- data_gex_h5 %>% select(library_id, molecule_h5, everything())
    data_gex_ss <- data_gex %>% select(sampleName, indices) %>% rename(Sample = sampleName, Index = indices) 
    data_gex_ss <- data_gex_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    data_vdj_ss <- data_vdj %>% select(sampleName, indices) %>% rename(Sample = sampleName, Index = indices) 
    data_vdj_ss <- data_vdj_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    write_csv(data_gex_h5, paste("${study}", "GEX_h5_samplesheet.csv", sep="_"))
    write_csv(data_gex_ss, paste("${study}", "GEX_samplesheet.csv", sep="_"))
    write_csv(data_vdj_ss, paste("${study}", "VDJ_samplesheet.csv", sep="_"))
    write_csv(data_vdj, paste("${study}", "VDJ_analysis_samplesheet.csv", sep="_"))
    """
  else if ("$params.count.fastq_type" == 'demux')
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    data <- read_csv('metadata.csv')
    data <- data %>% filter(platform == '10XGenomics')
    data <- data %>% mutate(library= str_extract(library, "(GEX|VDJ)"))
    data <- data %>% mutate(fastqSample = paste(library, sampleName, sep="_"))
    data <- data %>% mutate(fastqSample = paste(fastqSample, "demux_id", sep="/"))
    data <- data %>% mutate(fastqFolder = "${fastq_path}")
    data <- data %>% mutate(fastqPath =  paste(fastqFolder, fastqSample, sep="/${study}/"))
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    data <- data %>% mutate(molecule_h5 = paste("$params.output.folder", "Counts", sampleName, "outs/molecule_info.h5", sep="/"))
    data_gex <- data %>% filter(library == 'GEX')
    data_vdj <- data %>% filter(library == 'VDJ')
    data_gex_h5 <- data_gex %>% rename(library_id = sampleName)
    data_gex_h5 <- data_gex_h5 %>% select(library_id, molecule_h5, everything())
    data_gex_ss <- data_gex %>% select(sampleName, indices) %>% rename(Sample = sampleName, Index = indices) 
    data_gex_ss <- data_gex_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    data_vdj_ss <- data_vdj %>% select(sampleName, indices) %>% rename(Sample = sampleName, Index = indices) 
    data_vdj_ss <- data_vdj_ss %>% add_column(Lane = '1-2') %>% select(Lane, Sample, Index)
    write_csv(data_gex_h5, paste("${study}", "GEX_h5_samplesheet.csv", sep="_"))
    write_csv(data_gex_ss, paste("${study}", "GEX_samplesheet.csv", sep="_"))
    write_csv(data_vdj_ss, paste("${study}", "VDJ_samplesheet.csv", sep="_"))
    write_csv(data_vdj, paste("${study}", "VDJ_analysis_samplesheet.csv", sep="_"))
    """

}

//data_gex_h5 <- data_gex_h5 %>% select(library_id, molecule_h5, tissueType, HIVstatus)
gex_h5sheet.into { count_gex_h5sheet; aggr_gex_h5sheet; vdj_h5sheet}
process TENX_COUNT {
  echo false 

  publishDir "$params.output.folder/Counts"

  label 'gizmo_largenode'
  
  module 'cellranger'
  
  scratch "/fh/scratch/delete30/warren_h/sravisha/"

  input:
    each sample from count_gex_h5sheet.splitCsv(header: true)
    val run_count from count_gex
  
  output:
    path "${sample.library_id}" into count_path
    val task.exitStaus into count_status

  when:
    run_count == 1

  """
  cellranger count --id=$sample.library_id --transcriptome=$params.input.gex_reference \
      --fastqs=$sample.fastqPath --expect-cells=$params.count.cellcount --chemistry=$params.count.chemistry
  """
    
}

process TENX_AGGR {
  echo false

  publishDir "$params.output.folder/Counts"

  label 'gizmo_largenode'

  module 'cellranger'
  
  scratch "/fh/scratch/delete30/warren_h/sravisha/"

  input:
    path samplesheet from aggr_gex_h5sheet
    each mode from modes
    val status from count_status.collect()
    val run_aggr from aggr_gex
  
  output:
    path "Aggregate_${mode}_normalized" into aggr_path

  when:
    run_aggr == 1

  """
  cellranger aggr --id=Aggregate_${mode}_normalized --csv=${samplesheet} --normalize=${mode}
  """
}

process TENX_MATRIX {
  echo false

  publishDir "$params.output.folder/Counts/${aggr_out}/out/filtered_feature_bc_matrix"

  module 'cellranger'

  label 'gizmo'

  scratch "/fh/scratch/delete30/warren_h/sravisha/"

  input:
    path aggr_out from aggr_path
    val run_matrix from mat_gex

  output:
    path "Filtered_expression_matrix.csv" into mtx_path

  when:
    run_matrix == 1

  """
  cellranger mat2csv ${aggr_out}/outs/filtered_feature_bc_matrix Filtered_expression_matrix.csv
  """
}

process TENX_VDJ {
  echo false

  publishDir "$params.output.folder/VDJ"

  module 'cellranger'
  
  label 'gizmo_largenode'

  scratch "/fh/scratch/delete30/warren_h/sravisha/"

  input:
    val sample from vdj_analysissheet.splitCsv(header: true) 
    val run_vd from run_vdj
  output:
    path "${sample.sampleName}"

  when:
    run_vd == 1

  """
  cellranger vdj --id=$sample.sampleName --reference=$params.input.vdj_reference --fastqs=$sample.fastqFolder --sample=$sample.fastqSample
  """
}
