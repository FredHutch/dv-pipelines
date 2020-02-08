#!/usr/bin/env nextflow

count_dir = Channel.from(params.input.count_dir)
sample_list = Channel.from(params.input.sample_list)
nmad_val = Channel.from(params.preprocess.nmad)

process QCFilter {
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  input:
    val count_dir from count_dir
    each sample from sample_list
    val nmad from nmad_val

  output:
    path "${sample}_filtered_cds.rds" into filtered_cds

  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  ## Load in sample data
  cds <- read10xCounts(paste("${count_dir}", "${sample}", "outs", "filtered_feature_bc_matrix", sep="/"))
  cds\$Sample <- "${sample}"
  ## Quantify number of cells with low library size, gene count and high mitochondrial expression
  is_mito <- grepl(rowData(cds)\$Symbol, pattern= "^MT-")
  qc_df <- perCellQCMetrics(cds, subsets=list(mitochondrial= is_mito))
  discard <- quickPerCellQC(qc_df, percent_subsets=c("subsets_mitochondrial_percent"), nmad=${nmad}) 
  discard_summary <- discard %>% as_tibble() %>% summarize(low_lib_size = sum(low_lib_size), low_n_features = sum(low_n_features), 
                            high_subsets_mitochondrial_percent = sum(high_subsets_mitochondrial_percent), discard = sum(discard)) %>% add_column(sample = "${sample}")
  qc_tibble <- qc_df %>% as_tibble() %>% rowid_to_column()
  discard_tibble <- discard %>% as_tibble() %>% rowid_to_column() %>% dplyr::select(rowid, discard)
  qc_tibble <- left_join(qc_tibble, discard_tibble, by = "rowid") %>% add_column(sample = "${sample}")
  filtered_cds <- cds[,!discard\$discard]
  metadata(filtered_cds) <- list(raw_qc_metrics = qc_tibble, qc_fail_summary = discard_summary, nmad=${nmad})
  saveRDS(filtered_cds, paste(paste("${sample}", "filtered", "cds", sep="_"), "rds", sep="."))
  """    
}

process PlotQC {
  echo false
  publishDir "$params.output.folder/Metadata"
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from filtered_cds.collect()
  output:
    path "QC_fail_summary.csv" into discard_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(tidyverse)
  cds_list <- c("${cds_list.join('\",\"')}")
  getDiscardSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    discard_summary <- metadata(cds)\$qc_fail_summary
    return(discard_summary)
  }
  discard_summary <- cds_list %>% map(getDiscardSummary) %>% reduce(rbind)
  write_csv(discard_summary, "QC_fail_summary.csv")
  """
}