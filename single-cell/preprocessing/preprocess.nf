#!/usr/bin/env nextflow

count_dirs = Channel.fromPath(params.input.count_dir)
sample_list = Channel.from(params.input.sample_list)
nmad_val = Channel.from(params.preprocess.nmads)
gex_metadata = Channel.from(params.processmeta.gex_metadata)
vdj_metadata = Channel.from(params.processmeta.vdj_metadata)
vdj_meta = Channel.from(params.processmeta.vdj)

process QCFilter {
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/Preprocess/CDS" , mode : 'copy'
  label 'gizmo'
  input:
    val count_dir from count_dirs
    each sample from sample_list
    val meta_file from gex_metadata
    val nmad from nmad_val

  output:
    path "${sample}_${uuid}_filtered_cds.rds" into filtered_cds
    path "${sample}_${uuid}_filtered_monocle3_cds.rds" into mon_filtered_cds
    path "${sample}_${uuid}_filtered_seurat_cds.rds" into seu_filtered_cds
    path "${sample}_${uuid}_raw_cds.rds" into raw_cds
    path "${sample}_${uuid}_filtered_matrix" into matrix_path

  script:
    uuid = UUID.randomUUID().toString().substring(0,7)
    """
    #!/usr/bin/env Rscript
    library(scran)
    library(scater)
    library(monocle3)
    library(Seurat)
    library(tidyverse)
    library(DropletUtils)
    set.seed(12357)
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
    ## Write raw CDS object for reanalysis
    metadata(cds) <- list(quickPerCellQC = discard, perCellQCMetrics = qc_df, qc_fail_summary = discard_summary, Sample = "${sample}")
    saveRDS(cds, paste(paste("${sample}", "${uuid}", "raw", "cds", sep="_"), "rds", sep="."))
    ## Filter low QC cells
    filtered_cds <- cds[,!discard\$discard]
    is_mito <- grepl(rowData(filtered_cds)\$Symbol, pattern= "^MT-")
    filtered_qc_metrics <- perCellQCMetrics(filtered_cds, subsets=list(mitochondrial= is_mito)) %>% as_tibble() %>% rowid_to_column()
    metadata(filtered_cds)\$run_info <- list(uuid = "${uuid}", Sample = "${sample}", cds_file = paste(paste("${sample}", "${uuid}", "filtered", "cds", sep="_"), "rds", sep="."), 
                                    raw_qc_metrics = qc_tibble, qc_fail_summary = discard_summary, nmad=${nmad})
    metadata(filtered_cds)\$perCellQCMetrics <- filtered_qc_metrics
    ## Load metadata and grab metadata information
    meta <- read_csv("${meta_file}") %>% filter(library_id == "${sample}" & locus == "5primeGEX")
    lib_info <- colData(cds) %>% as_tibble()
    lib_info <- left_join(lib_info, meta, by = c("Sample"  = "library_id"))
    metadata(filtered_cds)\$study_info <- lib_info
    saveRDS(filtered_cds, paste(paste("${sample}", "${uuid}", "filtered", "cds", sep="_"), "rds", sep="."))
    write10xCounts(paste("${sample}", "${uuid}", "filtered", "matrix", sep="_"), counts(filtered_cds), barcodes = filtered_cds\$Barcode, gene.id=rownames(filtered_cds), gene.symbol=rowData(filtered_cds)\$Symbol)

    #Write monocle CDS
    cell_metadata <- colData(filtered_cds)
    row.names(cell_metadata) <- cell_metadata\$Barcode
    gene_metadata <- rowData(filtered_cds)
    colnames(gene_metadata) <- c("ID", "gene_short_name", "Type")
    matrix <- counts(filtered_cds)
    mon_cds <- new_cell_data_set(matrix, cell_metadata=cell_metadata, gene_metadata=gene_metadata)
    metadata(mon_cds) <- metadata(filtered_cds)
    saveRDS(mon_cds, paste(paste("${sample}", "${uuid}", "filtered", "monocle3", "cds", sep="_"), "rds", sep="."))

    #Write Seurat CDS
    counts <- assay(filtered_cds, "counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(filtered_cds) <- as.matrix(log2(t(t(counts)/size.factors) + 1))
    colnames(filtered_cds) <- filtered_cds\$Barcode
    seu_cds <- as.Seurat(filtered_cds, counts = "counts", data = "logcounts")
    saveRDS(seu_cds, paste(paste("${sample}", "${uuid}", "filtered", "seurat", "cds", sep="_"), "rds", sep="."))
    """    
}

filtered_cds.into{pqc_filtered_cds; vdj_filtered_cds}

process PlotQC {
  echo false
  publishDir "$params.output.folder/Preprocess/QCReports", mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val raw_cds_list from raw_cds.collect()
    val flt_cds_list from pqc_filtered_cds.collect()
  output:
    path "QC_fail_summary.csv" into discard_report
    path "Knee_plot.png" into knee_grid
    path "Percent_mito_plot.png" into mito_grid
    path "Raw_count_plot.png" into raw_count
    path "Filtered_count_plot.png" into flt_count
    path "Raw_gex_plot.png" into raw_gex
    path "Filtered_gex_plot.png" into flt_gex
    path "Study_summary.png" into study_summary
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  raw_cds_list <- c("${raw_cds_list.join('\",\"')}")
  filtered_cds_list <- c("${flt_cds_list.join('\",\"')}")
  getDiscardSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    discard_summary <- metadata(cds)\$qc_fail_summary
    return(discard_summary)
  }
  discard_summary <- raw_cds_list %>% map(getDiscardSummary) %>% reduce(rbind)
  write_csv(discard_summary, "QC_fail_summary.csv")

  getKneePlots <- function(cds_file) {
    cds <- readRDS(cds_file)
    sample <- metadata(cds)\$Sample
    bcrank <- barcodeRanks(counts(cds))
    uniq <- !duplicated(bcrank\$rank)
    kneeplot <- ggplot(as_tibble(bcrank[uniq, ]), aes(x = rank, y = total)) + geom_point(shape = 21) + coord_trans(x="log2", y="log2") + geom_hline(aes(yintercept = metadata(bcrank)\$inflection, linetype = "Inflection"),  colour = "darkgreen") + geom_hline(aes(yintercept = metadata(bcrank)\$knee, linetype = "Knee"),  colour = "dodgerblue") + theme_classic() + scale_linetype_manual(name = "Threshold", values = c(2, 2), guide = guide_legend(override.aes = list(color = c("darkgreen", "dodgerblue")))) + xlab("Rank") + ylab("Total UMI count") + ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5))
    return(kneeplot)
  }
  knee_grid <- raw_cds_list %>% map(getKneePlots) %>% wrap_plots() + plot_layout(guides = 'collect')
  ggsave("Knee_plot.png", plot = knee_grid, device = "png", width = 42, height = 42, units = "cm", dpi="retina")

  getMitoPlots <- function(raw_file, filtered_file) {
    raw_cds <- readRDS(raw_file)
    flt_cds <- readRDS(filtered_file)
    sample <- metadata(raw_cds)\$Sample
    raw_cds.qc <- metadata(raw_cds)\$perCellQCMetrics %>% as_tibble()
    raw_discard <- isOutlier(raw_cds.qc\$subsets_mitochondrial_percent, type="higher")
    rawplot <-  ggplot(raw_cds.qc, aes(x=`sum`, y= `subsets_mitochondrial_percent`)) + geom_point(shape=21) + coord_trans(x = "log2") + geom_hline(aes(yintercept = attr(raw_discard, "thresholds")["higher"]), colour = "red") + theme_classic() + xlab("Total count") + ylab("Mitochondrial %") + ggtitle(paste(sample, "raw")) + theme(plot.title = element_text(hjust = 0.5))
    flt_cds.qc <- metadata(flt_cds)\$perCellQCMetrics %>% as_tibble()
    flt_discard <- isOutlier(flt_cds.qc\$subsets_mitochondrial_percent, type="higher")
    fltplot <-  ggplot(flt_cds.qc, aes(x=`sum`, y= `subsets_mitochondrial_percent`)) + geom_point(shape=21) + coord_trans(x = "log2") + geom_hline(aes(yintercept = attr(flt_discard, "thresholds")["higher"]), colour = "red") + theme_classic() + xlab("Total count") + ylab("Mitochondrial %") + ggtitle(paste(sample, "filtered")) + theme(plot.title = element_text(hjust = 0.5))
    patched <- rawplot + fltplot
    return(patched)
  }
  
  cds_tibble <- tibble(raw = raw_cds_list, filtered = filtered_cds_list)
  mito_grid <-  map2(cds_tibble\$raw, cds_tibble\$filtered, getMitoPlots) %>% wrap_plots() 
  ggsave("Percent_mito_plot.png", plot = mito_grid, device = "png", width = 42, height = 42, units = "cm", dpi="retina")

  getAvgCounts <- function(cds_file) {
    cds <- readRDS(cds_file)
    cds_qc <- metadata(cds)\$perCellQCMetrics %>% as_tibble()
    median_count <- median(cds_qc\$sum)
    median_genes <- median(cds_qc\$detected)
    n_cells <- nrow(cds_qc)
    sample <- metadata(cds)\$Sample
    count_tibble = tibble(sample = sample, ncells = n_cells, median_count = median_count, median_genes = median_genes)
    return(count_tibble)
  }

  study_raw_tibble <- raw_cds_list %>% map(getAvgCounts) %>% reduce(rbind)
  study_flt_tibble <- filtered_cds_list %>% map(getAvgCounts) %>% reduce(rbind)
  raw_count <- ggplot(study_raw_tibble, aes(x=median_count, y=ncells, color = sample)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median count per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw count distribution") + xlim(0,5000) + ylim(0,5000) + theme(plot.title = element_text(hjust = 0.5))
  flt_count <- ggplot(study_flt_tibble, aes(x=median_count, y=ncells, color = sample)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median count per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Filtered count distribution") + xlim(0,5000) + ylim(0,5000) + theme(plot.title = element_text(hjust = 0.5))
  raw_gene <- ggplot(study_raw_tibble, aes(x=median_genes, y=ncells, color = sample)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median genes expressed per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw gene expression distribution") + xlim(0,5000) + ylim(0,5000) + theme(plot.title = element_text(hjust = 0.5))
  flt_gene <- ggplot(study_flt_tibble, aes(x=median_genes, y=ncells, color = sample)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median genes expressed per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw gene expression distribution") + xlim(0,5000) + ylim(0,5000) + theme(plot.title = element_text(hjust = 0.5))
  grid_plot <- (raw_count | flt_count) / (raw_gene | flt_gene)
  ggsave("Raw_count_plot.png", plot = raw_count, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Filtered_count_plot.png", plot = flt_count, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Raw_gex_plot.png", plot = raw_gene, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Filtered_gex_plot.png", plot = flt_gene, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Study_summary.png", plot =grid_plot, device = "png", width = 30, height = 30, units = "cm", dpi="retina") 
  """
}

process VDJMetadata {
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  scratch '/fh/scratch/delete30/warren_h/sravisha'
  publishDir "$params.output.folder/Preprocess/CDS", mode  : 'copy'
  
  when:
    vdj_data == true
  input:
    path cds_path from vdj_filtered_cds.mix(mon_filtered_cds)
    val meta_file from vdj_meta
  
  output:
    
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)  
  cds <- readRDS("${cds_path}")
  vdj <- read_csv("${meta_file}") %>% filter(library_id == metadata(cds)\$sample & locus == "5primeVDJ")
  vdj_table <- colData(cds) %>% as_tibble()

  collapseClonotype <- function(clonotypes, vdj_type) {
    if (vdj_type == NULL) {
    clonotypes <- clonotypes %>% filter(high_confidence == TRUE & is_cell == TRUE & productive == TRUE) %>% 
            group_by(barcode) %>% summarize(is_cell = first(is_cell), contig_id = paste(contig_id, collapse = ";"), 
                                            high_confidence = first(high_confidence), `length` = paste(`length`, collapse = ";"),
                                            chain = paste(chain, collapse = ";"), v_gene = paste(v_gene, collapse = ";"),
                                            j_gene = paste(j_gene, collapse = ";"), c_gene = paste(c_gene, collapse = ";"),
                                            full_length = first(full_length), productive = first(productive),
                                            cdr3 = paste(cdr3, collapse = ";"), cdr3_nt = paste(cdr3_nt, collapse=";"),
                                            reads = paste(reads, collapse = ";"), umis = paste(umis, collapse = ";"), 
                                            raw_clonotype_id = first(raw_clonotype_id), cell_type = paste(cell_type, collapse = ";"))
    } else {
    clonotypes <- clonotypes %>% filter(high_confidence == TRUE & is_cell == TRUE & productive == TRUE) %>% 
            group_by(barcode) %>% summarize(is_cell = first(is_cell), contig_id = paste(contig_id, collapse = ";"), 
                                            high_confidence = first(high_confidence), `length` = paste(`length`, collapse = ";"),
                                            chain = paste(chain, collapse = ";"), v_gene = paste(v_gene, collapse = ";"),
                                            j_gene = paste(j_gene, collapse = ";"), c_gene = paste(c_gene, collapse = ";"),
                                            full_length = first(full_length), productive = first(productive),
                                            cdr3 = paste(cdr3, collapse = ";"), cdr3_nt = paste(cdr3_nt, collapse=";"),
                                            reads = paste(reads, collapse = ";"), umis = paste(umis, collapse = ";"), 
                                            raw_clonotype_id = first(raw_clonotype_id)) %>% add_column(cell_type = vdj_type)
    }
    return(clonotypes)
  }

  mergeVDJ <- function(vdj_file, vdj_type) {
    clonotypes <- read_csv(vdj_file)
    clonotypes <- collapseClonotype(clonotypes, vdj_type)
    return(clonotypes)
  }
  clonotypes <- map2(vdj\$vdj_sequences, vdj\$VDJType, mergeVDJ) %>% reduce(rbind) %>% collapseClonotype(,vdj_type=NULL)
  vdj_table <- left_join(vdj_table, clonotypes, by=c("Barcode" = "barcode"))
  doublet_barcode <- vdj_table %>% filter(str_detect(cell_type, "B cell") & str_detect(cell_type, "T Cell")) %>% select(Barcode) %>% add_column(Method = "VDJ data")
  metadata(cds)\$doublet_barcodes <- doublet_barcode
  cds\$vdj_table <- vdj_table
  vdj_raw <- read_csv(vdj_file)
  cds\$vdj_raw <- vdj_raw %>% group_split(chain)
  cds\$vdj_raw_keys <- vdj_raw %>% group_keys(chain)
  saveRDS(cds, ${cds_path})
  """
  
}
