#!/usr/bin/env nextflow
//Define input channels
count_dirs = Channel.fromPath(params.input.count_dir)
sample_list = Channel.from(params.input.sample_list)
nmad_val = Channel.from(params.preprocess.nmads)
gex_metadata = Channel.from(params.processmeta.gex_metadata)
vdj_metadata = Channel.from(params.processmeta.vdj_metadata)
vdj_meta = Channel.from(params.processmeta.vdj)
//Load counts from cellranger output and filter out low quality cells
process QCFilter {
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  label 'gizmo'
  scratch "$task.scratch"
  input:
    val count_dir from count_dirs
    each sample from sample_list
    val meta_file from gex_metadata
    val nmad from nmad_val

  output:
    path "${sample}_sce.rds" into sce_obj
    path "${sample}_sce_raw.rds" into sce_raw

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
    sce <- read10xCounts(paste("${count_dir}", "${sample}", "outs", "filtered_feature_bc_matrix", sep="/"))
    sce\$Sample <- "${sample}"
    saveRDS(sce, "${sample}_sce_raw.rds")
    ## Quantify number of cells with low library size, gene count and high mitochondrial expression
    is_mito <- grepl(rowData(sce)\$Symbol, pattern= "^MT-")
    qc_df <- perCellQCMetrics(sce, subsets=list(mitochondrial= is_mito))
    discard <- quickPerCellQC(qc_df, percent_subsets=c("subsets_mitochondrial_percent"), nmad=${nmad}) 
    #discard_summary <- discard %>% as_tibble() %>% summarize(low_lib_size = sum(low_lib_size), low_n_features = sum(low_n_features), 
    #                          high_subsets_mitochondrial_percent = sum(high_subsets_mitochondrial_percent), discard = sum(discard)) %>% add_column(sample = "${sample}")
    qc_tibble <- qc_df %>% as_tibble() %>% rowid_to_column()
    discard_tibble <- discard %>% as_tibble() %>% rowid_to_column() %>% dplyr::select(rowid, discard)
    qc_tibble <- left_join(qc_tibble, discard_tibble, by = "rowid") %>% add_column(sample = "${sample}")
    ## Filter low QC cells
    sce <- sce[,!discard\$discard]
    is_mito <- grepl(rowData(sce)\$Symbol, pattern= "^MT-")
    flt_qc_df <- perCellQCMetrics(sce, subsets=list(mitochondrial= is_mito))
    ## Load metadata and grab metadata information
    meta <- read_csv("${meta_file}") %>% filter(library_id == "${sample}" & str_detect(locus, "GEX"))
    lib_info <- colData(sce) %>% as_tibble()
    lib_info <- left_join(lib_info, meta, by = c("Sample"  = "library_id"))
    ## Add metadata to CDS object for reanalysis
    metadata(sce)\$perCellQCMetrics_raw <- qc_df
    metadata(sce)\$quickPerCellQC_raw <- discard
    metadata(sce)\$perCellQCMetrics_filtered <- flt_qc_df
    metadata(sce)\$study_info <- lib_info
    metadata(sce)\$Sample <- "${sample}"
    ## Save RDS file
    saveRDS(sce, "${sample}_sce.rds")
    """
}    
//Add VDJ data when available
process VDJMetadata {
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  scratch "$task.scratch"
  label "gizmo"

  input:
    each sce from sce_obj
    //val vdj from vdj_meta
    val meta_file from vdj_metadata
    
  output:
    path "${sce.getName()}" into asce_obj

  script:
  sample = sce.getSimpleName() - ~/_\w+$/
  if( "$params.processmeta.vdj" != false )
    """
    #!/usr/bin/env Rscript
    library(scran)
    library(scater)
    library(tidyverse)
    library(DropletUtils)
    set.seed(12357)  
    sce <- readRDS("${sce}")
    vdj <- read_csv("${meta_file}") %>% filter(library_id == metadata(cds)\$sample & locus == "5primeVDJ")
    vdj_table <- colData(sce) %>% as_tibble()

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
    metadata(sce)\$doublet_barcodes <- doublet_barcode
    metadata(sce)\$vdj_table <- vdj_table
    vdj_raw <- read_csv(vdj_file)
    metadata(sce)\$vdj_raw <- vdj_raw 
    saveRDS(sce, "${sce.getName()}")
    """

  else if( "$params.processmeta.vdj" == false )
    """
    #!/usr/bin/env Rscript
    library(scran)
    library(scater)
    library(tidyverse)
    library(DropletUtils)
    set.seed(12357)  
    sce <- readRDS("${sce}")
    metadata(sce)\$doublet_barcode <- NULL
    metadata(sce)\$vdj_raw <- NULL
    saveRDS(sce, "${sce.getName()}")
    """
}
//Duplicate channels that contain annotated SCE objects
asce_obj.into{sce_write; sce_plotQC; sce_plotVDJ}
// Wtite SCE objects into different formats including monocle3, Seurat and 10X count matrix
process writeSCE{
  echo false
  module 'R/3.6.1-foss-2016b-fh2'
  scratch "$task.scratch"
  publishDir "$params.output.folder/Preprocess/CDS", mode : "move"
  label 'gizmo_largenode'
  input:
    each sce from sce_write
  output:
    path "${sample}_sce.rds" into sce_out
    path "${sample}_filtered_matrix" into mtx_obj
    path "${sample}_filtered_monocle3_cds.rds" into mon_obj
    path "${sample}_filtered_seurat_cds.rds" into seu_obj

  script:
    sample = sce.getSimpleName() - ~/_sce$/
    """
    #!/usr/bin/env Rscript
    library(scran)
    library(scater)
    library(monocle3)
    library(Seurat)
    library(tidyverse)
    library(DropletUtils)
    set.seed(12357)

    sce <- readRDS("${sce}")
    #Write SCE matrix
    saveRDS(sce, "${sample}_sce.rds")
    #Write count matrix
    write10xCounts("${sample}_filtered_matrix", counts(sce), barcodes = sce\$Barcode, 
                   gene.id=rownames(sce), gene.symbol=rowData(sce)\$Symbol)
    #Write monocle CDS
    cell_metadata <- colData(sce)
    row.names(cell_metadata) <- cell_metadata\$Barcode
    gene_metadata <- rowData(sce)
    colnames(gene_metadata) <- c("ID", "gene_short_name", "Type")
    matrix <- counts(sce)
    cds <- new_cell_data_set(matrix, cell_metadata=cell_metadata, gene_metadata=gene_metadata)
    metadata(cds) <- metadata(sce)
    saveRDS(cds, "${sample}_filtered_monocle3_cds.rds")
    #Write Seurat CDS
    counts <- assay(sce, "counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(sce) <- as.matrix(log2(t(t(counts)/size.factors) + 1))
    colnames(sce) <- sce\$Barcode
    Misc(sce, slot="perCellQCMetrics") <- metadata(sce)\$perCellQCMetrics_filtered
    Misc(sce, slot="study_info") <- metadata(sce)\$study_info
    Misc(sce, slot="doublet_barcode") <- metadata(sce)\$doublet_barcode
    Misc(sce, slot="vdj_raw") <- metadata(sce)\$vdj_raw
    Misc(sce, slot="vdj_raw_keys") <- metadata(sce)\$vdj_raw_keys
    seu_cds <- as.Seurat(sce, counts = "counts", data = "logcounts")
    saveRDS(seu_cds, "${sample}_filtered_seurat_cds.rds")
    """    
}
// Summarize and plot sample QC results from the study
process PlotQC {
  echo false
  publishDir "$params.output.folder/Preprocess/QCReports", mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    val raw_cds_list from sce_raw.collect()
    val flt_cds_list from sce_plotQC.collect()
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
    discard <- metadata(cds)\$quickPerCellQC_raw
    discard_summary <- discard %>% as_tibble() %>% 
                       summarize(low_lib_size = sum(low_lib_size), low_n_features = sum(low_n_features), 
                                 high_subsets_mitochondrial_percent = sum(high_subsets_mitochondrial_percent), discard = sum(discard)) %>% 
                       add_column(sample = metadata(cds)\$Sample)
    return(discard_summary)
  }
  discard_summary <- filtered_cds_list %>% map(getDiscardSummary) %>% reduce(rbind)
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

  getMitoPlots <- function(filtered_file) {
    flt_cds <- readRDS(filtered_file)
    sample <- metadata(flt_cds)\$Sample
    raw_cds.qc <- metadata(flt_cds)\$perCellQCMetrics_raw %>% as_tibble()
    raw_discard <- isOutlier(raw_cds.qc\$subsets_mitochondrial_percent, type="higher")
    rawplot <-  ggplot(raw_cds.qc, aes(x=`sum`, y= `subsets_mitochondrial_percent`)) + geom_point(shape=21) + coord_trans(x = "log2") + geom_hline(aes(yintercept = attr(raw_discard, "thresholds")["higher"]), colour = "red") + theme_classic() + xlab("Total count") + ylab("Mitochondrial %") + ggtitle(paste(sample, "raw")) + theme(plot.title = element_text(hjust = 0.5))
    flt_cds.qc <- metadata(flt_cds)\$perCellQCMetrics_filtered %>% as_tibble()
    flt_discard <- isOutlier(flt_cds.qc\$subsets_mitochondrial_percent, type="higher")
    fltplot <-  ggplot(flt_cds.qc, aes(x=`sum`, y= `subsets_mitochondrial_percent`)) + geom_point(shape=21) + coord_trans(x = "log2") + geom_hline(aes(yintercept = attr(flt_discard, "thresholds")["higher"]), colour = "red") + theme_classic() + xlab("Total count") + ylab("Mitochondrial %") + ggtitle(paste(sample, "filtered")) + theme(plot.title = element_text(hjust = 0.5))
    patched <- rawplot + fltplot
    return(patched)
  }
  mito_grid <-  filtered_cds_list %>% map(getMitoPlots) %>% wrap_plots() 
  ggsave("Percent_mito_plot.png", plot = mito_grid, device = "png", width = 42, height = 42, units = "cm", dpi="retina")

  getAvgCounts <- function(cds_file) {
    cds <- readRDS(cds_file)
    raw_cds_qc <- metadata(cds)\$perCellQCMetrics_raw %>% as_tibble()
    raw_median_count <- median(raw_cds_qc\$sum)
    raw_median_genes <- median(raw_cds_qc\$detected)
    raw_n_cells <- nrow(raw_cds_qc)
    flt_cds_qc <- metadata(cds)\$perCellQCMetrics_filtered %>% as_tibble()
    flt_median_count <- median(flt_cds_qc\$sum)
    flt_median_genes <- median(flt_cds_qc\$detected)
    flt_n_cells <- nrow(flt_cds_qc)
    sample <- metadata(cds)\$Sample
    raw_count_tibble = tibble(sample = sample, ncells = raw_n_cells, median_count = raw_median_count, median_genes = raw_median_genes, type = "raw")
    flt_count_tibble = tibble(sample = sample, ncells = flt_n_cells, median_count = flt_median_count, median_genes = flt_median_genes, type = "flt")
    count_tibble <- bind_rows(raw_count_tibble, flt_count_tibble)
    return(count_tibble)
  }

  study_table <- filtered_cds_list %>% map(getAvgCounts) %>% reduce(rbind)
  study_raw_tibble <- study_table %>% filter(type == "raw")
  study_flt_tibble <- study_table %>% filter(type == "flt")
  raw_count <- ggplot(study_raw_tibble, aes(x=median_count, y=ncells, color = `sample`)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median count per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw count distribution") + expand_limits(x = 0, y = 0) + theme(plot.title = element_text(hjust = 0.5))
  flt_count <- ggplot(study_flt_tibble, aes(x=median_count, y=ncells, color = `sample`)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median count per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Filtered count distribution") + expand_limits(x = 0, y = 0) + theme(plot.title = element_text(hjust = 0.5))
  raw_gene <- ggplot(study_raw_tibble, aes(x=median_genes, y=ncells, color = `sample`)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median genes expressed per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw gene expression distribution") + expand_limits(x = 0, y = 0) + theme(plot.title = element_text(hjust = 0.5))
  flt_gene <- ggplot(study_flt_tibble, aes(x=median_genes, y=ncells, color = `sample`)) + geom_point(alpha=0.8, size=3, shape="square") + coord_flip() + theme_classic() + xlab("Median genes expressed per cell") + ylab("Estimated number of cells") + labs(colour = "Sample") + ggtitle("Raw gene expression distribution") + expand_limits(x = 0, y = 0) + theme(plot.title = element_text(hjust = 0.5))
  grid_plot <- (raw_count | flt_count) / (raw_gene | flt_gene)
  ggsave("Raw_count_plot.png", plot = raw_count, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Filtered_count_plot.png", plot = flt_count, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Raw_gex_plot.png", plot = raw_gene, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Filtered_gex_plot.png", plot = flt_gene, device = "png", width = 15, height = 15, units = "cm", dpi="retina") 
  ggsave("Study_summary.png", plot =grid_plot, device = "png", width = 30, height = 30, units = "cm", dpi="retina") 
  """
}


/*process plotVDJ {
  echo false
  publishDir "$params.output.folder/Preprocess/QCReports", mode : 'move'
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    val sce from sce_plotVDJ
  output:

  script:  
    """
    #!/usr/bin/env Rscript
    library(scran)
    library(scater)
    library(ggplot2)
    library(patchwork)
    library(tidyverse)
    library(DropletUtils)
    set.seed(12357)


}*/