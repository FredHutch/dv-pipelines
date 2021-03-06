#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cds_path + "*/*_filtered_cds.rds")
adownsample = Channel.from(params.osca.abinitio.nsample)
anfeatures = Channel.from(params.osca.abinitio.nfeatures)
acenters = Channel.from(params.osca.abinitio.centers)
pdownsample = Channel.from(params.osca.priori.nsample)
pcenters = Channel.from(params.osca.priori.centers)


cds_path.into{ab_cds_path; pr_cds_path}

process OSCA_ABINTIO_CLUSTER {

  echo false
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  publishDir "$params.output.folder/OSCA/Cluster/Abinitio", mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each cds from ab_cds_path
    val nsample from adownsample
    each nfeature from anfeatures
    each cnumber from acenters

  output:
    file "Abinitio_${uuid}.rds" into abinitio_cds
    file "Abinitio_${uuid}_tSNE.png" into abinitio_tsne
    file "Abinitio_${uuid}_UMAP.png" into abinitio_umap
    
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(leiden)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    variable_genes <- getTopHVGs(gene_var, n=${nfeature})
    cds <- runPCA(cds, subset_row=variable_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(pcs)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metadata(choices)\$chosen]
    cds <- runTSNE(cds, dimred="PCA_chosen") 
    cds <- runUMAP(cds, dimred="PCA_chosen")
    cluster <- buildSNNGraph(cds, k=${cnumber}, use.dimred = 'PCA')
    cluster_adj <- igraph::as_adjacency_matrix(cluster)
    cluster_leiden <- leiden(cluster_adj)
    cluster_louvain <- igraph::cluster_louvain(cluster)
    #cds\$cluster <- factor(cluster_louvain\$membership)
    cds\$cluster <- factor(cluster_leiden)
    abi_run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Abinitio", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = metadata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "UMAP/tSNE", 
                            k = ${cnumber}, cluster_method = "louvain", sample = "${sample}")
    metadata(cds) <- list(run_condition = abi_run_condition)
    saveRDS(cds, paste(paste("Abinitio", "${uuid}", sep="_"), "rds", sep="."))
    abinitio_tsne <- plotReducedDim(cds, dimred="TSNE", colour_by="cluster")
    abinitio_umap <- plotReducedDim(cds, dimred="UMAP", colour_by="cluster")
    ggsave(paste(paste("Abinitio", "${uuid}", "tSNE", sep="_"), "png", sep="."), abinitio_tsne, units = "in", dpi = "retina")
    ggsave(paste(paste("Abinitio", "${uuid}", "UMAP", sep="_"), "png", sep="."), abinitio_umap, units = "in", dpi = "retina")
    """
}

process ABI_GATHER_METADATA {

  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Abinitio" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from abinitio_cds.collect()
  output:
    path "Abinitio_cluster_metadata.csv" into abi_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  cds_list <- c("${cds_list.join('\",\"')}")
  getClusterSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    cluster_summary <- metadata(cds)\$run_condition
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Abinitio_cluster_metadata.csv")
  """
}


process OSCA_PRIORI_CLUSTER {
  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Priori"
  module 'R/3.6.1-foss-2016b-fh2'
  module 'Python/3.7.4-foss-2016b'
  input:
    each cds from pr_cds_path
    val nsample from pdownsample
    each cnumber from pcenters

  output:
    file "Priori_${uuid}.rds" into priori_cds
    file "Priori_${uuid}_tSNE.png" into priori_tsne
    file "Priori_${uuid}_UMAP.png" into priori_umap
  
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      
    sample = cds.getSimpleName()
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    library(leiden)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(pcs)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metadata(choices)\$chosen]
    cds <- runTSNE(cds, dimred="PCA_chosen") 
    cds <- runUMAP(cds, dimred="PCA_chosen")
    cluster <- buildSNNGraph(cds, k=${cnumber}, use.dimred = 'PCA')
    cluster_adj <- igraph::as_adjacency_matrix(cluster)
    cluster_leiden <- leiden(cluster_adj)
    cluster_louvain <- igraph::cluster_louvain(cluster)
    #cds\$cluster <- factor(cluster_louvain\$membership)
    cds\$cluster <- factor(cluster_leiden)
    pri_run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = metadata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "UMAP/tSNE", 
                            k = ${cnumber}, cluster_method = "louvain", sample = "${sample}")
    metadata(cds) <- list(run_condition = pri_run_condition)
    saveRDS(cds, paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."))
    priori_tsne <- plotReducedDim(cds, dimred="TSNE", colour_by="cluster")
    priori_umap <- plotReducedDim(cds, dimred="UMAP", colour_by="cluster")
    ggsave(paste(paste("Priori", "${uuid}", "tSNE", sep="_"), "png", sep="."), priori_tsne, units = "in", dpi = "retina")
    ggsave(paste(paste("Priori", "${uuid}", "UMAP", sep="_"), "png", sep="."), priori_umap, units = "in", dpi = "retina")
    """

}

process PRI_GATHER_METADATA {

  echo false
  publishDir "$params.output.folder/OSCA/Cluster/Priori" , mode : 'copy'
  module 'R/3.6.1-foss-2016b-fh2'

  input:
    val cds_list from priori_cds.collect()
  output:
    path "Priori_cluster_metadata.csv" into pri_report
  
  """
  #!/usr/bin/env Rscript
  library(scran)
  library(scater)
  library(tidyverse)
  library(DropletUtils)
  set.seed(12357)
  cds_list <- c("${cds_list.join('\",\"')}")
  getClusterSummary <- function(cds_file) {
    cds <- readRDS(cds_file)
    cluster_summary <- metadata(cds)\$run_condition
    return(cluster_summary)
  }
  cluster_summary <- cds_list %>% map(getClusterSummary) %>% reduce(rbind)
  write_csv(cluster_summary, "Priori_cluster_metadata.csv")
  """
}


/*process GAT_UMAP_CLUSTER {
  echo false
  publishDir "$params.output.folder/Gattardo/Cluster/UMAPOpt"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    each cds from um_cds_path
    val nsample from udownsample


  output:
    file "Monocle_${uuid}.rds" into monocle_cds
    
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$gs_name\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(cds)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metdata(choices)\$chosen]
    cds <- runUMAP(cds, dimred="PCA_chosen", n_neighbors = ${uneighbor}, min_dist = ${udist})
    cluster <- buildSNNGraph(cds, k=${cnumber}, use.dimred = 'PCA')
    cluster_louvain <- igraph::cluster_louvain(cluster)
    cds\$cluster <- factor(cluster_louvain)
    run_condition <- tibble(uuid = "${uuid}", cds_file = paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."), preprocess_num_dim = metdata(choices)\$chosen, 
                            preprocess_norm_method = "Log normal", preprocess_method = "PCA", 
                            reduction_method = "UMAP/tSNE", 
                            k = ${cnumber}, cluster_method = "louvain", sample = cell_metadata\$Sample[1])
    metadata(cds) <- list(run_condition = run_condition)
    saveRDS(cds, paste(paste("Priori", "${uuid}", sep="_"), "rds", sep="."))
    priori_tsne <- plotReducedDim(cds, dimred="TSNE", colour_by="cluster")
    priori_umap <- plotReducedDim(cds, dimred="UMAP", colour_by="cluster")
    ggsave(paste(paste("Priori", "${uuid}", "tSNE", sep="_"), "png", sep="."), priori_tsne, units = "in", dpi = "retina")
    ggsave(paste(paste("Priori", "${uuid}", "UMAP", sep="_"), "png", sep="."), priori_umap, units = "in", dpi = "retina")
    """

}

process GAT_TSNE_CLUSTER {
  echo false
  publishDir "$params.output.folder/Gattardo/Cluster/TSNEOpt"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    each cds from ts_cds_path
    val core from num_cores
    each pdim from num_dim
    val pnorm from norm_method
    val preduce from dimred_method
    each rreduce from reduce_method
    each cmethod from clustering_method
    each cnumber from cluster_number

  output:
    file "Monocle_${uuid}.rds" into monocle_cds
    
  script:
    uuid = UUID.randomUUID().toString().substring(0,7)      

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(scater)
    library(scran)
    library(msigdbr)
    set.seed(12357)

    cds <- readRDS("${cds}")
    cluster_cds <- quickCluster(cds)
    cds <- computeSumFactors(cds, cluster=cluster_cds, min.mean=0.1)
    cds <- logNormCounts(cds, downsample=${nsample})
    gene_var <- modelGeneVar(cds)
    #gene_var_cv2 <- modelGeneCV2(cds)
    immuno_sig <- msigdbr(species = "Homo sapiens", category = "C7")
    immuno_sig_genes <- rowData(cds)\$Symbol %in% immuno_sig\$gs_name\$human_gene_symbol
    cds <- runPCA(cds, subset_row=immuno_sig_genes)
    pcs <- reducedDim(cds)
    choices <- getClusteredPCs(cds)
    reducedDim(cds, "PCA_chosen") <- reducedDim(cds, "PCA")[,1:metdata(choices)\$chosen]
    cds <- runTSNE(cds, dimred="PCA_chosen", perplexity = ${perplexity}) 
    """

}*/
