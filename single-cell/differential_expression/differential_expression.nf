#!/usr/bin/env nextflow

cds_path = Channel.fromPath(params.input.cdspath)
models = Channel.from(params.model.formula)
model_names = Channel.from(params.model.formula_name)

model_list = models.merge(model_names)

process MON_DIFF_EXP {
  echo true
  publishDir "$params.output.folder/differential_expression"
  label 'gizmo_meganode'
  scratch "/fh/scratch/delete30/warren_h/sravisha/"
  module 'R/3.6.1-foss-2016b-fh2'
  input:
    val cds from cds_path
    each model from model_list
    
  output:
    path "${params.input.study_id}_${model[1]}_model.rds" into de_model
    path "${params.input.study_id}_${model[1]}_fits_coefficients.csv" into de_results
    
  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  library(monocle3)
  library(ggplot2)
  set.seed(12357)

  cds <- readRDS("${cds}")
  gex_model <- fit_models(cds, model_formula_str = "${model[0]}", cores=10)
  cluster_fit_coefs <- coefficient_table(gex_model)
  cluster_fit_coefs <- cluster_fit_coefs %>% filter(term != "(Intercept)" & q_value < 0.05) %>%  select(gene_short_name, term, q_value, estimate)

  saveRDS(gex_model, paste(paste("${params.input.study_id}", "${model[1]}", "model", sep="_"), "rds", sep="."))
  write_csv(cluster_fit_coefs, paste(paste("${params.input.study_id}", "${model[1]}", "fits_coefficients", sep="_"), "csv", sep="."))
  """
}
