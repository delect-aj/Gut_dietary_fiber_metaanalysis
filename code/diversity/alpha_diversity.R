library(tidyverse)
library(qiime2R)
library(lme4)
library(lmerTest)
library(esc)  
library(metafor)
setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/diversity")

all_results_list <- list()
metadata <- read.delim('metadata.tsv', sep = '\t', stringsAsFactors = FALSE, row.names = 1) %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165")

all_study_ids <- unique(metadata$study)

for (j in c("shannon", "faith_pd", "observed_features")){
  if (j != "faith_pd"){
    alpha_df <- read_qza(paste0("core-metrics-results/", j, "_vector.qza"))$data
  } else {
    alpha_df <- read_qza(paste0("core-metrics-results/", j, "_vector.qza"))$data %>%
      column_to_rownames("V1")
  }
  colnames(alpha_df) <- c("alpha_index")
  common_samples <- intersect(rownames(metadata), rownames(alpha_df))
  metadata_j <- metadata_full[common_samples, ]
  alpha_df_j <- alpha_df[common_samples, , drop = FALSE]
  
  metadata_j$alpha_index <- alpha_df_j$alpha_index
  metadata_j <- metadata_j %>% filter(alpha_index > 0)
  study_stats_list <- list()
  for (i in all_study_ids){
    temp <- metadata_j %>% filter(study == i)
    
    p_val_ttest <- NA
    fc_simple <- NA
    es_g <- NA  
    var_g <- NA 
    
    alpha_0 <- temp[temp$group == 0, "alpha_index"]
    alpha_1 <- temp[temp$group == 1, "alpha_index"]
    
    m1 <- mean(alpha_1); sd1 <- sd(alpha_1); n1 <- length(alpha_1)
    m0 <- mean(alpha_0); sd0 <- sd(alpha_0); n0 <- length(alpha_0)
    
    if (sd1 > 0 && sd0 > 0) {
      es_calc <- esc_mean_sd(grp1m = m1, grp1sd = sd1, grp1n = n1,
                             grp2m = m0, grp2sd = sd0, grp2n = n0,
                             es.type = "g") 
      es_g <- es_calc$es    
      var_g <- es_calc$se^2 
      p_val_ttest <- t.test(alpha_1, alpha_0)$p.value
      fc_simple <- m1 / m0
    }
    
    study_stats_list[[i]] <- data.frame(
      study_id = i,
      p_value = p_val_ttest,
      fold_change = fc_simple,
      alpha_diversity = j,
      es_g = es_g,
      var_g = var_g
    )
  }
  per_study_df <- do.call(rbind, study_stats_list)

  all_results_list[[paste0(j, "_studies")]] <- per_study_df

  meta_data <- per_study_df %>% filter(!is.na(es_g) & !is.na(var_g))
  
  if (nrow(meta_data) > 1) {
    rma_res <- rma(yi = es_g, vi = var_g, data = meta_data, method = "REML")
    
    combine_df <- data.frame(
      study_id = "Combine",
      p_value = rma_res$pval,     
      fold_change = NA,           
      alpha_diversity = j,
      es_g = rma_res$b,          
      var_g = rma_res$se^2        
    )
    
    all_results_list[[paste0(j, "_combine")]] <- combine_df
  }
  
} 
df_res <- do.call(rbind, all_results_list)

write.csv(df_res, "alpha_res.csv", row.names = FALSE)
