library(pROC)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(forcats)   
library(purrr)    
library(tibble)   
library(qiime2R)
library(LinDA)
library(metafor)
library(readxl)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/difference_analysis/")

metadata <- read_excel("../metadata.xls") %>% column_to_rownames("Run") %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165")
abundance_table <- read.csv("../picrust2/pathways_out/path_abun_unstrat.tsv", 
                            sep="\t", row.names = 1) %>% t() %>% as.data.frame()
inter_sid <- intersect(rownames(abundance_table), rownames(metadata))
metadata <- metadata[inter_sid, ]
abundance_table <- abundance_table[inter_sid, ]
metadata$group <- factor(metadata$group)



### single study
study <- unique(metadata$study)
for (i in  study){
  metadata_pick <- metadata %>% filter(study == i)
  abundance_table_pick <- abundance_table[rownames(metadata_pick), ]
  abundance_table_filtered <- abundance_table_pick[, colSums(abundance_table_pick) != 0]
  if (i %in% c("PRJNA385004", "PRJNA560950")){
    linda_res <- linda(t(abundance_table_pick), metadata_pick, formula = '~ group', alpha = 0.05,
                       prev.cut = 0.01, lib.cut = 1000, p.adj.method = "fdr")
  } else {
    linda_res <- linda(t(abundance_table_pick), metadata_pick, formula = '~ group + (1 | subject_id)', alpha = 0.05,
                       prev.cut = 0.01, lib.cut = 1000, p.adj.method = "fdr")
  }
  res <- linda_res$output$group
  write.csv(res, paste0("linda_single_study_pathway_", i, ".csv"))
}


### combine results
df_list <- list()
for (i in study){
  df <- read.csv(paste0("linda_single_study_pathway_", i, ".csv"))
  df$study <- i
  df_list[[i]] <- df
}
results <- do.call(rbind, df_list)
rownames(results) <- NULL
r <- results %>% filter(padj < 0.1) %>% group_by(X) %>% summarize(count = n()) %>% subset(count>=3)
per_otu <- data.frame()
otu <- c()
for (i in r$X){
  pick_df <- results %>% subset(X == i)
  ma_model_1 <- rma(log2FoldChange, lfcSE, data = pick_df, method="HS")
  p <- ma_model_1$pval
  
  otu <- c(otu,i)
  per_otu_combind <- data.frame(baseMean='', log2FoldChange=ma_model_1$b[1],
                                lfcSE='', stat='',pvalue='',
                                padj=p, diff_otu=i, Study ='Combined', CI_low=ma_model_1$ci.lb, CI_hig=ma_model_1$ci.ub)
  per_otu <- rbind(per_otu, per_otu_combind)
}
write.csv(per_otu, paste0("combine_results_pathway.csv"))
