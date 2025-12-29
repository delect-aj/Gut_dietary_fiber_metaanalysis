library(tidyverse)
library(coin)
library(pROC)
library(qiime2R)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/difference_analysis")

### OTUs
feat_all <- read_q2biom("../diversity/exported-rel-table/feature-table.biom")
p_vals_all <- read.csv("combine_results.csv", row.names = 1)

meta <- read.delim("../diversity/metadata.tsv") %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165") %>%
  column_to_rownames("X.SampleID")
inter_sid <- intersect(colnames(feat_all), rownames(meta))
feat_all <- feat_all[, inter_sid]
meta <- meta[inter_sid, ]
studies <- unique(meta$study)

auc_values <- c()
fid <- c()
study <- c()
for (f in p_vals_all$diff_otu) {
  for (s in studies){
    x <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='1') %>% rownames()]
    y <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='0') %>% rownames()]
    fid <- c(fid, f)
    study <- c(study, s)
    auc_values <- c(auc_values, c(roc(controls=y, cases=x, 
                                      direction='<', ci=TRUE, auc=TRUE)$ci)[2])
  }
}
df_auc <- data.frame(fid = fid, AUC = auc_values, study = study)
write.table(df_auc, file = "aucs.tsv",quote=FALSE, sep='\t')


### KOs
feat_all <- read.delim("../picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv", row.names = 1)
feat_all <- as.matrix(apply(feat_all, 2, function(x) x / sum(x, na.rm = TRUE)))
p_vals_all <- read.csv("combine_results_ko.csv", row.names = 1)
meta <- read.delim("../diversity/metadata.tsv") %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165") %>%
  column_to_rownames("X.SampleID")
inter_sid <- intersect(colnames(feat_all), rownames(meta))
feat_all <- feat_all[, inter_sid]
meta <- meta[inter_sid, ]
studies <- unique(meta$study)

auc_values <- c()
fid <- c()
study <- c()
for (f in p_vals_all$diff_otu) {
  for (s in studies){
    x <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='1') %>% rownames()]
    y <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='0') %>% rownames()]
    fid <- c(fid, f)
    study <- c(study, s)
    auc_values <- c(auc_values, c(pROC::roc(controls = y, cases = x, direction='<', ci=TRUE, auc=TRUE)$ci)[2])
  }
}
df_auc <- data.frame(fid = fid, AUC = auc_values, study = study)
write.table(df_auc, file = "aucs_ko.tsv",quote=FALSE, sep='\t')


### CAZY
feat_all <- read.delim("../picrust2/CAZY_metagenome_out/pred_metagenome_unstrat.tsv", row.names = 1)
feat_all <- as.matrix(apply(feat_all, 2, function(x) x / sum(x, na.rm = TRUE)))
p_vals_all <- read.csv("combine_results_cazy.csv", row.names = 1)
meta <- read.delim("../diversity/metadata.tsv") %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165") %>%
  column_to_rownames("X.SampleID")
inter_sid <- intersect(colnames(feat_all), rownames(meta))
feat_all <- feat_all[, inter_sid]
meta <- meta[inter_sid, ]
studies <- unique(meta$study)

auc_values <- c()
fid <- c()
study <- c()
for (f in p_vals_all$diff_otu) {
  for (s in studies){
    x <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='1') %>% rownames()]
    y <- feat_all[f, meta %>% filter(study == s) %>% filter(group=='0') %>% rownames()]
    fid <- c(fid, f)
    study <- c(study, s)
    auc_values <- c(auc_values, c(pROC::roc(controls = y, cases = x, direction='<', ci=TRUE, auc=TRUE)$ci)[2])
  }
}
df_auc <- data.frame(fid = fid, AUC = auc_values, study = study)
write.table(df_auc, file = "aucs_cazy.tsv",quote=FALSE, sep='\t')
