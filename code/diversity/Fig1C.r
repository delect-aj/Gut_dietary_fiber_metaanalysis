library(qiime2R)
library(vegan)
library(ggplot2)
library(tidyverse)

col <- c("#E69F00",
         "#0072B2",
         "#D6604D",
         "#7570B3")
col <- c("#E69F00", "#0072B2", "#D6604D", "#7570B3")
setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/diversity")
bc <- as.data.frame(as.matrix(read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")$data))
jaccard <- as.data.frame(as.matrix(read_qza("core-metrics-results/jaccard_distance_matrix.qza")$data))
unweightUF <- as.data.frame(as.matrix(read_qza("core-metrics-results/unweighted_unifrac_distance_matrix.qza")$data))
weightUF <- as.data.frame(as.matrix(read_qza("core-metrics-results/weighted_unifrac_distance_matrix.qza")$data))

distance_metric <- list(bc, jaccard, weightUF, unweightUF)
distance_str <- c('Bray Curtis','Jaccard','Weighted UniFrac','Unweighted UniFrac')
md <- readxl::read_xls('../metadata.xls') %>% column_to_rownames("Run")
study <- unique(md$study)

beta_stat <- function(p, dis){
  As <- list()
  s_md <- subset(md, study == p)
  sample_ids <- intersect(row.names(s_md), row.names(dis))
  s_md <- s_md[sample_ids, ]
  distance <- dis[sample_ids, sample_ids]
  adonis2_res <- adonis2(distance ~ group, data = s_md)
  return(adonis2_res)
}

Project <- c('PRJNA428736','PRJNA780023','PRJNA560950','PRJNA293971','PRJNA306884','PRJNA891951','PRJEB41443')
PerStudy_beta <- list()
for (i in 1:4){
  dis <- distance_metric[[i]]
  dis_str <- distance_str[i]
  for(p in Project){        
    PerStudy_beta[[dis_str]][[p]] <- beta_stat(p, dis)
  }
}

data <- data.frame()
for (dis_str in distance_str){
  diag_res <- data.frame(row.names=Project) 
  for (p in Project){       
    diag_res[p,1] <- PerStudy_beta[[dis_str]][[p]]['group','R2']
    diag_res[p,2] <- PerStudy_beta[[dis_str]][[p]]['group','Pr(>F)'] 
  }   
  names(diag_res) <- c('R2','pvalue')
  diag_res['study']<-row.names(diag_res)
  diag_res['metric'] <- dis_str
  data <- rbind(data, diag_res)
}

mtd <- subset(md,study == 'PRJNA428736' | study == 'PRJNA780023' | study == 'PRJNA560950'| study == 'PRJEB2165' | study == 'PRJNA293971'| study == 'PRJNA306884' | study == 'PRJNA891951' | study == 'PRJEB41443') %>% 
  subset(group == '0'| group == '1')

df <- data.frame()
for (i in 1:4){
  dis <- distance_metric[[i]]
  dis_str <- distance_str[i]
  sample_ids = intersect(row.names(mtd), row.names(dis))
  s_md <- mtd[sample_ids, ]
  distance <- as.dist(dis[sample_ids, sample_ids])
  adonis2_res <- adonis2(distance ~ group, strata = s_md$study, data = s_md) 
  R2 <- adonis2_res['group', 'R2']
  pvalue <- adonis2_res['group', 'Pr(>F)']
  l <- data.frame('R2' = R2, 'pvalue' = pvalue, 
                  'study' = 'Combined', 'metric' = dis_str)
  df <- rbind(df, l)
}
hg <- rbind(data, df) %>% mutate(Compare='Control vs. Fiber') %>% mutate(Pvalue=case_when(
  pvalue < 0.05 ~ 'P < 0.05', TRUE ~ "ns"))
write.csv(hg, "adonis_group.csv", row.names = FALSE)
