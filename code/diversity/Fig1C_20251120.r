library(qiime2R)
library(vegan)
library(ggplot2)
library(tidyverse)
library(dplyr)
biom_data <- read_q2biom("feature-table.biom")
md <- read.delim('metadata.tsv', sep = '\t', stringsAsFactors = FALSE, row.names = 1) %>%
  filter(study != "PRJNA385004" & study != "PRJEB2165")
biom_data <- biom_data[,rownames(md)]
dim(biom_data_filtered)
biom_data_filtered <- biom_data[rowSums(biom_data) != 0, ]
biom_data_filtered_t <- t(biom_data_filtered)
#dune_bray <- adonis2(data ~ group, data = md, permutations = 999, method="bray")
#dune_bray
distance_str <- c('Bray Curtis','Jaccard')
method_distance <- c("bray","jaccard")


beta_stat <- function(p, c1, c2, m){
  As <- list()
  s_md <- subset(md, study == p) %>% subset(group == c1 | group == c2)
  sample_ids <- intersect(row.names(s_md), row.names(biom_data_filtered_t))
  s_md <- s_md[sample_ids, ]
  distance <- biom_data_filtered_t[sample_ids,]
  adonis2_res <- adonis2(distance ~ group,data = s_md, permutations = 999,method = m)
  return(adonis2_res)
}
Project <- c('PRJNA428736','PRJNA780023','PRJNA560950','PRJNA293971','PRJNA306884','PRJNA891951','PRJEB41443')
PerStudy_beta <- list()
for (i in 1:2){
  dis_str <- distance_str[i]
  for(p in Project){
    PerStudy_beta[[dis_str]][[p]] <- beta_stat(p, '0', '1',method_distance[i])
  }
}
data <- data.frame()
for (dis_str in distance_str){
  diag_res <- data.frame(row.names=Project) 
  for (p in Project){       
    diag_res[p,1] <- PerStudy_beta[[dis_str]][[p]]['Model','R2']
    diag_res[p,2] <- PerStudy_beta[[dis_str]][[p]]['Model','Pr(>F)'] 
  }   
  names(diag_res) <- c('R2','pvalue')
  diag_res['study']<-row.names(diag_res)
  diag_res['metric'] <- dis_str
  data <- rbind(data, diag_res)      
}
mtd <- subset(md,study == 'PRJNA428736' | study == 'PRJNA780023' | study == 'PRJNA560950' | study == 'PRJNA293971'| study == 'PRJNA306884' | study == 'PRJNA891951' | study == 'PRJEB41443') %>%
  subset(group == '0'| group == '1')

df <- data.frame()
for (i in 1:2){
  dis_str <- distance_str[i]
  sample_ids = intersect(row.names(mtd), row.names(biom_data_filtered_t))
  s_md <- mtd[sample_ids, ]
  distance <- biom_data_filtered_t[sample_ids,]
  adonis2_res <- adonis2(distance ~ group, strata = s_md$study, data = s_md,permutations = 999,method = method_distance[i])
  R2 <- adonis2_res['Model', 'R2']
  pvalue <- adonis2_res['Model', 'Pr(>F)']
  l <- data.frame('R2' = R2, 'pvalue' = pvalue,
                  'study' = 'Combined', 'metric' = dis_str)
  df <- rbind(df, l)
}

unweightUF <- as.matrix(read.table("unifrac_unweight.dm", header = TRUE, row.names = 1))
weightUF <- as.matrix(read.table("unifrac_weight.dm", header = TRUE, row.names = 1))
unweightUF <- unweightUF[rownames(md), rownames(md)]
weightUF <- weightUF[rownames(md), rownames(md)]
distance_metric1 <- list(weightUF, unweightUF)
distance_str1 <- c('Weighted UniFrac','Unweighted UniFrac')
study <- unique(md$study)
beta_stat1 <- function(p, dis , c1, c2){
  As <- list()
  s_md <- subset(md, study == p) %>% subset(group == c1 | group == c2)
  sample_ids <- intersect(row.names(s_md), row.names(dis))
  s_md <- s_md[sample_ids, ]
  distance <- as.dist(dis[sample_ids, sample_ids])
  adonis2_res <- adonis2(distance ~ group, permutations = 999,data = s_md)
  return(adonis2_res)
}
Project <- c('PRJNA428736','PRJNA780023','PRJNA560950','PRJNA293971','PRJNA306884','PRJNA891951','PRJEB41443')
PerStudy_beta1 <- list()
for (i in 1:2){
  dis <- distance_metric1[[i]]
  dis_str <- distance_str1[i]
  for(p in Project){
    PerStudy_beta1[[dis_str]][[p]] <- beta_stat1(p, dis, '0', '1')
  }
}

data1 <- data.frame()
for (dis_str in distance_str1){
  diag_res <- data.frame(row.names=Project) 
  for (p in Project){       
    diag_res[p,1] <- PerStudy_beta1[[dis_str]][[p]]['Model','R2']
    diag_res[p,2] <- PerStudy_beta1[[dis_str]][[p]]['Model','Pr(>F)'] 
  }   
  names(diag_res) <- c('R2','pvalue')
  diag_res['study']<-row.names(diag_res)
  diag_res['metric'] <- dis_str
  data1 <- rbind(data1, diag_res)      
}

mtd <- subset(md,study == 'PRJNA428736' | study == 'PRJNA780023' | study == 'PRJNA560950' | study == 'PRJNA293971'| study == 'PRJNA306884' | study == 'PRJNA891951' | study == 'PRJEB41443') %>%
  subset(group == '0'| group == '1')

df1 <- data.frame()
for (i in 1:2){
  dis <- distance_metric1[[i]]
  dis_str <- distance_str1[i]
  sample_ids = intersect(row.names(mtd), row.names(dis))
  s_md <- mtd[sample_ids, ]
  distance <- as.dist(dis[sample_ids, sample_ids])
  adonis2_res <- adonis2(distance ~ group, strata = s_md$study, data = s_md)
  R2 <- adonis2_res['Model', 'R2']
  pvalue <- adonis2_res['Model', 'Pr(>F)']
  l <- data.frame('R2' = R2, 'pvalue' = pvalue,
                  'study' = 'Combined', 'metric' = dis_str)
  df1 <- rbind(df1, l)
}
hg <- bind_rows(data, data1, df, df1) %>% mutate(Compare='Control vs. Fiber') %>% mutate(Pvalue=case_when(
  pvalue < 0.05 ~ 'P < 0.05', TRUE ~ "ns"))
hg$metric <- as.factor(hg$metric)
col <- c("#E69F00", "#0072B2", "#D6604D", "#7570B3")
p_2 <- ggplot(data = hg, aes(x = R2, y = study, color = metric))+
  geom_point(aes(shape =Pvalue),size=2)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10,hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.margin = unit(rep(2,4), 'lines')) + scale_color_manual(values = col)+
  labs(title = 'Control vs Fiber', y = '')
p_2
ggsave("Figure_1C_new2_20251119.pdf", p_2, width = 4.5, height = 3.5)
