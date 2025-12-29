library(tidyverse)
library(qiime2R)
library(lme4)
library(vegan)
library(lmerTest)
library(esc)  
library(readxl)
library(metafor)
setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/diversity")

metadata <- read_excel("../metadata.xls") %>% column_to_rownames("Run") %>%
  filter(study != "PRJNA385004") %>% filter(study != "PRJEB2165")
r2 <- c()
p_value <- c()
for (i in c("bray_curtis", "jaccard", "weighted_unifrac", "unweighted_unifrac")){
  distance_matrix <- read_qza(paste0("core-metrics-results/", i, "_distance_matrix.qza"))
  distance_matrix <- as.matrix(distance_matrix$data)
  inter_id <- intersect(rownames(metadata), rownames(distance_matrix))
  metadata <- metadata[inter_id, ]
  distance_matrix <- as.dist(distance_matrix[inter_id, inter_id])
  adonis_res <- adonis2(distance_matrix ~ study + group + intervention, 
                        data = metadata, permutations = 999, by="terms")
  adonis_res <- as.data.frame(adonis_res)
  r2 <- c(r2, adonis_res$R2[1:3])
  p_value <- c(p_value, adonis_res$`Pr(>F)`[1:3])
}

adonis_res <- data.frame(r2 = r2, p_value = p_value, 
                         group = rep(c("study", "group", "fiber type"), 4),
                         method = rep(c("bray_curtis", "jaccard", "weighted_unifrac", "unweighted_unifrac"), each=3))

adonis_res$group <- factor(adonis_res$group, levels = c("study", "group", "fiber type"))

p <- ggplot(adonis_res, aes(x = group, y = r2)) +
  geom_col() +
  geom_text(aes(label = paste0(round(r2, 3), "*")),
            vjust = -0.5, size = 3) +
  facet_wrap(~method, nrow = 1) +
  labs(x = '', y = 'R2') +
  theme_bw()

ggsave("supplementary_adonis.pdf", p, width = 8, height = 4)

