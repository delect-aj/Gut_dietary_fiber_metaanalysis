library(qiime2R)
library(vegan)
library(ggExtra)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(cowplot)
brewer.pal(9, "Set3")
meta <- read.delim('metadata.tsv', sep = '\t', stringsAsFactors = FALSE)
table(meta$study,meta$group)
bray_curtis_distance_matrix <- read_qza("betadistance/bray_curtis_distance_matrix.qza")
bray_curtis_distance_matrix_data <- bray_curtis_distance_matrix$data
head(bray_curtis_distance_matrix_data)
pcoa <- cmdscale(bray_curtis_distance_matrix_data, k = 2, eig = T)
pcoa_data <- data.frame(pcoa$points)
names(pcoa_data)[1:2] <- paste0("PCoA", 1:2)
pcoa_data$SampleID <- rownames(pcoa_data)
df.plot <- pcoa_data %>%
  left_join(meta[c("X.SampleID", "study", "group")], by = c("SampleID" = "X.SampleID"))
axis.1.title <- paste('PCo1 [', round((pcoa$eig[1]/sum(pcoa$eig))*100,1),'%]', sep='')
axis.2.title <- paste('PCo2 [', round((pcoa$eig[2]/sum(pcoa$eig))*100,1),'%]', sep='')
study.cols <- c("PRJEB2165" = "#8DD3C7", "PRJEB41443" = "#FFFFB3", "PRJNA293971" = "#BEBADA", "PRJNA306884" = "#FB8072", 
                "PRJNA385004" = "#80B1D3", "PRJNA428736" = "#FDB462", "PRJNA560950" = "#B3DE69", "PRJNA780023" = "#FCCDE5", 
                "PRJNA891951" = "#D9D9D9")
dune.div_study <- adonis2(bray_curtis_distance_matrix_data ~ study, data = meta, permutations = 999)
dune.div_study
dune_adonis_study <- paste0("adonis R2: ",round(dune.div_study$R2,2), "; p_value: ", dune.div_study$`Pr(>F)`)
dune_adonis_study
dune.div_group <- adonis2(bray_curtis_distance_matrix_data ~ group, data = meta, permutations = 999)
dune.div_group
dune_adonis_group <- paste0("adonis R2: ",round(dune.div_group$R2,2), "; p_value: ", dune.div_group$`Pr(>F)`)
dune_adonis_group

df.plot$group <- as.factor(df.plot$group)
g.main <- df.plot %>% 
  ggplot(aes(x=PCoA1, y=PCoA2, shape=group, col=study)) +
  geom_point() + 
  scale_colour_manual(values=study.cols) + 
  scale_shape_manual(values=c(19, 1),guide=F) + 
  scale_x_continuous(position='top') +
  xlab(axis.1.title) + ylab(axis.2.title) +
  theme(legend.position = c(0.10,0.13),panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks=element_blank(), axis.text = element_blank(),legend.text = element_text(size=10),
        panel.grid = element_blank(),legend.key = element_blank(),legend.key.size = unit(0.5, "cm"),legend.title = element_blank())
g.main
result_study_PCoA1 <- kruskal.test(PCoA1 ~ study, data = df.plot)
result_study_PCoA1$p.value
g.s.1 <- df.plot %>% 
  mutate(study=factor(study, levels=names(study.cols))) %>% 
  ggplot(aes(y=PCoA1, x=study, fill=study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study.cols, guide=FALSE) +
  xlab(paste0("Study\np = ", format(result_study_PCoA1$p.value, scientific = TRUE, digits = 2))) +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank()) + 
  coord_flip()
g.s.1
result_study_PCoA2 <- kruskal.test(PCoA2 ~ study, data = df.plot)
g.s.2 <- df.plot %>% 
  mutate(study=factor(study, levels=names(study.cols))) %>% 
  ggplot(aes(y=PCoA2, x=study, fill=study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study.cols, guide=FALSE) +
  scale_x_discrete(position='top') +
  xlab(paste0("Study\np = ", format(result_study_PCoA2$p.value, scientific = TRUE, digits = 2))) +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
g.s.2
group.cols <- c("0"="#377EB8","1"="#E41A1C")
result_group1_PCoA1 <- wilcox.test(PCoA1 ~ group, data = df.plot, alternative = "two.sided")
g.g.1 <- df.plot %>% 
  ggplot(aes(x=group, y=PCoA1, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols, guide=FALSE) + 
  ylab(axis.1.title) + 
  xlab(paste0("Group\np = ", format(result_group1_PCoA1$p.value, scientific = TRUE, digits = 2))) +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank()) + 
  coord_flip()
g.g.1
result_group1_PCoA2 <- wilcox.test(PCoA2 ~ group, data = df.plot, alternative = "two.sided")
g.g.2 <- df.plot %>% 
  ggplot(aes(x=group, y=PCoA2, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=group.cols, guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  ylab(axis.2.title) + 
  xlab(paste0("Group\np = ", format(result_group1_PCoA2$p.value, scientific = TRUE, digits = 2))) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())
g.g.2
p <- plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL, NULL, g.g.1, NULL, NULL,
          nrow=3,
          rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))
p
ggsave("Fig1B_bary.pdf",p,width = 12,height = 12)

Beita_diversity <- function(file,qza){
  meta <- read.delim(file, sep = '\t', stringsAsFactors = FALSE)
  bray_curtis_distance_matrix <- read_qza(paste0("betadistance/",qza))
  bray_curtis_distance_matrix_data <- bray_curtis_distance_matrix$data
  head(bray_curtis_distance_matrix_data)
  pcoa <- cmdscale(bray_curtis_distance_matrix_data, k = 2, eig = T)
  pcoa_data <- data.frame(pcoa$points)
  names(pcoa_data)[1:2] <- paste0("PCoA", 1:2)
  pcoa_data$SampleID <- rownames(pcoa_data)
  df.plot <- pcoa_data %>%
    left_join(meta[c("X.SampleID", "study", "group")], by = c("SampleID" = "X.SampleID"))
  axis.1.title <- paste('PCo1 [', round((pcoa$eig[1]/sum(pcoa$eig))*100,1),'%]', sep='')
  axis.2.title <- paste('PCo2 [', round((pcoa$eig[2]/sum(pcoa$eig))*100,1),'%]', sep='')
  study.cols <- c("PRJEB2165" = "#8DD3C7", "PRJEB41443" = "#FFFFB3", "PRJNA293971" = "#BEBADA", "PRJNA306884" = "#FB8072", "PRJNA385004" = "#80B1D3", "PRJNA428736" = "#FDB462", "PRJNA560950" = "#B3DE69", "PRJNA780023" = "#FCCDE5", "PRJNA891951" = "#D9D9D9")
  df.plot$group <- as.factor(df.plot$group)
  g.main <- df.plot %>% 
    ggplot(aes(x=PCoA1, y=PCoA2, shape=group, col=study)) +
    geom_point() + 
    scale_colour_manual(values=study.cols) + 
    scale_shape_manual(values=c(19, 1),guide=F) + 
    scale_x_continuous(position='top') +
    xlab(axis.1.title) + ylab(axis.2.title) +
    theme(legend.position = c(0.10,0.13),panel.background = element_rect(fill='white', color = 'black'),
          axis.ticks=element_blank(), axis.text = element_blank(),legend.text = element_text(size=10),
          panel.grid = element_blank(),legend.key = element_blank(),legend.key.size = unit(0.5, "cm"),legend.title = element_blank())
  result_study_PCoA1 <- kruskal.test(PCoA1 ~ study, data = df.plot)
  p_value <- result_study_PCoA1$p.value
  g.s.1 <- df.plot %>% 
    mutate(study=factor(study, levels=names(study.cols))) %>% 
    ggplot(aes(y=PCoA1, x=study, fill=study)) + 
    geom_boxplot() + 
    scale_fill_manual(values=study.cols, guide=FALSE) +
    xlab(paste0("Study\np = ", format(result_study_PCoA1$p.value, scientific = TRUE, digits = 2))) +
    theme(axis.ticks = element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(), axis.title.x = element_blank(),
          panel.grid = element_blank()) + 
    coord_flip()
  result_study_PCoA2 <- kruskal.test(PCoA2 ~ study, data = df.plot)
  g.s.2 <- df.plot %>% 
    mutate(study=factor(study, levels=names(study.cols))) %>% 
    ggplot(aes(y=PCoA2, x=study, fill=study)) + 
    geom_boxplot() + 
    scale_fill_manual(values=study.cols, guide=FALSE) +
    scale_x_discrete(position='top') +
    xlab(paste0("Study\np = ", format(result_study_PCoA2$p.value, scientific = TRUE, digits = 2))) +
    theme(axis.ticks=element_blank(), 
           panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank())
  group.cols <- c("0"="#377EB8","1"="#E41A1C")
  result_group1_PCoA1 <- wilcox.test(PCoA1 ~ group, data = df.plot, alternative = "two.sided")
  g.g.1 <- df.plot %>% 
    ggplot(aes(x=group, y=PCoA1, fill=group)) +
    geom_boxplot() +
    scale_fill_manual(values=group.cols, guide=FALSE) + 
    ylab(axis.1.title) + 
    xlab(paste0("Group\np = ", format(result_group1_PCoA1$p.value, scientific = TRUE, digits = 2))) +
    theme(axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          panel.background = element_rect(fill='white', color='black'),
          panel.grid = element_blank()) + 
    coord_flip()
  result_group1_PCoA2 <- wilcox.test(PCoA2 ~ group, data = df.plot, alternative = "two.sided")
  g.g.2 <- df.plot %>% 
    ggplot(aes(x=group, y=PCoA2, fill=group)) +
    geom_boxplot() +
    scale_fill_manual(values=group.cols, guide=FALSE) + 
    scale_x_discrete(position='top') + 
    scale_y_continuous(position = 'right') +
    ylab(axis.2.title) + 
    xlab(paste0("Group\np = ", format(result_group1_PCoA2$p.value, scientific = TRUE, digits = 2))) +
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          panel.background = element_rect(fill='white', color='black'),
          panel.grid = element_blank())
  p <- plot_grid(g.main, g.s.2, g.g.2, g.s.1, NULL, NULL, g.g.1, NULL, NULL,
                 nrow=3,rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))
  ggsave(paste0("Fig1B_",qza,".pdf"),p,width = 12,height = 12)
}
Beita_diversity("metadata.tsv","jaccard_distance_matrix.qza")
Beita_diversity("metadata.tsv","unweighted_unifrac_distance_matrix.qza")
Beita_diversity("metadata.tsv","weighted_unifrac_distance_matrix.qza")