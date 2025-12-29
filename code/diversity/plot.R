library(qiime2R)
library(vegan)
library(ggExtra)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/diversity")

### PCoA
meta <- read.delim('metadata.tsv', sep = '\t', stringsAsFactors = FALSE, row.names = 1) %>%
  filter(study != "PRJNA385004")
distance_matrix <- as.matrix(read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")$data)
inter_id <- intersect(rownames(meta), rownames(distance_matrix))
distance_matrix <- as.dist(distance_matrix[inter_id, inter_id])

meta <- meta[inter_id, ]
# adonis2
adonis_res <- adonis2(distance_matrix ~ study + group, data = meta, permutations = 999, by = "terms")
adonis_res <- as.data.frame(adonis_res)
adonis_study <- paste0("Study R2: ", round(adonis_res$R2[1], 4), " ", if_else(adonis_res$`Pr(>F)`[1] < 0.05, "*", ""))
adonis_group <- paste0("Group R2: ", round(adonis_res$R2[2], 4), " ", if_else(adonis_res$`Pr(>F)`[2] < 0.05, "*", ""))

pcoa <- cmdscale(distance_matrix, k = 2, eig = T)
pcoa_data <- data.frame(pcoa$points)
colnames(pcoa_data) <- paste0("PC", 1:2)
plot_data <- cbind(pcoa_data, meta[rownames(pcoa_data), ])

axis_1_title <- paste('PCoA1 [', round((pcoa$eig[1]/sum(pcoa$eig))*100,1),'%]', sep='')
axis_2_title <- paste('PCoA2 [', round((pcoa$eig[2]/sum(pcoa$eig))*100,1),'%]', sep='')
study_cols <- c("PRJEB2165" = "#8DD3C7", "PRJEB41443" = "#FFFFB3", "PRJNA293971" = "#BEBADA", 
                "PRJNA306884" = "#FB8072", "PRJNA385004" = "#80B1D3", "PRJNA428736" = "#FDB462", 
                "PRJNA560950" = "#B3DE69", "PRJNA780023" = "#FCCDE5", "PRJNA891951" = "#D9D9D9")
plot_data$group <- as.factor(plot_data$group)
g_main <- plot_data %>% 
  ggplot(aes(x = PC1, y = PC2, shape = group, col = study)) +
  geom_point() + 
  scale_colour_manual(name = "Study",values = study_cols, guide=FALSE) + 
  scale_shape_manual(name = "Group",values = c(19, 1), labels = c("Control", "Fiber"), guide=FALSE) + 
  scale_x_continuous(position = 'top') +
  xlab(paste0(axis_1_title, "\n")) + ylab(paste0(axis_2_title, "\n")) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste(adonis_study, adonis_group, sep = "\n"),
           hjust = -0.05, vjust = -0.5, 
           size = 3, color = "black") +
  theme_bw() +
  theme(panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), legend.text = element_text(size=10),
        panel.grid = element_blank(), legend.key = element_blank(), 
        legend.key.size = unit(0.5, "cm"), legend.title = element_blank()
  )
result_study_PCoA1 <- kruskal.test(PC1 ~ study, data = plot_data)
if (result_study_PCoA1$p.value < 2.2e-16){
  x_lab <- "Study\np<2.2e-10"
} else {
  x_lab <- paste0("Group\np=", format(result_study_PCoA1$p.value, scientific = TRUE, digits = 2))
}
g_s_1 <- plot_data %>% 
  mutate(study=factor(study, levels = names(study_cols))) %>% 
  ggplot(aes(y = PC1, x = study, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study_cols, guide=FALSE) +
  xlab(x_lab) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank()) + 
  coord_flip()
result_study_PCoA2 <- kruskal.test(PC2 ~ study, data = plot_data)
if (result_study_PCoA2$p.value < 2.2e-16){
  x_lab <- "Study\np<2.2e-10"
} else {
  x_lab <- paste0("Group\np=", format(result_study_PCoA2$p.value, scientific = TRUE, digits = 2))
}
g_s_2 <- plot_data %>% 
  mutate(study = factor(study, levels = names(study_cols))) %>% 
  ggplot(aes(y = PC2, x = study, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values = study_cols, guide = FALSE) +
  scale_x_discrete(position='top') +
  xlab(x_lab) +
  theme_bw() +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
group_cols <- c("0" = "#377EB8", "1" = "#E41A1C")
result_group_PCoA1 <- wilcox.test(PC1 ~ group, data = plot_data, alternative = "two.sided")
if (result_group_PCoA1$p.value < 2.2e-16){
  x_lab <- "Study\np<2.2e-10"
} else {
  x_lab <- paste0("Group\np=", format(result_group_PCoA1$p.value, scientific = TRUE, digits = 2))
}
g_g_1 <- plot_data %>% 
  ggplot(aes(x=group, y=PC1, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(name = "Type",values=group_cols, labels = c("Control", "Fiber"), guide=FALSE) + 
  xlab(x_lab) +
  ylab("PCoA1") +
  theme_bw() +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank(),legend.key = element_blank()) + 
  coord_flip()
result_group_PCoA2 <- wilcox.test(PC2 ~ group, data = plot_data, alternative = "two.sided")
if (result_group_PCoA2$p.value < 2.2e-16){
  x_lab <- "Study\np<2.2e-10"
} else {
  x_lab <- paste0("Group\np=", format(result_group_PCoA2$p.value, scientific = TRUE, digits = 2))
}
g_g_2 <- plot_data %>% 
  ggplot(aes(x=group, y=PC2, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=group_cols, guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  xlab(x_lab) +
  ylab("PCoA2") +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())
p_1 <- plot_grid(g_main, g_s_2, g_g_2, g_s_1, NULL, NULL, g_g_1, NULL, NULL,
               nrow=3, rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))

### adonis in each study
col <- c("#E69F00", "#0072B2", "#D6604D", "#7570B3")
hg <- read.csv("adonis_group.csv")
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

### plot alpha
df_res <- read.csv("alpha_res.csv")
df_plot <- df_res %>%
  mutate(
    alpha_diversity = recode(alpha_diversity, 
                             "shannon" = "Shannon", 
                             "observed_features" = "observed_OTUs",
                             "faith_pd" = "Faith's PD"),
    es_label = sprintf("%.2f", es_g),
    significance = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE           ~ ""  
    ),
    plot_label = paste(es_label, significance, sep = "\n")
  )

p_3 <- ggplot(df_plot, aes(x = alpha_diversity, y = study_id)) +
  geom_tile(aes(fill = es_g), color = "white", linewidth = 1) +
  geom_text(aes(label = plot_label), color = "black", size = 3.5, vjust = 0.5) +
  scale_fill_gradient2(low = "blue", high = "white",
                       midpoint = 1.0, name = "Effect Size") +
  theme_bw() + 
  labs(title = 'Control vs Fiber', y = '', x = '') +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1.2), 
    panel.grid = element_blank(), 
    axis.text = element_text(size = 10, colour = 'black'), 
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.margin = unit(rep(2,4), 'lines'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank() 
  )
plot_bottom <- ggarrange(p_2, p_3, labels = c("b", "c"), common.legend = TRUE, legend="right", align ='hv')

p <- plot_grid(p_1, plot_bottom, align="hv", labels = c("a", ""),
               nrow = 2, ncol=1, plot=FALSE)
ggsave("Figure_1.pdf", p, width = 8, height = 13)

