library(qiime2R)
library(vegan)
library(ggExtra)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(cowplot)
meta <- read.delim('metadata.tsv', sep = '\t', stringsAsFactors = FALSE, row.names = 1)
pcoa <- read_qza("betadistance/bray_curtis_pcoa_results.qza")$data
pcoa_df <- pcoa$Vectors %>% column_to_rownames("SampleID") %>% select(PC1, PC2)
plot_data <- cbind(pcoa_df, meta[rownames(pcoa_df), ])
axis_1_title <- paste('PCo1 [', round((pcoa$Eigvals[1]/sum(pcoa$Eigvals))*100,1),'%]', sep='')
axis_2_title <- paste('PCo2 [', round((pcoa$Eigvals[2]/sum(pcoa$Eigvals))*100,1),'%]', sep='')
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
  theme(panel.background = element_rect(fill='white', color = 'black'),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), legend.text = element_text(size=10),
        panel.grid = element_blank(), legend.key = element_blank(), 
        legend.key.size = unit(0.5, "cm"), legend.title = element_blank()
        )
g_main
result_study_PCoA1 <- kruskal.test(PC1 ~ study, data = plot_data)
g_s_1 <- plot_data %>% 
  mutate(study=factor(study, levels = names(study_cols))) %>% 
  ggplot(aes(y = PC1, x = study, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values=study_cols, guide=FALSE) +
  xlab("Study\np < *e-10") +
  theme(axis.ticks = element_blank(),
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(), axis.title.x = element_blank(),
        panel.grid = element_blank()) + 
  coord_flip()
g_s_1
result_study_PCoA2 <- kruskal.test(PC2 ~ study, data = plot_data)
g_s_2 <- plot_data %>% 
  mutate(study = factor(study, levels = names(study_cols))) %>% 
  ggplot(aes(y = PC2, x = study, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values = study_cols, guide = FALSE) +
  scale_x_discrete(position='top') +
  xlab("Study\np < *e-10") +
  theme(axis.ticks=element_blank(), 
        panel.background = element_rect(fill='white', color = 'black'),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
g_s_2
group_cols <- c("0" = "#377EB8", "1" = "#E41A1C")
result_group1_PCoA1 <- wilcox.test(PC1 ~ group, data = plot_data, alternative = "two.sided")
g_g_1 <- plot_data %>% 
  ggplot(aes(x=group, y=PC1, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(name = "Type",values=group_cols, labels = c("Control", "Fiber"), guide=FALSE) + 
  xlab(paste0("Group\np = ", format(result_group1_PCoA1$p.value, scientific = TRUE, digits = 2))) +
  ylab("PCoA1") +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank(),legend.key = element_blank()) + 
  coord_flip()
g_g_1
result_group1_PCoA2 <- wilcox.test(PC2 ~ group, data = plot_data, alternative = "two.sided")
g_g_2 <- plot_data %>% 
  ggplot(aes(x=group, y=PC2, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=group_cols, guide=FALSE) + 
  scale_x_discrete(position='top') + 
  scale_y_continuous(position = 'right') +
  xlab(paste0("Group\np = ", format(result_group1_PCoA2$p.value, scientific = TRUE, digits = 2))) +
  ylab("PCo2") +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='white', color='black'),
        panel.grid = element_blank())
g_g_2
p <- plot_grid(g_main, g_s_2, g_g_2, g_s_1, NULL, NULL, g_g_1, NULL, NULL,
               nrow=3, rel_widths = c(0.8, 0.2, 0.2), rel_heights = c(0.8, 0.2, 0.2))
p
ggsave("Fig1B_bary.pdf",p,width = 12,height = 12)

