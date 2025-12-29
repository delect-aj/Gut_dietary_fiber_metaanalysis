library(cowplot)
library(tidyverse)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/difference_analysis")

df_auc <- read.table('aucs_ko.tsv',sep='\t', quote='',stringsAsFactors = FALSE, check.names = FALSE)
p.vals_PRJNA891951 <- read.csv("linda_single_study_ko_PRJNA891951.csv",row.names = 1)
p.vals_PRJEB41443 <- read.csv("linda_single_study_ko_PRJEB41443.csv",row.names = 1)
p.vals_PRJNA293971 <- read.csv("linda_single_study_ko_PRJNA293971.csv",row.names = 1)
p.vals_PRJNA306884 <- read.csv("linda_single_study_ko_PRJNA306884.csv",row.names = 1)
p.vals_PRJNA385004 <- read.csv("linda_single_study_ko_PRJNA385004.csv",row.names = 1)
p.vals_PRJNA428736 <- read.csv("linda_single_study_ko_PRJNA428736.csv",row.names = 1)
p.vals_PRJNA560950 <- read.csv("linda_single_study_ko_PRJNA560950.csv",row.names = 1)
p.vals_PRJNA780023 <- read.csv("linda_single_study_ko_PRJNA780023.csv",row.names = 1)
p.vals_all <- read.csv("combine_results_ko.csv", row.names = 1)
rownames(p.vals_all) <- p.vals_all$diff_otu

p.vals_all1 <- p.vals_all
species.heatmap <- rownames(p.vals_all)[which(p.vals_all$padj < 0.1)]
#str(as.numeric(fc.mat))
#fc.sign <- sign(fc.mat)
#fc.sign[fc.sign == 0] <- 1
p.vals_all$log2FoldChange <- sign(p.vals_all$log2FoldChange)
p.val.signed <- -log10(p.vals_all[species.heatmap,"padj", drop=FALSE])*p.vals_all[species.heatmap, 'log2FoldChange']
species.heatmap.orderd <- rownames(p.val.signed[order(p.val.signed$padj),,drop=FALSE])
p.vals_PRJNA891951 <- p.vals_PRJNA891951[species.heatmap.orderd,]
p.vals_PRJEB41443 <- p.vals_PRJEB41443[species.heatmap.orderd,]
p.vals_PRJNA293971 <- p.vals_PRJNA293971[species.heatmap.orderd,]
p.vals_PRJNA306884 <- p.vals_PRJNA306884[species.heatmap.orderd,]
p.vals_PRJNA385004 <- p.vals_PRJNA385004[species.heatmap.orderd,]
p.vals_PRJNA428736 <- p.vals_PRJNA428736[species.heatmap.orderd,]
p.vals_PRJNA560950 <- p.vals_PRJNA560950[species.heatmap.orderd,]
p.vals_PRJNA780023 <- p.vals_PRJNA780023[species.heatmap.orderd,]
p.vals_all1 <- p.vals_all1[species.heatmap.orderd,]
p.vals.plot <- data.frame(
  PRJNA891951 = p.vals_PRJNA891951$padj,
  PRJEB41443 = p.vals_PRJEB41443$padj,
  PRJNA293971 = p.vals_PRJNA293971$padj,
  PRJNA306884 = p.vals_PRJNA306884$padj,
  PRJNA385004 = p.vals_PRJNA385004$padj,
  PRJNA428736 = p.vals_PRJNA428736$padj,
  PRJNA560950 = p.vals_PRJNA560950$padj,
  PRJNA780023 = p.vals_PRJNA780023$padj,
  all = p.vals_all1$padj
)
rownames(p.vals.plot) <- species.heatmap.orderd
fc.mat.plot <- data.frame(
  PRJNA891951 = p.vals_PRJNA891951$log2FoldChange,
  PRJEB41443 = p.vals_PRJEB41443$log2FoldChange,
  PRJNA293971 = p.vals_PRJNA293971$log2FoldChange,
  PRJNA306884 = p.vals_PRJNA306884$log2FoldChange,
  PRJNA385004 = p.vals_PRJNA385004$log2FoldChange,
  PRJNA428736 = p.vals_PRJNA428736$log2FoldChange,
  PRJNA560950 = p.vals_PRJNA560950$log2FoldChange,
  PRJNA780023 = p.vals_PRJNA780023$log2FoldChange,
  all = p.vals_all1$log2FoldChange
)
rownames(fc.mat.plot) <- species.heatmap.orderd

mx <- max(abs(range(fc.mat.plot, na.rm=TRUE)))
mx <- ifelse(round(mx, digits = 1) < mx, 
             round(mx, digits = 1) + 0.1, 
             round(mx, digits = 1))
brs = seq(-mx, mx, by=0.3)
num.col.steps = length(brs)-1
n = floor(0.45*num.col.steps)
col.hm <- c("#0B559F","#1967AD","#2A7AB9","#3D8DC3","#539ECC","#6BAED6","#88BEDC","#A3CCE3",
            "#BAD6EB","#FFFFFF","#FFFFFF","#FCAF93","#FC9778","#FB8060","#FB6A4A","#F44F38",
            "#E93629","#D52221","#C0151A","#AA1016")
p.vals.plot[is.na(p.vals.plot)] <- 1
fc.mat.plot[is.na(fc.mat.plot)] <- 0
summary(p.vals.plot)
alpha.breaks=c(1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01)
p.vals.bin <- data.frame(apply(p.vals.plot, 2, FUN=.bincode, 
                               breaks = c(0, alpha.breaks, 1), 
                               include.lowest = TRUE),
                         check.names = FALSE)
p.val.greys <-  c(paste0('grey', 
                         round(seq(from=10, to=80, 
                                   length.out = length(alpha.breaks)))), 
                  'white')
names(p.val.greys) <- as.character(1:7)

plot.single.study.heatmap <- function(x){
  df.plot <- tibble(species=factor(rownames(p.vals.plot), 
                                   levels=rev(rownames(p.vals.plot))),
                    p.vals=as.factor(p.vals.bin[[x]]),
                    fc=fc.mat.plot[[x]])
  g1 <- df.plot %>% 
    ggplot(aes(x=species, y=1, fill=fc)) + 
    geom_tile() + theme_bw() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 8), 
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), guide=FALSE) +
    labs(x="", y="")
}
a <- tibble(species=factor(rownames(p.vals.plot), 
                            levels=rev(rownames(p.vals.plot))),
             p.vals=-log10(p.vals.plot$all),
             colour=p.vals > 1)
g1 <- tibble(species=factor(rownames(p.vals.plot), 
                            levels=rev(rownames(p.vals.plot))),
             p.vals=-log10(p.vals.plot$all),
             colour=p.vals > 2) %>% 
  ggplot(aes(x=species, y=p.vals, fill=colour)) + 
  geom_bar(stat='identity') + 
  theme_bw(base_size = 10) + 
  xlab('Meta-analysis significance') + 
  ylab('-log10(padj)') + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10, margin = margin(t = 5)),
        panel.background = element_rect(fill=NULL, colour = 'black')) + 
  scale_y_continuous(limits=c(0, 4), expand = c(0, 0)) + 
  scale_x_discrete(position='top',expand = c(0, 0)) + 
  scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide=FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype="solid", color="red") +  # 添加红线
  annotate("text", x = 11.5, y = 1.7, label = "FDR < 0.05", color = "red", size = 5)
g1
ref.studies <- c("PRJNA891951", "PRJEB41443", "PRJNA293971", "PRJNA306884", "PRJNA428736", "PRJNA560950", "PRJNA780023")
g.lst <- lapply(ref.studies, plot.single.study.heatmap)
p_1 <- plot_grid(g1, g.lst[[1]], g.lst[[2]], g.lst[[3]], g.lst[[4]],g.lst[[5]], g.lst[[6]],g.lst[[7]],
          ncol=1, align = 'v', rel_heights = c(0.3,rep(0.02, 8)))
data_legend1 <- data.frame(
  value = c(1:6),
  p.vals = c(1:6)
)
data_legend1$p.vals <- as.factor(data_legend1$p.vals)
legend1 <- ggplot(data_legend1,aes(x=value, y=1, fill=p.vals)) +
  geom_tile() + theme_bw() + 
  theme(axis.text.x = element_text(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(expand = c(0, 0),breaks = seq(1, 6, length.out = 6),labels = c("1e-05", "1e-04", "1e-03", "1e-02", "1e-01", "0")) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=p.val.greys, na.value='white', guide=FALSE) + labs(title = "Per-study significance\n(adjusted P value)")
legend1
ggsave("legend1_KO.png", legend1, width = 5, height = 1)
data_legend2 <- data.frame(
  value = c(1:57),
  fc = seq(from = -2.8, to = 2.8, by = 0.1)
)
legend2 <- ggplot(data_legend2,aes(x=value, y=1, fill=fc)) +
  geom_tile() + theme_bw() + 
  theme(axis.text.x = element_text(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NULL, colour='black'),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(expand = c(0, 0),breaks = seq(1, 57, by = 7),labels = seq(-2.8, 2.8, by = 0.7)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=col.hm, na.value='white', guide=FALSE) + labs(title = "Generalized fold change\n(log relative abundance)")
legend2
ggsave("legend2_KO.png",legend2,width = 5, height = 1)

## Forest plot
ko_description <- read.delim("/home/dongbiao/software/picrust2-master/picrust2/default_files/description_mapfiles/ko_info.tsv.gz", 
                             sep = "\t", stringsAsFactors = FALSE, header = FALSE) %>%
  column_to_rownames("V1")
colnames(ko_description) <- c("description")
df_auc$ko_id <- sapply(strsplit(df_auc$fid, ":"), function(x) x[2])

marker.set <- rownames(p.val.signed)[
  abs(p.val.signed$padj) > 1.5]
p.val.signed.red <- p.val.signed[marker.set, ,drop=FALSE]
marker.set.orderd <- rev(rownames(p.val.signed.red[order(p.val.signed.red$padj),,
                                                   drop=FALSE]))
df.plot1 <- tibble()
for (i in marker.set.orderd){
  for (s in ref.studies){
    temp <- df_auc[[i]][[s]]
    df.plot1 <- bind_rows(df.plot1, tibble(
      species=i, study=s,
      low=temp[1], auc=temp[2], high=temp[3]
    ))
  }
}
pick_fid <- rownames(p.vals_all)[which(p.vals_all$padj < 0.002)]
df_auc <- df_auc %>% mutate(ko_desc = ko_description[df_auc$ko_id, "description"])
df_auc <- df_auc %>% filter(fid %in% pick_fid)
df_auc <- df_auc %>% drop_na()
df_auc <- df_auc %>% filter(ko_desc != "K09009; uncharacterized protein")

ko_name <- c("ATP-dependent RNA \nhelicase HelY", "putative thioredoxin", "periplasmic iron \nbinding protein",          
                   "two-component system", "gluB; glutamate \ntransport system", "gluC; glutamate \ntransport system", 
                   "gluD; glutamate \ntransport system")
names(ko_name) <- unique(df_auc$ko_desc)
df_auc$ko <- ko_name[df_auc$ko_desc]
g <- df_auc %>% 
  ggplot(aes(x=study, y=AUC)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray", alpha = 0.7) +
  geom_point(pch=23, aes(fill=study), size=3) + 
  facet_grid(~ko, scales = 'free_x', space='free', switch = "x") +
  theme_classic(base_size = 10) + 
  
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(
      ~ ., 
      breaks = c(0.25, 0.75), 
      labels = c("Depleted in Fiber", "Enriched in Fiber") 
    ),
    expand = c(0, 0)
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    
    axis.title.y.right = element_blank(), 
    axis.text.y.right = element_text(angle = 90, hjust = 0.5, size = 10),
    
    panel.grid = element_blank(), 
    strip.placement = "outside",
    panel.spacing.x = unit(0, "lines"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(angle=90, hjust = 1, size = 10)
  ) +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5"),
                    guide = FALSE) + 
  ylab('AUROC') + ggtitle('Effect size of core associations')
  # guides(fill = guide_legend(override.aes = list(size = 5)))

p <- plot_grid(p_1, g, align="hv", labels = c("a", "b"),
               nrow = 2, ncol=1, plot=FALSE)
ggsave(p, filename = "Figure_3.pdf",width = 10, height = 10, useDingbats=FALSE)

