library(dplyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(cowplot)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/classification/")

### OTUS, Within-study cross validation，Cross-study validation ，Leave-one-study-out validation
df_cv <- read.csv("results_cv.csv")
df_cv$`Within-study cross validation (CV)` <- c("Within-study cross validation (CV)")
df_csv <- read.csv("results_csv.csv")
df_csv$`Cross-study validation (CSV)` <- c("Cross-study validation (CSV)")
df_loso <- read.csv("results_loso.csv")
df_loso$`Leave-one-study-out validation (LOSO)` <- c("Leave-one-study-out validation (LOSO)")
plot_results <- data.frame(test_id = c(df_cv$study, df_csv$test_study, df_loso$test_study),
                           group = c(df_cv$`Within-study cross validation (CV)`,
                                     df_csv$`Cross-study validation (CSV)`,
                                     df_loso$`Leave-one-study-out validation (LOSO)`),
                           auc = c(df_cv$auc, df_csv$auc, df_loso$auc))
plot_results <- plot_results %>%
  mutate(group = factor(group, 
                        levels = c("Within-study cross validation (CV)",
                                   "Cross-study validation (CSV)",
                                   "Leave-one-study-out validation (LOSO)"),
                        labels = c("CV", "CSV", "LOSO")))
plot_results$test_id <- as.factor(plot_results$test_id)

color_values <- c("CV" = "#0072B2", "CSV" = "grey30", "LOSO" = "#E69F00")
shape_values <- c("CV" = 19, "CSV" = 19, "LOSO" = 17)
legend_labels <- c("CV" = "Within-study cross validation (CV)",
                   "CSV" = "Cross-study validation (CSV)",
                   "LOSO" = "Leave-one-study-out validation (LOSO)")

overall_stats <- plot_results %>%
  group_by(group) %>%
  summarize(overall_mean = mean(auc),
            overall_sd = sd(auc))
cv_stats <- filter(overall_stats, group == "CV")
loso_stats <- filter(overall_stats, group == "LOSO")

p_1 <- ggplot(plot_results, aes(x = test_id, y = auc)) +
  geom_hline(yintercept = cv_stats$overall_mean, color = "#0072B2", linewidth = 1) +
  geom_hline(yintercept = loso_stats$overall_mean, color = "#E69F00", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_violin(
    data = . %>% filter(group == "CSV"),
    aes(fill = "CSV"), 
    color = NA, 
    alpha = 0.5,
    trim = TRUE
  ) +
  stat_summary(
    aes(color = group, shape = group), 
    fun = mean, 
    geom = "point", 
    size = 3,
    stroke = 1.2 
  ) +
  scale_fill_manual(values = c("CSV" = "grey80"), guide = "none") + 
  scale_color_manual(
    name = NULL, 
    values = color_values,
    labels = legend_labels
  ) +
  scale_shape_manual(
    name = NULL, 
    values = shape_values,
    labels = legend_labels
  ) +
  scale_y_continuous(limits = c(0.45, 0.85)) +
  labs(x = NULL, y = "AUC") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 10)
  )

### subset of OTUs
df_loso <- read.csv("results_loso.csv")
df_loso$group <- "Full"
df_biomarker <- read.csv("results_loso_biomark.csv")
df_biomarker$group <- "Biomarker"
plot_results <- rbind(df_loso, df_biomarker)


# df_plot$test_study <- factor(df_plot$test_study, levels = sort(unique(df_plot$test_study)))

avg_aucs <- plot_results %>%
  group_by(group) %>%
  summarize(avg_auc_val = mean(auc))

legend_labels <- c(
  "Full" = paste0("Full\n[avg. AUC: ", sprintf("%.3f", avg_aucs$avg_auc_val[avg_aucs$group == "Full"]), "]"),
  "Biomarker" = paste0("Biomarker\n[avg. AUC: ", sprintf("%.3f", avg_aucs$avg_auc_val[avg_aucs$group == "Biomarker"]), "]")
)

p_2 <- ggplot(plot_results, aes(x = test_study, y = auc)) +
  geom_line(
    aes(group = test_study), 
    color = "grey60", 
    linewidth = 1.5
  ) +
  geom_point(
    aes(shape = group, fill = group),
    size = 4,    
    color = "black"  
  ) +
  labs(y = "AUC", x = NULL) +
  scale_shape_manual(
    name = "Feature set", 
    values = c("Full" = 21, "Biomarker" = 24),
    labels = legend_labels 
  ) +
  scale_fill_manual(
    name = "Feature set", 
    values = c("Full" = "yellow", "Biomarker" = "yellow"),
    labels = legend_labels
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 10)
  )

### KOs, Within-study cross validation，Cross-study validation ，Leave-one-study-out validation
df_cv <- read.csv("results_ko_cv.csv")
df_cv$`Within-study cross validation (CV)` <- c("Within-study cross validation (CV)")
df_csv <- read.csv("results_ko_csv.csv")
df_csv$`Cross-study validation (CSV)` <- c("Cross-study validation (CSV)")
df_loso <- read.csv("results_ko_loso.csv")
df_loso$`Leave-one-study-out validation (LOSO)` <- c("Leave-one-study-out validation (LOSO)")
plot_results <- data.frame(test_id = c(df_cv$study, df_csv$test_study, df_loso$test_study),
                           group = c(df_cv$`Within-study cross validation (CV)`,
                                     df_csv$`Cross-study validation (CSV)`,
                                     df_loso$`Leave-one-study-out validation (LOSO)`),
                           auc = c(df_cv$auc, df_csv$auc, df_loso$auc))
plot_results <- plot_results %>%
  mutate(group = factor(group, 
                        levels = c("Within-study cross validation (CV)",
                                   "Cross-study validation (CSV)",
                                   "Leave-one-study-out validation (LOSO)"),
                        labels = c("CV", "CSV", "LOSO")))
plot_results$test_id <- as.factor(plot_results$test_id)

color_values <- c("CV" = "#0072B2", "CSV" = "grey30", "LOSO" = "#E69F00")
shape_values <- c("CV" = 19, "CSV" = 19, "LOSO" = 17)
legend_labels <- c("CV" = "Within-study cross validation (CV)",
                   "CSV" = "Cross-study validation (CSV)",
                   "LOSO" = "Leave-one-study-out validation (LOSO)")

overall_stats <- plot_results %>%
  group_by(group) %>%
  summarize(overall_mean = mean(auc),
            overall_sd = sd(auc))
cv_stats <- filter(overall_stats, group == "CV")
loso_stats <- filter(overall_stats, group == "LOSO")

p_3 <- ggplot(plot_results, aes(x = test_id, y = auc)) +
  geom_hline(yintercept = cv_stats$overall_mean, color = "#0072B2", linewidth = 1) +
  geom_hline(yintercept = loso_stats$overall_mean, color = "#E69F00", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_violin(
    data = . %>% filter(group == "CSV"),
    aes(fill = "CSV"), 
    color = NA, 
    alpha = 0.5,
    trim = TRUE
  ) +
  stat_summary(
    aes(color = group, shape = group), 
    fun = mean, 
    geom = "point", 
    size = 3,
    stroke = 1.2 
  ) +
  scale_fill_manual(values = c("CSV" = "grey80"), guide = "none") + 
  scale_color_manual(
    name = NULL, 
    values = color_values,
    labels = legend_labels
  ) +
  scale_shape_manual(
    name = NULL, 
    values = shape_values,
    labels = legend_labels
  ) +
  labs(x = NULL, y = "AUC") +
  scale_y_continuous(limits = c(0.45, 0.85)) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 10)
  )

### subsets KOs
df_loso <- read.csv("results_ko_loso.csv")
df_loso$group <- "Full"
df_biomarker <- read.csv("results_ko_loso_biomark.csv")
df_biomarker$group <- "Biomarker"
plot_results <- rbind(df_loso, df_biomarker)


# df_plot$test_study <- factor(df_plot$test_study, levels = sort(unique(df_plot$test_study)))

avg_aucs <- plot_results %>%
  group_by(group) %>%
  summarize(avg_auc_val = mean(auc))

legend_labels <- c(
  "Full" = paste0("Full\n[avg. AUC: ", sprintf("%.3f", avg_aucs$avg_auc_val[avg_aucs$group == "Full"]), "]"),
  "Biomarker" = paste0("Biomarker\n[avg. AUC: ", sprintf("%.3f", avg_aucs$avg_auc_val[avg_aucs$group == "Biomarker"]), "]")
)

p_4 <- ggplot(plot_results, aes(x = test_study, y = auc)) +
  geom_line(
    aes(group = test_study), 
    color = "grey60", 
    linewidth = 1.5
  ) +
  geom_point(
    aes(shape = group, fill = group),
    size = 4,    
    color = "black"  
  ) +
  labs(y = "AUC", x = NULL) +
  scale_shape_manual(
    name = "Feature set", 
    values = c("Full" = 21, "Biomarker" = 24),
    labels = legend_labels 
  ) +
  scale_fill_manual(
    name = "Feature set", 
    values = c("Full" = "yellow", "Biomarker" = "yellow"),
    labels = legend_labels
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom", 
    legend.text = element_text(size = 10)
  )

p <- plot_grid(p_1, p_3, p_2, p_4, align="hv", labels = c("a", "c", "b", "d"),
               nrow = 2, ncol=2, plot=FALSE)
ggsave(p, filename = "Figure_4.pdf",width = 10, height = 7, useDingbats=FALSE)
