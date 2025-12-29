library(qiime2R)
library(tidyverse)

setwd("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/network/")

abundance_table <- read_q2biom("../projects/table_6721_2378.biom") %>% as.data.frame()
metadata <- read.csv("../diversity/metadata.tsv", sep = "\t", row.names = 1)
metadata <- metadata %>% filter(!(study %in% c("PRJEB2165")))

relative_abundance <- function(df) {
  sweep(df, 2, colSums(df), "/") * 100
}

study <- unique(metadata$study)
for (i in study){
  for (j in c("fiber", "control")){
    if (j == "fiber") {
      metadata_pick <- metadata %>% filter(study == i) %>% filter(group == 1)
    } else {
      metadata_pick <- metadata %>% filter(study == i) %>% filter(group == 0)
    }
    abundance_table_pick <- abundance_table[, rownames(metadata_pick)]
    abundance_table_pick <- abundance_table_pick[rowSums(abundance_table_pick) != 0, ]
    
    abundance_table_pick <- relative_abundance(abundance_table_pick)
    min_samples <- ceiling(0.1 * ncol(abundance_table_pick))  
    threshold <- 0.1 
    filter_criteria <- rowSums(abundance_table_pick >= threshold) >= min_samples
    abundance_table_pick <- abundance_table_pick[filter_criteria, ] %>%
      rownames_to_column(var = "X")
    if (j == "fiber") {
      write.table(abundance_table_pick, paste0("fiber_dir/", i, ".txt"), 
                  sep="\t", row.names = FALSE)
    } else {
      write.table(abundance_table_pick, paste0("control_dir/", i, ".txt"), 
                  sep="\t", row.names = FALSE)
    }
  }
}

