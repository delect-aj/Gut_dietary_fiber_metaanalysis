### load package
library(NetMoss2)
library(rsparcc)
setwd('/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/network')

### read directory
fiber_dir = paste0(getwd(),"/fiber_dir")
control_dir = paste0(getwd(),"/control_dir")

# construct networks
netBuild(case_dir = fiber_dir,
        control_dir = control_dir,
        method = "sparcc")

net_case_dir = paste0(getwd(),"/net_case_dir")
net_control_dir = paste0(getwd(),"/net_control_dir")
#calculate NetMoss score
nodes_result = NetMoss(case_dir = fiber_dir,    
                       control_dir = control_dir,    
                       net_case_dir = net_case_dir,   
                       net_control_dir = net_control_dir) 
# saveRDS(nodes_result, file = "nodes_result.rds")
nodes_result <- readRDS("nodes_result.rds")
result <- nodes_result[[1]]
netPlot(result = nodes_result,
        num.top = 5,
        num.score = 30,
        e.th = 0.1,
        my.layout = layout_components,
        my.label = TRUE)


# write.csv(result, "netmoss.csv", row.names = FALSE)
# plot networks

### Fiber 
importance_fid <- read.csv("../difference_analysis/combine_results.csv", row.names = 1) %>%
  filter(padj < 0.05)
importance_fid <- importance_fid$diff_otu
fiber_cor <- nodes_result[[2]]
fiber_cor <- fiber_cor %>%
  as.data.frame() %>%
  set_names(rownames(fiber_cor)) %>%
  rownames_to_column("id_1") %>%
  pivot_longer(cols = -id_1, names_to = "id_2",
    values_to = "cor") %>%
  filter(match(id_1, rownames(fiber_cor)) < match(id_2, rownames(fiber_cor)) &
      cor != 0)
fiber_cor <- fiber_cor %>% filter(abs(cor) >= 0.05)
temp <- fiber_cor %>% filter(id_1 %in% importance_fid)
keep_fid_1 <- unique(temp$id_2)
temp <- fiber_cor %>% filter(id_2 %in% importance_fid)
keep_fid_2 <- unique(temp$id_1)

keep_fid <- intersect(result$taxon_names, c(importance_fid, keep_fid_1, keep_fid_2))
fiber_cor <- nodes_result[[2]]
fiber_cor <- fiber_cor[keep_fid, keep_fid]
fiber_cor <- fiber_cor %>%
  as.data.frame() %>%
  set_names(rownames(fiber_cor)) %>%
  rownames_to_column("id_1") %>%
  pivot_longer(cols = -id_1, names_to = "id_2",
               values_to = "cor") %>%
  filter(match(id_1, rownames(fiber_cor)) < match(id_2, rownames(fiber_cor)) &
           cor != 0)
fiber_cor <- fiber_cor %>% filter(abs(cor) >= 0.05)
node_name <- data.frame(fid = unique(c(fiber_cor$id_1, fiber_cor$id_2)))
fiber_cor <- fiber_cor %>% mutate(label = ifelse(cor > 0 , "positive", "negtive"))
fiber_cor$cor <- abs(fiber_cor$cor)
write.csv(fiber_cor, "fiber_cor.csv", row.names = FALSE, quote = FALSE)

### node relative abundance
tax <- read.delim("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv", 
                  sep = "\t", stringsAsFactors = FALSE, row.names = 1)
tax_df <- tax %>% .[node_name$fid, , drop = FALSE] %>%
  mutate(Taxon_split = str_split(Taxon, ";")) %>%
  mutate(Taxon_split = map(Taxon_split, ~ {
    if (length(.x) < 7) c(.x, rep(NA, 7 - length(.x))) else .x[1:7]
  })) %>% separate_wider_delim(Taxon, delim = ";", 
                               names = c('Kingdom', 'Phylum', 'Class', 'Order', 
                                         'Family', 'Genus', 'Species'),
                               too_few = "align_start",
                               too_many = "drop") %>% select(-Taxon_split)

rownames(tax_df) <- node_name$fid
node_name$species <- tax_df$Species
node_name$phylum <- tax_df$Phylum
node_name$netmoss_score <- result[node_name$fid, "NetMoss_Score"]
write.csv(node_name, "fiber_node.csv", row.names = FALSE, quote = FALSE)


### control
importance_fid <- read.csv("../difference_analysis/combine_results.csv", row.names = 1) %>%
  filter(padj < 0.05)
importance_fid <- importance_fid$diff_otu
control_cor <- nodes_result[[3]]
control_cor <- control_cor[keep_fid, keep_fid]
control_cor <- control_cor %>%
  as.data.frame() %>%
  set_names(rownames(control_cor)) %>%
  rownames_to_column("id_1") %>%
  pivot_longer(cols = -id_1, names_to = "id_2",
               values_to = "cor") %>%
  filter(match(id_1, rownames(control_cor)) < match(id_2, rownames(control_cor)) &
           cor != 0)
control_cor <- control_cor %>% filter(abs(cor) >= 0.05)
node_name <- data.frame(fid = unique(c(control_cor$id_1, control_cor$id_2)))
control_cor <- control_cor %>% mutate(label = ifelse(cor > 0 , "positive", "negtive"))
control_cor$cor <- abs(control_cor$cor)
write.csv(control_cor, "control_cor.csv", row.names = FALSE, quote = FALSE)

### node relative abundance
tax <- read.delim("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv", 
                  sep = "\t", stringsAsFactors = FALSE, row.names = 1)
tax_df <- tax %>% .[node_name$fid, , drop = FALSE] %>%
  mutate(Taxon_split = str_split(Taxon, ";")) %>%
  mutate(Taxon_split = map(Taxon_split, ~ {
    if (length(.x) < 7) c(.x, rep(NA, 7 - length(.x))) else .x[1:7]
  })) %>% separate_wider_delim(Taxon, delim = ";", 
                               names = c('Kingdom', 'Phylum', 'Class', 'Order', 
                                         'Family', 'Genus', 'Species'),
                               too_few = "align_start",
                               too_many = "drop") %>% select(-Taxon_split)
rownames(tax_df) <- node_name$fid
node_name$species <- tax_df$Species
node_name$phylum <- tax_df$Phylum
node_name$netmoss_score <- result[node_name$fid, "NetMoss_Score"]
write.csv(node_name, "control_node.csv", row.names = FALSE, quote = FALSE)
