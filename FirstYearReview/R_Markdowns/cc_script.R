
library(tidyverse)
library(ConsensusClusterPlus)

###############################################################################################################################
###############################################################################################################################
# ALL GENES

# Gandal dataset
load('./../Data/Gandal/preprocessed_data.RData')
datExpr = datExpr %>% data.frame
DE_info = DE_info %>% data.frame

# SFARI Genes
SFARI_genes = read_csv('./../Data/SFARI/SFARI_genes_with_ensembl_IDs.csv')

# GO Annotations
GO_annotations = read.csv('./../Data/GO_annotations/genes_GO_annotations.csv')
GO_neuronal = GO_annotations %>% filter(grepl('neuron', go_term)) %>% 
  mutate('ID' = as.character(ensembl_gene_id)) %>% 
  dplyr::select(-ensembl_gene_id) %>% distinct(ID) %>%
  mutate('Neuronal' = 1)

# Add SFARI scores and Neuronal functionality to DE_info
DE_info = DE_info %>% mutate('ID'=rownames(.)) %>% left_join(SFARI_genes, by='ID') %>% 
  mutate(`gene-score`=ifelse(is.na(`gene-score`), 'None', `gene-score`)) %>%
  distinct(ID, .keep_all = TRUE) %>% left_join(GO_neuronal, by='ID') %>%
  mutate(Neuronal=ifelse(is.na(Neuronal), 0, Neuronal)) %>%
  mutate(gene.score=ifelse(`gene-score`=='None' & Neuronal==1, 'Neuronal', `gene-score`))

rm(GO_annotations)

plot_data = data.frame('id'=rownames(datExpr), 'mean_expression' = rowMeans(datExpr), 'SFARI_score'=DE_info$`gene-score`!='None')

datExpr = datExpr[plot_data$mean_expression>6,]
DE_info = DE_info[plot_data$mean_expression>6,]
datProbes = datProbes[plot_data$mean_expression>6,]

pca = prcomp(datExpr)
datExpr_redDim = pca$x %>% data.frame %>% dplyr::select(PC1:PC32)


###############################################################################################################################
# CONSENSUS CLUSTERING

# Level 1
cc_output = datExpr_redDim %>% as.matrix %>% t %>% 
            ConsensusClusterPlus(maxK=10, reps=5, seed=123, title='./../Data/Gandal/consensusClustering/probes/l1', plot='png')
save(cc_output, file='./../Data/Gandal/consensusClustering/probes/l1/cc_output.RData')

# # Level 2
# best_k = 4
# cc_output_c1 = datExpr_redDim %>% filter(cc_output[[best_k]]$consensusClass==1) %>% as.matrix %>% t %>%
#   ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='./../Data/Gandal/consensusClustering/probes/l2_1', plot='png')
# save(cc_output_c1, file='./../Data/Gandal/consensusClustering/probes/l2_1/cc_output.RData')
# 
# cc_output_c1 = datExpr_redDim %>% filter(cc_output[[best_k]]$consensusClass==2) %>% as.matrix %>% t %>%
#   ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='./../Data/Gandal/consensusClustering/probes/l2_2', plot='png')
# save(cc_output_c1, file='./../Data/Gandal/consensusClustering/probes/l2_2/cc_output.RData')
