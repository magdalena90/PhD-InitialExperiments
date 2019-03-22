
setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/Gandal/RNAseq')

library(tidyverse) ; library(reshape2) ; library(glue) ; library(plotly)
library(ConsensusClusterPlus)

##### Load and transform data ###########################################################

#load('./working_data/RNAseq_ASD_4region_normalized.Rdata')
load('./working_data/RNAseq_ASD_4region_DEgenes_adj_pval_lfc.Rdata')

reduce_dim_datExpr = function(datExpr, datMeta, var_explained=0.98){
  
  datExpr = data.frame(datExpr)
  datExpr_pca = prcomp(datExpr, scale=TRUE)
  last_pc = data.frame(summary(datExpr_pca)$importance[3,]) %>% rownames_to_column(var='ID') %>% 
    filter(.[[2]] >= var_explained) %>% top_n(-1, ID)
  
  print(glue('Keeping top ', substr(last_pc$ID, 3, nchar(last_pc$ID)), ' components that explain ',
             var_explained*100, '% of the variance'))
  
  datExpr_top_pc = datExpr_pca$x %>% data.frame %>% dplyr::select(PC1:last_pc$ID)
  
  return(list('datExpr'=datExpr_top_pc, 'pca_output'=datExpr_pca))
}

reduce_dim_output = reduce_dim_datExpr(datExpr, datMeta)
datExpr_redDim = reduce_dim_output$datExpr
pca_output = reduce_dim_output$pca_output


rm(datSeq, datProbes, reduce_dim_datExpr, reduce_dim_output, datExpr)

##### Clustering ########################################################################

# # Consensus Clustering
# cc_gene_list = sample(rownames(datExpr_redDim), 10000)
# cc_output = datExpr_redDim[cc_gene_list,] %>% as.matrix %>% t %>%
#   ConsensusClusterPlus(maxK=10, reps=5, seed=123, title='R_markdowns/clustering_genes_03_08/cc_l1', plot='png')
# 
# save(cc_output, cc_gene_list, file='R_markdowns/clustering_genes_03_08/cc_l1/cc_output.RData')
# 
# load('R_markdowns/clustering_genes_03_08/cc_l1/cc_output.RData')
# cc_clusters = cc_output[[4]]$consensusClass
# 
# cc_output_c1 = datExpr_redDim[cc_gene_list,] %>% filter(cc_clusters==1) %>% as.matrix %>% t %>%
#    ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='R_markdowns/clustering_genes_03_08/cc_l2_1/', plot='png')
# save(cc_output_c1, file='R_markdowns/clustering_genes_03_08/cc_l2_1/cc_output.RData')
# 
# cc_output_c2 = datExpr_redDim[cc_gene_list,] %>% filter(cc_clusters==2) %>% as.matrix %>% t %>%
#   ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='R_markdowns/clustering_genes_03_08/cc_l2_2/', plot='png')
# save(cc_output_c2, file='R_markdowns/clustering_genes_03_08/cc_l2_2/cc_output.RData')

# Removing first PC
# new_colnames = colnames(datExpr_redDim)[-ncol(datExpr_redDim)]
# datExpr_redDim = datExpr_redDim %>% dplyr::select(-PC1)
# colnames(datExpr_redDim) = new_colnames
# 
# # Consensus Clustering
# cc_output = datExpr_redDim %>% as.matrix %>% t %>%
#   ConsensusClusterPlus(maxK=10, reps=5, seed=123, title='R_markdowns/clustering_genes_03_21/cc_l1', plot='png')
# 
# save(cc_output, file='R_markdowns/clustering_genes_03_21/cc_l1/cc_output.RData')


load('R_markdowns/clustering_genes_03_21/cc_l1/cc_output.RData')
cc_clusters = cc_output[[6]]$consensusClass

cc_output_c1 = datExpr_redDim %>% filter(cc_clusters==1) %>% as.matrix %>% t %>%
   ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='R_markdowns/clustering_genes_03_21/cc_l2_1/', plot='png')
save(cc_output_c1, file='R_markdowns/clustering_genes_03_21/cc_l2_1/cc_output.RData')

cc_output_c2 = datExpr_redDim %>% filter(cc_clusters==2) %>% as.matrix %>% t %>%
  ConsensusClusterPlus(maxK=8, reps=50, seed=123, title='R_markdowns/clustering_genes_03_21/cc_l2_2/', plot='png')
save(cc_output_c2, file='R_markdowns/clustering_genes_03_21/cc_l2_2/cc_output.RData')
