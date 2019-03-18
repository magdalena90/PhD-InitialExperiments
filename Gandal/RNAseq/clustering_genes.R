
#setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/Gandal/RNAseq')
setwd('/home/magdalena/PhD/initialExperiments/Gandal/RNAseq')

library(tidyverse); library(reshape2); library(glue); library(plotly); library(RColorBrewer); library(gplots)
library(ConsensusClusterPlus)
library(JADE) ; library(MineICA) ; library(moments) ; library(fdrtool)
library(WGCNA)
library(pdfCluster)
library(ClusterR)

##### Load and transform data ###########################################################

load('./working_data/RNAseq_ASD_4region_DEgenes.Rdata')

reduce_dim_datExpr = function(datExpr, datMeta, var_explained=0.98){
  
  datExpr = data.frame(datExpr)
  
  datExpr_pca = prcomp(datExpr, scale.=TRUE)
  plot(summary(datExpr_pca)$importance[2,], type='b')
  plot(summary(datExpr_pca)$importance[3,], type='b')
  abline(h=var_explained, col='blue')
  
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

datExpr_redDim %>% ggplot(aes(x=PC1, y=PC2)) + geom_point() + theme_minimal()


rm(datSeq, datProbes, reduce_dim_datExpr, reduce_dim_output, datExpr)

##### Clustering ########################################################################

# K-Means Clustering
set.seed(123)
wss = sapply(1:10, function(k) kmeans(datExpr_redDim, k, nstart=25, iter.max=50)$tot.withinss)
plot(wss, type='b', main='K-Means Clustering')
best_k = 4
abline(v=best_k, col='blue')

datExpr_k_means = kmeans(datExpr_redDim, best_k, nstart=25)
km_clusters = datExpr_k_means$cluster


# Hierarchical Clustering
h_clusts = datExpr_redDim %>% dist %>% hclust# %>% as.dendrogram
plot(h_clusts, hang = -1, cex = 0.6, labels=FALSE)
best_k = 5
hc_clusters = cutree(h_clusts, best_k)


# Consensus Clustering
# cc_output = datExpr_redDim %>% as.matrix %>% t %>% ConsensusClusterPlus(maxK=8, reps=5, seed=123)
load('cc_DE/cc_output.RData')
best_k = 5 # 2 clusters and 3 outliers # check mean gene expression for each cluster
cc_clusters = cc_output[[best_k]]$consensusClass
cc_clusters_l1 = cc_clusters

# cc_output_c1 = datExpr_redDim %>% filter(cc_clusters==1) %>% as.matrix %>% t %>% 
#   ConsensusClusterPlus(maxK=8, reps=50, seed=123)
load('cc_DE_c1/cc_output.RData')
best_k = 4
cc_clusters[cc_clusters==1] = cc_output_c1[[best_k]]$consensusClass %>% sapply(function(x) glue('1_', x))

# cc_output_c2 = datExpr_redDim %>% filter(cc_clusters==2) %>% as.matrix %>% t %>% 
#   ConsensusClusterPlus(maxK=8, reps=50, seed=123)
load('cc_DE_c2/cc_output.RData')
best_k = 6
cc_clusters[cc_clusters==2] = cc_output_c2[[best_k]]$consensusClass %>% sapply(function(x) glue('2_', x))


# Independent Component Analysis (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002367)
ICA_output = datExpr_redDim %>% runICA(nbComp=ncol(datExpr_redDim), method='JADE')
signals_w_kurtosis = ICA_output$S %>% data.frame %>% dplyr::select(names(apply(ICA_output$S, 2, kurtosis)>3))
ICA_clusters = apply(signals_w_kurtosis, 2, function(x) fdrtool(x, plot=F)$qval<0.01) %>% data.frame
ICA_clusters = ICA_clusters[colSums(ICA_clusters)>1]

ICA_clusters_min = rep(NA, nrow(ICA_clusters))
ICA_clusters_names = colnames(ICA_clusters) %>% rev
for(c in ICA_clusters_names) ICA_clusters_min[ICA_clusters[,c]] = c

ICA_clusters %>% rowSums %>% table
#    0    1    2    3    4    5    6    7
# 4851  991  390  131   37   13    3    1


# WGCNA
best_power = datExpr_redDim %>% t %>% pickSoftThreshold(powerVector = seq(30, 50, by=1))
# best_power=40 but blockwiseModules only accepts powers up to 30, using 30 instead
network = datExpr_redDim %>% t %>% blockwiseModules(power=30, numericLabels=TRUE)
wgcna_clusters = network$colors
names(wgcna_clusters) = rownames(datExpr_redDim)


# GMM (hard thresholding)
n_clust = datExpr_redDim %>% Optimal_Clusters_GMM(max_clusters=50, criterion='BIC', plot_data=FALSE)
plot(n_clust, type='l')
best_k = n_clust %>% which.min
gmm = datExpr_redDim %>% GMM(best_k)
gmm_clusters = gmm$Log_likelihood %>% apply(1, which.max)

gmm_points = rbind(datExpr_redDim, setNames(data.frame(gmm$centroids), names(datExpr_redDim)))
gmm_labels = c(gmm_clusters, rep(NA, best_k)) %>% as.factor
gmm_points %>% ggplot(aes(x=PC1, y=PC2, color=gmm_labels)) + geom_point() + theme_minimal()


# Manual clustering
manual_clusters = as.factor(as.numeric(0.08*datExpr_redDim$PC1 + 0.2 > datExpr_redDim$PC2))
datExpr_redDim %>% ggplot(aes(x=PC1, y=PC2, color=manual_clusters)) + geom_point() + 
  geom_abline(slope=0.08, intercept=0.2, color='gray') + theme_minimal()
names(manual_clusters) = rownames(datExpr_redDim)

manual_clusters_data = cbind(apply(datExpr_redDim, 1, mean), apply(datExpr_redDim, 1, sd), 
                             manual_clusters) %>% data.frame
colnames(manual_clusters_data) = c('mean','sd','cluster')
manual_clusters_data = manual_clusters_data %>% mutate('cluster'=as.factor(cluster))
manual_clusters_data %>% ggplot(aes(x=mean, color=cluster, fill=cluster)) + 
  geom_density(alpha=0.4) + theme_minimal()
manual_clusters_data %>% ggplot(aes(x=sd, color=cluster, fill=cluster)) + 
  geom_density(alpha=0.4) + theme_minimal()

# # This doesn't work, no matter what separation you use, it always says its significant 
# wilcox.test(x = as.numeric(manual_clusters_data$cluster), y = manual_clusters_data$mean,
#             alternative = 'two.sided')
# wilcox.test(x = as.numeric(manual_clusters_data$cluster), y = manual_clusters_data$sd,
#             alternative = 'two.sided')
#
# wilcox.test(x=c(rep(0,3200), rep(1,3217)), y=manual_clusters_data$sd, alternative='two.sided')

rm(wss, datExpr_k_means, h_clusts, cc_output, cc_output_c1, cc_output_c2, best_k, ICA_output, 
   signals_w_kurtosis, best_power, network, n_clust, gmm, gmm_points, gmm_labels)

##### Plot cluterings ###################################################################

create_2D_plot = function(cat_var, filter_NA=TRUE, plotly=TRUE){
  
  if(filter_NA) plot_points = plot_points %>% dplyr::select(PC1, PC2, cat_var) %>% filter(complete.cases(.))
  
  p = plot_points %>% ggplot(aes_string(x='PC1', y='PC2', color=cat_var)) + 
    geom_point(alpha=0.5) + theme_minimal() + 
    xlab(paste0('PC1 (', round(summary(pca_output)$importance[2,1]*100,2),'%)')) +
    ylab(paste0('PC2 (', round(summary(pca_output)$importance[2,2]*100,2),'%)'))
  
  if(plotly) p = ggplotly(p)
  
  return(p)
}
create_3D_plot = function(cat_var, filter_NA=TRUE){
  
  if(filter_NA) plot_points = plot_points %>% dplyr::select(PC1, PC2, PC3, cat_var) %>% filter(complete.cases(.))
  
  plot_points %>% plot_ly(x=~PC1, y=~PC2, z=~PC3) %>% 
    add_markers(color=plot_points[,cat_var], size=0.5, opacity=0.8) %>% 
    layout(title = glue('Samples coloured by ', cat_var),
           scene = list(xaxis=list(title=glue('PC1 (',round(summary(pca_output)$importance[2,1]*100,2),'%)')),
                        yaxis=list(title=glue('PC2 (',round(summary(pca_output)$importance[2,2]*100,2),'%)')),
                        zaxis=list(title=glue('PC3 (',round(summary(pca_output)$importance[2,3]*100,2),'%)'))))  
}

# ICA
ICA_clusters %>% mutate(cl_sum=rowSums(.)) %>% as.matrix %>% melt %>% ggplot(aes(Var2, Var1)) + 
  geom_tile(aes(fill=value)) + xlab('Clusters') + ylab('Samples') + 
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + coord_flip()

# All clusterings
plot_points = datExpr_redDim %>% data.frame() %>% dplyr::select(PC1:PC3) %>%
  mutate(ID = rownames(.), cluster_km = as.factor(km_clusters), 
         cluster_cc_l1 = as.factor(cc_clusters_l1),
         cluster_cc = as.factor(cc_clusters), 
         cluster_hc = as.factor(hc_clusters),
         cluster_ica = as.factor(ICA_clusters_min),
         cluster_wgcna = as.factor(wgcna_clusters),
         cluster_gmm = as.factor(gmm_clusters),
         cluster_manual = as.factor(manual_clusters),
         n_cluster_ICA = as.factor(rowSums(ICA_clusters)))

create_2D_plot('cluster_km')
create_2D_plot('cluster_hc')
create_2D_plot('cluster_cc_l1')
create_2D_plot('cluster_cc')
create_2D_plot('cluster_ica')
create_2D_plot('cluster_wgcna')
create_2D_plot('cluster_gmm')

create_3D_plot('cluster_km') ; create_3D_plot('cluster_cc')
create_3D_plot('cluster_hc') ; create_3D_plot('cluster_ica') 
create_3D_plot('cluster_wgcna')


rm(create_2D_plot, create_3D_plot, ICA_clusters_names, c)

##### Compare Clusterings ###############################################################

# Adjusted Rand Index
clusterings = list(km_clusters, hc_clusters, cc_clusters_l1, cc_clusters, ICA_clusters_min, 
                   wgcna_clusters, gmm_clusters, manual_clusters)
cluster_sim = data.frame(matrix(nrow = length(clusterings), ncol = length(clusterings)))
for(i in 1:(length(clusterings))){
  cluster1 = clusterings[[i]]
  for(j in (i):length(clusterings)){
    cluster2 = clusterings[[j]]
    cluster_sim[i,j] = adj.rand.index(cluster1, cluster2)
  }
}
colnames(cluster_sim) = c('K-means','Hierarchical','Consensus L1', 'Consensus', 'ICA', 'WGCNA', 'GMM', 'Manual')
rownames(cluster_sim) = colnames(cluster_sim)

cluster_sim = cluster_sim %>% as.matrix %>% round(2)
heatmap.2(x = cluster_sim, Rowv = FALSE, Colv = FALSE, dendrogram = 'none', 
          cellnote = cluster_sim, notecol = 'black', trace = 'none', key = FALSE, 
          cexRow=1, cexCol=1, margins=c(7,7))


rm(i, j, cluster1, cluster2, clusterings, cluster_sim)

##### Enrichment Analysis ###############################################################

# Create topGO object
mod_n_genes = function(allModules) { return(allModules == 1) }
GOdata = new('topGOdata', ontology = 'BP', allGenes = gene_list, geneSel = mod_n_genes,
             nodeSize = 10, annot = annFUN.db, affyLib = 'illuminaHumanv4.db')

# Perform statistical tests
fisher = runTest(GOdata, statistic = 'fisher')
all_res = GenTable(GOdata, fisher=fisher, orderBy='fisher', ranksOf='fisher',
                   topNodes=10)

# Plot significant genes in the GO tree
print(showSigOfNodes(GOdata, score(fisher), firstSigNodes=5, useInfo='all'))

View(datMeta[gsub('X', '', names(sort(abs(pca_output$rotation[,1]^2), decreasing=TRUE))[1:10]),])
