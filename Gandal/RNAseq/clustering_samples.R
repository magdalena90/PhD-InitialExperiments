
#setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/Gandal/RNAseq')
setwd('/home/magdalena/PhD/initialExperiments/Gandal/RNAseq')

library(tidyverse) ; library(reshape2) ; library(glue) ; library(plotly) ; library(dendextend)
library(viridis)
library(ConsensusClusterPlus)
library(JADE) ; library(MineICA) ; library(moments) ; library(fdrtool)
library(ClusterR)
library(WGCNA)
library(pdfCluster) ; library(gplots)

##### Load and transform data ###########################################################

load('./working_data/RNAseq_ASD_4region_DEgenes.Rdata')

reduce_dim_datExpr = function(datExpr, datMeta, var_explained=0.8, filter_controls=FALSE){

  datExpr = data.frame(datExpr)
  
  if(filter_controls){
    datMeta = datMeta %>% filter(Diagnosis_=='ASD')
    datExpr = datExpr %>% select(paste0('X', datMeta_ASD$Dissected_Sample_ID))
  }
  
  datExpr_pca = prcomp(t(datExpr), scale=TRUE)
  plot(summary(datExpr_pca)$importance[2,], type='b')
  plot(summary(datExpr_pca)$importance[3,], type='b')
  abline(h=var_explained, col='blue')
  
  last_pc = data.frame(summary(datExpr_pca)$importance[3,]) %>% rownames_to_column(var='ID') %>% 
    filter(.[[2]] >= var_explained) %>% top_n(-1, ID)
  
  print(glue('Keeping top ', substr(last_pc$ID, 3, nchar(last_pc$ID)), ' components that explain ',
             var_explained*100, '% of the variance'))
  
  datExpr_top_pc = datExpr_pca$x %>% data.frame %>% dplyr::select(PC1:last_pc$ID)
  
  return(list('datExpr'=datExpr_top_pc, 'datMeta'=datMeta, 'pca_output'=datExpr_pca))
}

reduce_dim_output = reduce_dim_datExpr(datExpr, datMeta)
datExpr_redDim = reduce_dim_output$datExpr
datMeta_redDim = reduce_dim_output$datMeta
pca_output = reduce_dim_output$pca_output


rm(datSeq, datProbes, reduce_dim_datExpr, reduce_dim_output, datExpr, datMeta)

##### Clustering ########################################################################

# K-Means Clustering
set.seed(123)
wss = sapply(1:15, function(k) kmeans(datExpr_redDim, k, nstart=25)$tot.withinss)
plot(wss, type='b', main='K-Means Clustering')
best_k = 5

datExpr_k_means = kmeans(datExpr_redDim, best_k, nstart=25)
km_clusters = datExpr_k_means$cluster

# Hierarchical Clustering
  h_clusts = datExpr_redDim %>% dist %>% hclust %>% as.dendrogram
  h_clusts %>% plot
  best_k = 6
  hc_clusters = cutree(h_clusts, best_k)
  
  create_viridis_dict = function(age){
    min_age = datMeta_redDim$Age %>% min
    max_age = datMeta_redDim$Age %>% max
    viridis_age_cols = viridis(max_age - min_age + 1)
    names(viridis_age_cols) = seq(min_age, max_age)
    
    return(viridis_age_cols)
  }
  viridis_age_cols = create_viridis_dict()
  
  dend_meta = datMeta_redDim[match(gsub('X','',labels(h_clusts)), rownames(datMeta_redDim)),] %>% 
              mutate('Diagnosis' = ifelse(Diagnosis_=='CTL','#008080','#86b300'), # Blue control, Green ASD
                     'Sex' = ifelse(Sex=='F','#ff6666','#008ae6'),                # Pink Female, Blue Male
                     'Region' = case_when(Brain_lobe=='Frontal'~'#F8766D',        # ggplot defaults for 4 colours
                                        Brain_lobe=='Temporal'~'#7CAE00',
                                        Brain_lobe=='Parietal'~'#00BFC4',
                                        Brain_lobe=='Occipital'~'#C77CFF'),
                     'Age' = viridis_age_cols[as.character(Age)]) %>%            # Purple: young, Yellow: old
              dplyr::select(Age, Region, Sex, Diagnosis)
  h_clusts %>% set('labels', rep('', nrow(datMeta_redDim))) %>% set('branches_k_color', k=best_k) %>% plot
  colored_bars(colors=dend_meta)

# Consensus Clustering
cc_output = datExpr_redDim %>% as.matrix %>% t %>% ConsensusClusterPlus(maxK=5, reps=50, seed=123) # 2
best_k = 2
cc_clusters = cc_output[[best_k]]$consensusClass

cc_output_c1 = datExpr_redDim %>% filter(cc_clusters==1) %>% as.matrix %>% t %>% 
                                  ConsensusClusterPlus(maxK=8, reps=50, seed=123)
best_k = 4
cc_clusters[cc_clusters==1] = cc_output_c1[[best_k]]$consensusClass %>% sapply(function(x) glue('1_', x))

cc_output_c2 = datExpr_redDim %>% filter(cc_clusters==2) %>% as.matrix %>% t %>% 
                                  ConsensusClusterPlus(maxK=8, reps=50, seed=123)
best_k = 6
cc_clusters[cc_clusters==2] = cc_output_c2[[best_k]]$consensusClass %>% sapply(function(x) glue('2_', x))

# Independent Component Analysis (www.journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002367)
# Run ICA with that same number of nbComp to then filter them:
ICA_output = datExpr_redDim %>% runICA(nbComp=ncol(datExpr_redDim), method='JADE')
# Select components with kurtosis > 3
signals_w_kurtosis = ICA_output$S %>% data.frame %>% dplyr::select(names(apply(ICA_output$S, 2, kurtosis)>3))
# Assign obs to genes with FDR<0.001 using the fdrtool package
ICA_clusters = apply(signals_w_kurtosis, 2, function(x) fdrtool(x, plot=F)$qval<0.001) %>% data.frame
# Remove clusters with zero or one elements
ICA_clusters = ICA_clusters[colSums(ICA_clusters)>1]

ICA_clusters_min = rep(NA, nrow(ICA_clusters))
ICA_clusters_names = colnames(ICA_clusters)
for(c in ICA_clusters_names) ICA_clusters_min[ICA_clusters[,c]] = c


ICA_clusters %>% rowSums %>% table


# WGCNA
best_power = datExpr_redDim %>% t %>% pickSoftThreshold(powerVector = seq(10, 30, by=1)) # weird jumps
network = datExpr_redDim %>% t %>% blockwiseModules(power=best_power$powerEstimate, numericLabels=TRUE)
wgcna_clusters = network$colors
names(wgcna_clusters) = rownames(datExpr_redDim)

# GMM (hard thresholding)
n_clust = datExpr_redDim %>% Optimal_Clusters_GMM(max_clusters=80, criterion='BIC', plot_data=FALSE)
plot(n_clust, type='l') # ?!
best_k = 5
gmm = datExpr_redDim %>% GMM(best_k)
gmm_clusters = gmm$Log_likelihood %>% apply(1, which.max)

gmm_points = rbind(datExpr_redDim, setNames(data.frame(gmm$centroids), names(datExpr_redDim)))
gmm_labels = c(gmm_clusters, rep(NA, best_k)) %>% as.factor
gmm_points %>% ggplot(aes(x=PC1, y=PC2, color=gmm_labels)) + geom_point() + theme_minimal()


rm(wss, datExpr_k_means, h_clusts, cc_output, cc_output_c1, cc_output_c2, best_k, ICA_output, 
   ICA_clusters_names, signals_w_kurtosis, n_clust, gmm, gmm_points, gmm_labels, network, dend_meta, 
   best_power, c, viridis_age_cols, create_viridis_dict)

##### Plot cluterings ###################################################################

create_2D_plot = function(cat_var){
 ggplotly(plot_points %>% ggplot(aes_string(x='PC1', y='PC2', color=cat_var)) + 
          geom_point() + theme_minimal() + 
          xlab(paste0('PC1 (', round(summary(pca_output)$importance[2,1]*100,2),'%)')) +
          ylab(paste0('PC2 (', round(summary(pca_output)$importance[2,2]*100,2),'%)')))
}
create_3D_plot = function(cat_var){
  plot_points %>% plot_ly(x=~PC1, y=~PC2, z=~PC3) %>% add_markers(color=plot_points[,cat_var], size=1) %>% 
    layout(title = glue('Samples coloured by ', cat_var),
           scene = list(xaxis=list(title=glue('PC1 (',round(summary(pca_output)$importance[2,1]*100,2),'%)')),
                        yaxis=list(title=glue('PC2 (',round(summary(pca_output)$importance[2,2]*100,2),'%)')),
                        zaxis=list(title=glue('PC3 (',round(summary(pca_output)$importance[2,3]*100,2),'%)'))))  
}

# ICA
ICA_clusters %>% mutate(cl_sum=rowSums(.)) %>% as.matrix %>% melt %>% ggplot(aes(Var2,Var1)) + 
  geom_tile(aes(fill=value)) + xlab('Clusters') + ylab('Samples') + 
  theme_minimal() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + coord_flip()

# All clusterings
plot_points = datExpr_redDim %>% data.frame() %>% dplyr::select(PC1:PC3) %>%
                mutate(ID = rownames(.), subject_ID = datMeta_redDim$Subject_ID,
                cluster_km = as.factor(km_clusters), big_cluster_cc = as.factor(substr(cc_clusters,0,1)),
                cluster_cc = as.factor(cc_clusters), cluster_hc = as.factor(hc_clusters),
                cluster_ica = as.factor(ICA_clusters_min), n_cluster_ICA = as.factor(rowSums(ICA_clusters)),
                cluster_gmm = as.factor(gmm_clusters), cluster_wgcna = as.factor(wgcna_clusters),
                diagnosis = as.factor(datMeta_redDim$Diagnosis_), region = as.factor(datMeta_redDim$Brain_lobe), 
                sex = as.factor(datMeta_redDim$Sex), age = datMeta_redDim$Age)

create_2D_plot('cluster_km') ; create_2D_plot('big_cluster_cc') ; create_2D_plot('cluster_cc')
create_2D_plot('cluster_hc') ; create_2D_plot('cluster_ica')    ; create_2D_plot('cluster_gmm')
create_2D_plot(('cluster_wgcna'))
create_2D_plot('diagnosis')  ; create_2D_plot('age')
create_2D_plot('region')

create_3D_plot('cluster_km') ; create_3D_plot('cluster_cc')
create_3D_plot('cluster_hc') ; create_3D_plot('cluster_ica')
create_3D_plot('diagnosis')  ; create_3D_plot('age')


rm(create_2D_plot, create_3D_plot, ICA_clusters_names)

##### Compare Clusterings ###############################################################

# Adjusted Rand Index
clusterings = list(km_clusters, hc_clusters, substr(cc_clusters,0,1), cc_clusters, ICA_clusters_min, 
                   wgcna_clusters, gmm_clusters, datMeta_redDim$Diagnosis_)
cluster_sim = data.frame(matrix(nrow = length(clusterings), ncol = length(clusterings)))
for(i in 1:(length(clusterings))){
  cluster1 = clusterings[[i]]
  for(j in (i):length(clusterings)){
    cluster2 = clusterings[[j]]
    cluster_sim[i,j] = adj.rand.index(cluster1, cluster2)
  }
}
colnames(cluster_sim) = c('K-means','Hierarchical','Consensus L1', 'Consensus', 'ICA', 'WGCNA', 'GMM','ASD')
rownames(cluster_sim) = colnames(cluster_sim)

cluster_sim = cluster_sim %>% as.matrix %>% round(2)
heatmap.2(x = cluster_sim, Rowv = FALSE, Colv = FALSE, dendrogram = 'none', 
          cellnote = cluster_sim, notecol = 'black', trace = 'none', key = FALSE, 
          cexRow=1, cexCol=1, margins=c(7,7))


rm(i, j, cluster1, cluster2, clusterings, cluster_sim)

##### Enrichment Analysis ###############################################################
# Number of samples is too small for most clusters to be statistically significant

# Using Fisher's exact test for categorical data and Mann-Whitney for continuous (age)

enrichment_tests = function(cluster_col, signif_threshold=0.05){
  
  clusters = sort(unique(plot_points[,cluster_col]))
  
  # Diagnosis: one-tailed Fisher's exact test because we only care for ASD-related modules, not CTRL-related
  diagnosis_p_vals = clusters %>% sapply(function(c) {
    cont_table = table(plot_points[,cluster_col]==c, plot_points$diagnosis=='ASD') 
    fisher.test(cont_table, alternative='greater')$p.value
  })
  
  # Sex: two-tailed Fisher's exact test because we care both for female or male related modules
  sex_p_vals = clusters %>% sapply(function(c) {
    cont_table = table(plot_points[,cluster_col]==c, plot_points$sex) 
    fisher.test(cont_table)$p.value
  })
  
  # Region: two-tailed Fisher's exact test
  region_p_vals = clusters %>% sapply(function(c) {
    cont_table = table(plot_points[,cluster_col]==c, plot_points$region) 
    fisher.test(cont_table)$p.value
  })
  
  # Age: Mann-Whitney, like t-test (continuous vs categorical) but without the normality assumption
  age_p_vals = clusters %>% sapply(function(c) {
    wilcox.test(x = as.numeric(plot_points[,cluster_col]==c), y = plot_points$age)$p.value
  })
  
  output_df = cbind('method' = rep(cluster_col, length(age_p_vals)), 'cluster' = clusters, 
                    'diagnosis_pval' = diagnosis_p_vals, 'sex_pval' = sex_p_vals,
                    'region_pval' = region_p_vals, 'age_pval' = age_p_vals) %>% data.frame %>% 
              mutate(diagnosis_pval = as.numeric(as.character(diagnosis_pval)),
                     sex_pval = as.numeric(as.character(sex_pval)),
                     region_pval = as.numeric(as.character(region_pval)),
                     age_pval = as.numeric(as.character(age_pval)),
                     diagnosis_signif = diagnosis_pval < signif_threshold/length(clusters),
                     sex_signif = sex_pval < signif_threshold/length(clusters),
                     region_signif = region_pval < signif_threshold/length(clusters),
                     age_signif = age_pval < signif_threshold/length(clusters))
  
  return(output_df)
}
enrichment_tests_ICA = function(signif_threshold=0.05){
  
  clusters = ICA_clusters %>% colnames
  
  # Diagnosis:
  diagnosis_p_vals = clusters %>% sapply(function(c) {
    cont_table = ICA_clusters[,c] %>% table(datMeta_redDim$Diagnosis_=='ASD') 
    fisher.test(cont_table, alternative='greater')$p.value
  })
  
  # Sex:
  sex_p_vals = clusters %>% sapply(function(c) {
    cont_table = ICA_clusters[,c] %>% table(datMeta_redDim$Sex) 
    fisher.test(cont_table)$p.value
  })
  
  # Region:
  region_p_vals = clusters %>% sapply(function(c) {
    cont_table = ICA_clusters[,c] %>% table(datMeta_redDim$Brain_lobe)
    fisher.test(cont_table)$p.value
  })
  
  # Age:
  age_p_vals = clusters %>% sapply(function(c) {
    wilcox.test(x = as.numeric(ICA_clusters[,c]), y = datMeta_redDim$Age)$p.value
  })
  
  output_df = cbind('method' = rep('ICA_clusters', length(age_p_vals)), 'cluster' = clusters, 
                    'diagnosis_pval' = diagnosis_p_vals, 'sex_pval' = sex_p_vals,
                    'region_pval' = region_p_vals, 'age_pval' = age_p_vals) %>% data.frame %>% 
              mutate(diagnosis_pval = as.numeric(as.character(diagnosis_pval)),
                     sex_pval = as.numeric(as.character(sex_pval)),
                     region_pval = as.numeric(as.character(region_pval)),
                     age_pval = as.numeric(as.character(age_pval)),
                     diagnosis_signif = diagnosis_pval < signif_threshold/length(clusters),
                     sex_signif = sex_pval < signif_threshold/length(clusters),
                     region_signif = region_pval < signif_threshold/length(clusters),
                     age_signif = age_pval < signif_threshold/length(clusters))
            
  return(output_df)  
}

enrichment_results = rbind(enrichment_tests('cluster_km'), enrichment_tests('cluster_hc'), 
                           enrichment_tests('cluster_cc'), enrichment_tests_ICA())
