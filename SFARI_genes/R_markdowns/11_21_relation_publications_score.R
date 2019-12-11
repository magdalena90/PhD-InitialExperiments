
library(dplyr)
library(easyPubMed)
library(biomaRt)

# Load datExpr to get full list of Ensembl IDs
datExpr = read.csv('~/Documents/PhD-InitialExperiments/FirstYearReview/Data/Gandal/RNAseq_ASD_datExpr.csv')

# Get gene names from Ensembl IDs
getinfo = c('ensembl_gene_id','external_gene_id')
mart = useMart(biomart='ENSEMBL_MART_ENSEMBL',  dataset='hsapiens_gene_ensembl', host='feb2014.archive.ensembl.org')
datGenes = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), values=datExpr$X, mart=mart)
datGenes = datGenes[match(datExpr$X, datGenes$ensembl_gene_id),]

# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[10001:20000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_20000.csv', row.names=FALSE)




# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[25001:30000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_35000.csv', row.names=FALSE)




# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[30001:35000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_40000.csv', row.names=FALSE)




# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[35001:40000], 'Count' = 0)
for(i in 3895:3900){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_44000.csv', row.names=FALSE)





# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[38901:40000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_45000.csv', row.names=FALSE)





# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[40001:42500], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_50000.csv', row.names=FALSE)




# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[42501:45000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_52000.csv', row.names=FALSE)




# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[45001:50000], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_55000.csv', row.names=FALSE)






# Get pubmed count from gene names
pubmed_count_df = data.frame('Gene' = unique(datGenes$external_gene_id)[50001:length(unique(datGenes$external_gene_id))], 'Count' = 0)
for(i in 1:nrow(pubmed_count_df)){
  pubmed_count_df$Count[i] = get_pubmed_ids(paste0(pubmed_count_df$Gene[i],'[All Fields]'))$Count
  if(i%%100==0){
    cat(paste0(i, '/', nrow(pubmed_count_df),' genes processed\n'))
    cat(paste0(pubmed_count_df$Gene[i],': ', pubmed_count_df$Count[i],'\n\n'))
  }
}

# Save results
write.csv(pubmed_count_df, file='~/Documents/PhD-InitialExperiments/SFARI_genes/Data/pubmed_count_by_gene_60000.csv', row.names=FALSE)
