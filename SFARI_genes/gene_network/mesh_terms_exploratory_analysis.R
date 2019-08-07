
setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/SFARI_genes/gene_network')

library(tidyverse)
library(jsonlite)
library(viridisLite)

######################################################################################################
# Read json with MeSH terms by paper

mesh_terms_by_article_list = read_json('./Data/MeSH_terms_by_paper_.json', simplifyVector = TRUE)

######################################################################################################
# Get list of all mesh terms

all_mesh_terms = c()
n = 0

for(mesh_terms_in_article in mesh_terms_by_article_list){
  all_mesh_terms = unique(c(all_mesh_terms, mesh_terms_in_article))
  n = c(n, length(all_mesh_terms))
}

mesh_terms_counts = data.frame('articles' = seq(1,length(n)), 'mesh_terms'=n)

ggplotly(mesh_terms_counts %>% ggplot(aes(articles, mesh_terms, color='#434343')) + geom_line() +
         geom_smooth(method='lm', color='#666666', se=FALSE, size=0.5) + theme_minimal() + theme(legend.position='none') +
         ggtitle('Increase in number of MeSH terms by each new article'))


rm(mesh_terms_in_article, n, mesh_terms_counts)
######################################################################################################
# Create article-MeSH term (sparse) matrix

all_articles = names(mesh_terms_by_article_list)

article_mesh_term_mat = data.frame(matrix(0, nrow=length(all_articles), ncol=sum(!is.na(all_mesh_terms))))
rownames(article_mesh_term_mat) = all_articles
colnames(article_mesh_term_mat) = all_mesh_terms[!is.na(all_mesh_terms)]

for(article in all_articles){
  article_mesh_terms = mesh_terms_by_article_list[[article]]
  if(sum(is.na(article_mesh_terms))==0){
    article_mesh_term_mat[article, article_mesh_terms] = 1
  }
}

sparse_article_mesh_term_mat = Matrix(as.matrix(article_mesh_term_mat), sparse=TRUE)

mesh_terms_by_article = data.frame('id' = rownames(article_mesh_term_mat), 'n_mesh_terms' = rowSums(article_mesh_term_mat))
bins = max(mesh_terms_by_article$n_mesh_terms)-min(mesh_terms_by_article$n_mesh_terms)+1
ggplotly(mesh_terms_by_article %>% ggplot(aes(n_mesh_terms)) + geom_histogram(bins=bins, alpha=0.6) + 
              ggtitle('Distribution of number of MeSH terms by article') + theme_minimal())

articles_by_mesh_term = data.frame('id' = colnames(article_mesh_term_mat), 'n_articles' = colSums(article_mesh_term_mat))
bins = max(articles_by_mesh_term$n_articles)-min(articles_by_mesh_term$n_articles)+1
ggplotly(articles_by_mesh_term %>% ggplot(aes(n_articles)) + geom_histogram(bins=bins, alpha=0.6) + 
                ggtitle('Distribution of number of articles that contain each MeSH term') + theme_minimal())

top_n_counts = 10
top_mesh_terms = articles_by_mesh_term %>% top_n(top_n_counts, wt=n_articles)
ggplotly(top_mesh_terms %>% ggplot(aes(reorder(id, n_articles), n_articles, fill=n_articles)) + 
         geom_bar(stat='identity') + scale_fill_viridis_c() + coord_flip() + 
         ggtitle(paste0('MeSH terms contained in at least ', top_n_counts, ' articles')) + 
         xlab('') + ylab('No. Articles') + theme_minimal() + theme(legend.position='none'))


rm(article, article_mesh_terms, mesh_terms_by_article, articles_by_mesh_term)
######################################################################################################








