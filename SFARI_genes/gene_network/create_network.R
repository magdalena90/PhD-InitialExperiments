
setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/SFARI_genes/gene_network')

library(tidyverse)
library(plotly)
library(visNetwork)
library(Matrix)

#########################################################################################################
# Select if results should be filtered to ASD related or not

filter_ASD = FALSE

#########################################################################################################
# Get list of genes and pubmed IDs

# Genes
list_gene_files = list.files('./SFARI_scraping/genes/')

# Pubmed IDs
get_pubmed_IDs = function(autism=FALSE){
  all_pubmed_ids = c()
  n = 0
  
  for(gene_file in list_gene_files){
    file = read.delim(paste0('./SFARI_scraping/genes/', gene_file))
    
    if(autism) file = file %>% filter(Autism.Report == 'Yes')
    
    if(nrow(file)>0){
      pubmed_ids = unique(file$pubmed.ID)
      all_pubmed_ids = unique(c(all_pubmed_ids, pubmed_ids)) 
    }
    
    n = c(n, length(all_pubmed_ids))
  }
  
  all_pubmed_ids = as.character(all_pubmed_ids)
  
  return(list('n'=n, 'all_pubmed_ids'=all_pubmed_ids))
}

get_pubmed_IDs_output = get_pubmed_IDs(filter_ASD)
n = get_pubmed_IDs_output$n
all_pubmed_ids = get_pubmed_IDs_output$all_pubmed_ids

# See change in number of pubmed articles by each new gene (Seems like the number of new articles were beginning to decrease in the end)
pubmed_counts = data.frame('genes' = seq(1,length(n)), 'articles'=n)

pubmed_counts %>% ggplot(aes(genes, articles,color='#434343')) + geom_line() +
  geom_smooth(method='lm', color='#666666', size=0.5) + theme_minimal() + theme(legend.position='none') +
  ggtitle('Increase in number of pubmed articles by each new gene')

rm(n, get_pubmed_IDs_output)
#########################################################################################################
# Create gene - pubmedID sparse matrix

SFARI_genes_info = read.csv('./../Data/SFARI_genes_01-15-2019.csv') %>% mutate('ID' = gene.symbol)

create_gene_pmID_mat = function(autism=FALSE){
  
  # Create empty dataframe
  gene_pmID_mat = data.frame(matrix(0, nrow=length(list_gene_files), ncol=length(all_pubmed_ids)))
  rownames(gene_pmID_mat) = sapply(list_gene_files, function(x) gsub('.csv', '', x))
  colnames(gene_pmID_mat) = all_pubmed_ids
  
  # Fillout dataframe
  for(gene_file in list_gene_files){
    file = read.delim(paste0('./SFARI_scraping/genes/', gene_file))
    
    if(autism) file = file %>% filter(Autism.Report == 'Yes')
    
    if(nrow(file)>0){
      pubmed_ids = as.character(unique(file$pubmed.ID))
      gene_pmID_mat[gsub('.csv', '', gene_file), pubmed_ids] = 1 
    }
  }
  
  gene_names = rownames(gene_pmID_mat)
  sparse_gene_pmID_mat = Matrix(as.matrix(gene_pmID_mat), sparse=TRUE)
    
  return(list('gene_names'=gene_names, 'sparse_gene_pmID_mat'=sparse_gene_pmID_mat))
}

gene_pmID_mat_output = create_gene_pmID_mat(filter_ASD)
gene_names = gene_pmID_mat_output$gene_names
gene_pmID_mat = gene_pmID_mat_output$sparse_gene_pmID_mat

plot_data = data.frame('ID' = rownames(gene_pmID_mat), 'n_papers'=rowSums(gene_pmID_mat)) %>% 
            left_join(SFARI_genes_info, by='ID') %>% mutate(gene.score=as.factor(gene.score))

ggplotly(plot_data %>% ggplot(aes(n_papers)) + geom_histogram(bins=max(plot_data$n_papers), alpha=0.6) + 
         theme_minimal() + ggtitle('Distribution of number of papers per gene'))

ggplotly(plot_data %>% ggplot(aes(gene.score, n_papers, fill=gene.score)) + geom_boxplot() + theme_minimal() +
              ggtitle('Number of papers per gene'))

# Remove genes with no articles related to them
gene_pmID_mat = gene_pmID_mat[rowSums(gene_pmID_mat)>0,]

rm(plot_data, gene_pmID_mat_output)
#########################################################################################################
# Create edge and node lists 

# SFARI gene score colors:
# Syndromic: #DC2641, S1: #FF7631, S2: #FFB100, S3: #E8E328, S4: #8CC83F, S5: #62CCA6, S6: #59B9C9

# Nodes
nodes = SFARI_genes_info %>% dplyr::select(gene.symbol, gene.score, syndromic, number.of.reports) %>% 
        mutate(id = gene.symbol, title=gene.symbol, value=number.of.reports,
               shape=ifelse(syndromic==0, 'circle','square'),
               color=case_when(gene.score==1 ~ '#FF7631', gene.score==2 ~ '#FFB100',
                               gene.score==3 ~ '#E8E328', gene.score==4 ~ '#8CC83F',
                               gene.score==5 ~ '#62CCA6', gene.score==6 ~ '#59B9C9',
                               is.na(gene.score) ~ '#b3b3b3')) %>%
        filter(!duplicated(id))


# Edges
edges = tibble('from'=character(), 'to'=character(), 'weight'=integer())

gene_gene_mat = gene_pmID_mat %*% t(gene_pmID_mat)
edges = summary(gene_gene_mat) %>% rename('from'=i, 'to'=j, 'weight'=x) %>% 
        mutate(from = gene_names[from], to = gene_names[to]) %>%
        filter(from!=to)

nodes = nodes %>% filter(gene.symbol %in% unique(c(edges$from, edges$to)))

# save(nodes, edges, file='nodes_and_edges_ASD.RData')

rm(gene_gene_mat)
#########################################################################################################
# Plot network

visNetwork(nodes, edges, height = '700px', width='100%') %>% visOptions(selectedBy = 'gene.score') %>% 
  visEdges(smooth = TRUE) %>% visIgraphLayout() %>% visLayout(randomSeed=12)

edges_filtered = edges %>% filter(weight>1)
nodes_filtered = nodes %>% filter(gene.symbol %in% unique(c(edges_filtered$from, edges_filtered$to)))

visNetwork(nodes_filtered, edges_filtered, height = '700px', width='100%') %>% 
           visOptions(selectedBy='gene.score') %>% visEdges(smooth=TRUE) %>% 
           visIgraphLayout() %>% visLayout(randomSeed=12)

#########################################################################################################
# Network analysis

graph = graph_from_data_frame(edges, vertices=nodes, directed=FALSE)

# Eigencentrality
eigen_centrality_output = eigen_centrality(graph, weights=edges$weight)
eigencentr = eigen_centrality_output$vector

# Strong correlation (negative because 'higher' SFARI scores have actually lower numbers)
score_eigencentr_cor = cor(nodes$gene.score[!is.na(nodes$gene.score)], eigencentr[!is.na(nodes$gene.score)])

# Clustering





