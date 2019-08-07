
setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/FirstYearReview/SFARI_genes/')

library(tidyverse)
library(reutils)
library(jsonlite)

######################################################################################################
# Get list of PubMed IDs

# Get genes
list_gene_files = list.files('./SFARI_scraping_env/SFARI_scraping/genes/')

# Get Pubmed IDs
all_pubmed_ids = c()
  
for(gene_file in list_gene_files){
  file = read.delim(paste0('./SFARI_scraping_env/SFARI_scraping/genes/', gene_file))
  pubmed_ids = unique(file$pubmed.ID)
  all_pubmed_ids = unique(c(all_pubmed_ids, pubmed_ids)) 
}
  
all_pubmed_ids = as.character(all_pubmed_ids)


rm(file, gene_file)

######################################################################################################
# Get year of publication for each paper
year_of_paper = list()
n=1

print('Retrievieng year of publication from all pubmed IDs...')

for(pubmed_id in all_pubmed_ids){
  paper_info = efetch(pubmed_id, db='pubmed')
  year = paper_info$xmlValue('//PubDate') %>% substr(1,4) %>% as.numeric
  year_of_paper[[pubmed_id]] = year
  
  if(n%%10==0) print(n)
  n=n+1
}

year_of_paper_df = data.frame('pubmedID' = names(year_of_paper),
                              'year' = unlist(year_of_paper, use.names = FALSE))

rm(n, pubmed_id, paper_info, year)

write.csv(year_of_paper_df, './../Data/SFARI/year_by_paper.csv', row.names = FALSE)

######################################################################################################
# Get MeSH terms for each paper

print('Retrievieng MeSH terms from all pubmed IDs...')

MeSH_terms_by_paper = list()
n=1

for(pubmed_id in all_pubmed_ids){
  paper_info = efetch(pubmed_id, db='pubmed')
  paper_MeSH_terms = paper_info$xmlValue('//MeshHeading')
  MeSH_terms_by_paper[[pubmed_id]] = paper_MeSH_terms

  if(n%%10==0) print(n)
  n=n+1
}

rm(pubmed_id, paper_info)

# Write json with MeSH terms by paper
write_json(MeSH_terms_by_paper, './../Data/SFARI/MeSH_terms_by_paper.json')
