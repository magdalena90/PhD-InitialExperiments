
library(biomaRt)
library(GO.db)

setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/FirstYearReview/Gandal')

datExpr = read.csv('./../Data/Gandal/RNAseq_ASD_datExpr.csv', row.names=1)
datExpr = datExpr[rowSums(datExpr)>0,]
gene_ids = rownames(datExpr)[substr(rownames(datExpr),1,3)=='ENS']

mart = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl',
               host='feb2014.archive.ensembl.org') ## Gencode v19

datProbes = getBM(attributes=c('ensembl_gene_id','go_id'), filters=c('ensembl_gene_id'), 
                  values=gene_ids, mart=mart)
datProbes = datProbes %>% filter(go_id!='')

GO_id_term_map = data.frame('go_id' = unique(datProbes$go_id),
                            'go_term' = unname(Term(unique(datProbes$go_id))))

GO_annotations = datProbes %>% left_join(GO_id_term_map, by='go_id')

write.csv(GO_annotations, file='./../Data/GO_annotations/genes_GO_annotations.csv', row.names = FALSE)
