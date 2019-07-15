
setwd('/afs/inf.ed.ac.uk/user/s17/s1725186/Documents/PhD-InitialExperiments/FirstYearReview/R_Markdowns')

library(tidyverse) ; library(Rtsne)

load('./../Data/Gandal/preprocessed_data.RData')
datExpr = datExpr %>% data.frame
DE_info = DE_info %>% data.frame

perplexities = c(1,2,5,10,50,100)

for(i in 1:length(perplexities)){
  print(i)
  set.seed(123)
  tsne = datExpr %>% Rtsne(perplexity=perplexities[i])
  tsne_coords = tsne$Y
  write.csv(tsne_coords, paste0('./../Data/Gandal/Visualisations/tsne_perplexity_',perplexities[i],'.csv'), row.names=F)
}