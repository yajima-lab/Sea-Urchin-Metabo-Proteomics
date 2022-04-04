#Timecourse
library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(KEGGREST)
library(ggrepel)
library(Mfuzz)
library(stringr)
#BiocManager::install("marray")
library(marray)

#for PCA plotting
library(ggforce)
library(scales)

#heatmap creation
library(ComplexHeatmap)
library(pathview) #visualizationa


##https://stackoverflow.com/questions/15688758/r-stats-citation-for-a-scientific-paper
citations <- function(includeURL = TRUE, includeRStudio = TRUE) {
  if(includeRStudio == TRUE) {
    ref.rstudio <- RStudio.Version()$citation
    if(includeURL == FALSE) {
      ref.rstudio$url <- NULL;
    }
    print(ref.rstudio, style = 'text')
    cat('\n')
  }
  
  cit.list <- c('base', names(sessionInfo()$otherPkgs))
  for(i in 1:length(cit.list)) {
    ref <- citation(cit.list[i])
    if(includeURL == FALSE) {
      ref$url <- NULL;
    }
    print(ref, style = 'text')
    cat('\n')
  }
}


citations()