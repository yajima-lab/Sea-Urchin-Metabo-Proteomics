GSEA_stat <- function(timecourse_full, cond1, cond2, pval_cutoff){
  #GSEA stat works by cond2 - cond1
  #cond2 and cond1 must be string lists "character vector"
  #LogFC
  
  col <- c("Egg.1", "Egg.2", "Egg.3",
           "2cell.1", "2cell.2", "2cell.3",
           "16cell.1", "16cell.2", "16cell.3",
           "Morula.1", "Morula.2", "Morula.3",
           "Blastula.1", "Blastula.2", "Blastula.3",
           "Gastrula.1", "Gastrula.2", "Gastrula.3",
           "Prism.1", "Prism.2", "Prism.3",
           "Pluteus.1",  "Pluteus.2",  "Pluteus.3") 
  
  timecourse_full$logFC.cond2.m.cond1 <- apply(timecourse_full[,col], 1, function(x) log2(mean(x[cond2])) - log2(mean(x[cond1])))
  
  #ttest
  timecourse_full$ttest.cond2.m.cond1 <- apply(timecourse_full[,col], 1, function(x) {t.test(x[cond2], x[cond1], paired = TRUE)$p.value})
  
  #GSEA
  timecourse_full$GSEA.cond2.m.cond1 <- apply( timecourse_full[,c("ttest.cond2.m.cond1", "logFC.cond2.m.cond1")],1,function(x) -log10(x[c("ttest.cond2.m.cond1")])/sign(x[c("logFC.cond2.m.cond1")]))
  
  
  #Add GSEA to timecourse
  timecourse_full_pathway <- subset(timecourse_full, select = c("GSEA.cond2.m.cond1", "Protein.ID"))
  
  
  #Annotate the data, Conversion of IDs
  NCBItoGENEID <- read.delim(file = "Data/Anno/NCBItoGENEPAGE.txt", header = F)
  GENEIDtoGENEPAGE <- read.delim(file = "Data/Anno/GENEIDtoGENEPAGE.txt", header = F)
  NCBItoGENEID <- subset(NCBItoGENEID, select = c("V2","V3"))
  GENEIDtoGENEPAGE <- subset(GENEIDtoGENEPAGE, select = c("V1","V3"))
  
  timecourse_full$Protein.ID<- gsub("\\..*","",timecourse_full$Protein.ID)
  
  colnames(NCBItoGENEID) <- c("Protein.ID","GENEID")
  colnames(GENEIDtoGENEPAGE) <- c("GENEPAGE","GENEID")
  timecourse_full <- merge(timecourse_full, NCBItoGENEID, by = "Protein.ID" )
  timecourse_full <- merge(timecourse_full, GENEIDtoGENEPAGE, by = "GENEID")
  
  timecourse_full <- timecourse_full[!duplicated(timecourse_full$GENEPAGE), ]
  
  
  #Pathway annotation Files
  KEGGt2g <- read.csv(file = "Data/Anno/KEGGt2g.csv", row.names = 1)
  KEGGt2n <- read.csv(file = "Data/Anno/KEGGt2n.csv", row.names = 1)
  GOt2g <- read.csv(file = "Data/Anno/GOt2g.csv", row.names = 1)
  GOt2n <- read.csv(file = "Data/Anno/GOt2n.csv", row.names = 1)
  FUNCt2g <- read.csv(file = "Data/Anno/FUNCt2g.csv", row.names = 1)
  FUNCt2n <- read.csv(file = "Data/Anno/FUNCt2n.csv", row.names = 1)
  FUNCver2_t2g <- read.csv(file = "Data/Anno/FUNCver2_t2g.csv", row.names = 1)
  
  
  #make GSEA genelists
  cond2.m.cond1 <- timecourse_full$GSEA.cond2.m.cond1
  names(cond2.m.cond1) <- timecourse_full$GENEPAGE
  cond2.m.cond1[cond2.m.cond1 == Inf] <- 0
  cond2.m.cond1 <- sort(cond2.m.cond1, decreasing = TRUE)
  
  gsea.KEGG = GSEA(cond2.m.cond1, pvalueCutoff = pval_cutoff, TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, verbose = FALSE)
  gsea.FUNC = GSEA(cond2.m.cond1, pvalueCutoff = pval_cutoff,  TERM2GENE = FUNCt2g, verbose=FALSE)
  gsea.FUNCver2 = GSEA(cond2.m.cond1, pvalueCutoff = pval_cutoff,  TERM2GENE = FUNCver2_t2g, verbose=FALSE)
  gsea.GO = GSEA(cond2.m.cond1, pvalueCutoff = pval_cutoff,  TERM2GENE = GOt2g, TERM2NAME = GOt2n, verbose=FALSE)
  #Make sure NA's are deleted from set
  gsea.GO@result <- gsea.GO@result[complete.cases(gsea.GO@result), ]
  
  #Return list of results
  gsea.ALL = list(gsea.KEGG, gsea.FUNC, gsea.FUNCver2, gsea.GO)
  return(gsea.ALL)
}
gsea_timepoints <- list(c("Egg.1", "Egg.2", "Egg.3"), c("2cell.1", "2cell.2", "2cell.3"),
                        c("16cell.1", "16cell.2", "16cell.3"), c("Morula.1", "Morula.2", "Morula.3"),
                        c("Blastula.1", "Blastula.2", "Blastula.3"), c("Gastrula.1", "Gastrula.2", "Gastrula.3"),
                        c("Prism.1", "Prism.2", "Prism.3"), c("Pluteus.1",  "Pluteus.2",  "Pluteus.3"))
BlastulamMorula <- GSEA_stat(timecourse_full,gsea_timepoints[[4]], gsea_timepoints[[5]], 0.05)