library(KEGGREST)
library(clusterProfiler)
library(enrichplot)


###Pathway input files + Small Conversions###
EC <- keggList("enzyme")
EC <- as.data.frame(EC)
EC$EC <- gsub("\\;.*","",EC$EC)
EC$ID <- rownames(EC)
rownames(EC) <- 1:nrow(EC)
EC <- EC[,c("ID","EC")]
EC$ID <- gsub("ec:","EC:",EC$ID)
colnames(EC) <- c("GOID","TERM")
#GSEA uses all of the genes not just DE genes
FUNCt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/funcanno_t2g.csv", row.names = 1)
FUNCt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/funcanno_t2n.csv", row.names = 1)
GOt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GOanno_t2g.csv", row.names = 1)
GOt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GOanno_t2n.csv", row.names = 1)
KEGGt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGGanno_t2g.csv", row.names = 1)
KEGGt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGGanno_t2n.csv", row.names = 1)
KEGGt2n$PathwayID <- gsub("path:","",KEGGt2n$PathwayID)
GOt2n <- rbind(GOt2n, EC)
#Convert FUNC list to general FUNC
FUNCt2n <- unique(FUNCt2n)
FUNCver2_t2g <- merge(FUNCt2g, FUNCt2n, by = "Class.L3")
FUNCver2_t2g <- FUNCver2_t2g[,c("Class.L1","V1")]




#Get genelist in form of GENEPAGE-ID, use the merge or subset function if you have NCBI ID's 
annotation_matrix <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/annotation_matrix.txt")


#Only need to keep 4 lines of code below if you have mRNA transcript ID's (transcriptomics or other projects)
NcbiMrna <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/NcbiMrnaEchinobaseGene_Spur.txt", header = F)
colnames(NcbiMrna) <- c("gi","NcbiMrna.ID","GENE.ID","Symbol")
NcbiMrna <- subset(NcbiMrna, select = c("GENE.ID","NcbiMrna.ID"))
full_annotation_matrix <- merge(annotation_matrix, NcbiMrna, by = "GENE.ID")

#write.csv(full_annotation_matrix, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/VASA MO/annotations.csv")





#example of how to do ORA 
example_NCBIID <- c("XM_786251.5",	"XM_011672995.2",	"XM_030980878.1","XM_003729355.3",
                    "XM_030991022.1","XM_030972690.1","XM_030986685.1","XM_030987698.1") #example list for practice
example_NCBIID <- gsub("\\..*","",example_NCBIID)


#define genelist for GENEPAGE.ID
genelist <- unique(full_annotation_matrix[full_annotation_matrix$NcbiMrna.ID %in% example_NCBIID,]$GENEPAGE.ID)



#Make a function
 ORA_analysis <- function(genelist,FUNCt2g, FUNCver2_t2g, GOt2g, GOt2n,KEGGt2g, KEGGt2n){

   #You can change the pvalueCutoff, and pAdjustMethod
   funcver1 = enricher(genelist, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = FUNCt2g)
   funcver2 = enricher(genelist, pvalueCutoff = 1,  TERM2GENE = FUNCver2_t2g)
   GO =  enricher(genelist, pvalueCutoff = 1,  TERM2GENE = GOt2g, TERM2NAME = GOt2n)
   KEGG =  enricher(genelist, pvalueCutoff = 1,  TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n)
   list_ORA <- list(funcver1, funcver2, GO, KEGG)
   return(list_ORA)
 }

 
#Execute function
  x <- ORA_analysis(genelist,FUNCt2g, FUNCver2_t2g, GOt2g, GOt2n,KEGGt2g, KEGGt2n)
 
#Example of plotting, more plotting functions at clusterprofiler website
  #x[[1]] ~ funcver1
  #x[[2]] ~ funcver2
  #x[[3]] ~ GO
  #x[[4]] ~ KEGG
  
  dotplot(x[[1]]) + ggtitle("FUNC: Title INPUT HERE") #<- should have error b/c no enrichment in this example
  dotplot(x[[2]]) + ggtitle("FUNC: Title INPUT HERE") #<- should have error b/c no enrichment in this example
 dotplot(x[[3]]) + ggtitle("GO: Title INPUT HERE")
 dotplot(x[[4]]) + ggtitle("KEGG: Title INPUT HERE")