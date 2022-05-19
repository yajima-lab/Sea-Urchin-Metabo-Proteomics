library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(KEGGREST)
library(pathview)
#Something to consider is this really a paired t-test?
#Is the timecourse data a paired t-test




#####Proteomics Part#####
EC <- keggList("enzyme")
EC <- as.data.frame(EC)
EC$EC <- gsub("\\;.*","",EC$EC)
EC$ID <- rownames(EC)
rownames(EC) <- 1:nrow(EC)
EC <- EC[,c("ID","EC")]
EC$ID <- gsub("ec:","EC:",EC$ID)
colnames(EC) <- c("GOID","TERM")

mic_nm_file6 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Mic_NM_file6.txt")

#JUST TAKE THE RESULT from file 6 and perform analysis
#strsplit the NM and MIC columns
NM.df = do.call("rbind", strsplit(mic_nm_file6$Value.NM, ";"))
NM.df = data.frame(apply(NM.df, 2, as.numeric))
names(NM.df) = paste("NM.", 1:3, sep = "")
head(NM.df)

MIC.df = do.call("rbind", strsplit(mic_nm_file6$Value.MIC, ";"))
MIC.df = data.frame(apply(MIC.df, 2, as.numeric))
names(MIC.df) = paste("MIC.", 1:3, sep = "")
head(MIC.df)

mic_nm_full <- cbind(mic_nm_file6, MIC.df, NM.df)

#Get Log2FC for Mic-NM
col <- c("MIC.1","NM.1","MIC.2","NM.2","MIC.3","NM.3") 
#log(x/y) = log(x) - log(y)
mic_nm_full$logFC.Mic.m.NM <- apply(mic_nm_full[,col], 1, function(x) log2(mean(x[c("MIC.1","MIC.2","MIC.3")])) - log2(mean(x[c("NM.1","NM.2","NM.3")])))
mic_nm_full$logFC.NM.m.Mic <- apply(mic_nm_full[,col], 1, function(x) log2(mean(x[c("NM.1","NM.2","NM.3")])) - log2(mean(x[c("MIC.1","MIC.2","MIC.3")])))


#GSEA uses all of the genes not just DE genes
FUNCt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/funcanno_t2g.csv", row.names = 1)
FUNCt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/funcanno_t2n.csv", row.names = 1)
GOt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GOanno_t2g.csv", row.names = 1)
GOt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GOanno_t2n.csv", row.names = 1)
KEGGt2g <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGGanno_t2g.csv", row.names = 1)
KEGGt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGGanno_t2n.csv", row.names = 1)
KEGGt2n$PathwayID <- gsub("path:","",KEGGt2n$PathwayID)
GOt2n <- rbind(GOt2n, EC)


mic_nm_full$pval.NM.vs.MIC.ttest <- as.numeric(mic_nm_full$pval.NM.vs.MIC.ttest)
mic_nm_full$GSEA.Mic.m.NM <- apply(mic_nm_full[,c("pval.NM.vs.MIC.ttest", "logFC.Mic.m.NM")],1,function(x) -log10(x[c("pval.NM.vs.MIC.ttest")])/sign(x[c("logFC.Mic.m.NM")]))
mic_nm_full$GSEA.NM.m.Mic <- apply(mic_nm_full[,c("pval.NM.vs.MIC.ttest","logFC.NM.m.Mic")],1,function(x) -log10(x[c("pval.NM.vs.MIC.ttest")])/sign(x[c("logFC.NM.m.Mic")]))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

mic_nm_full <- completeFun(mic_nm_full, c("GSEA.Mic.m.NM", "GSEA.NM.m.Mic"))


NCBItoGENEID <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/NCBItoGENEPAGE.txt", header = F)
GENEIDtoGENEPAGE <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GENEIDtoGENEPAGE.txt", header = F)
NCBItoGENEID <- subset(NCBItoGENEID, select = c("V2","V3"))
GENEIDtoGENEPAGE <- subset(GENEIDtoGENEPAGE, select = c("V1","V3"))

mic_nm_full$Protein.ID<- gsub("\\..*","",mic_nm_full$Protein.ID)

colnames(NCBItoGENEID) <- c("Protein.ID","GENEID")
colnames(GENEIDtoGENEPAGE) <- c("GENEPAGE","GENEID")
mic_nm_full <- merge(mic_nm_full, NCBItoGENEID, by = "Protein.ID" )
mic_nm_full <- merge(mic_nm_full, GENEIDtoGENEPAGE, by = "GENEID")

mic_nm_full <- mic_nm_full[!duplicated(mic_nm_full$GENEPAGE), ]

#make GSEA genelist
MIC.m.NM <-  mic_nm_full$GSEA.Mic.m.NM
names(MIC.m.NM) <- mic_nm_full$GENEPAGE
MIC.m.NM <- sort(MIC.m.NM, decreasing = TRUE)

NM.m.MIC <- mic_nm_full$GSEA.NM.m.Mic
names(NM.m.MIC) <- mic_nm_full$GENEPAGE
NM.m.MIC <- sort(NM.m.MIC, decreasing = TRUE)
#####Proteomics ANNOTATION####
####PATHWAY MAPS Annotation
annotation_matrix <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/annotation_matrix.txt")
genename_to_ENTREZID <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GeneExternalRef.txt", header = F)
colnames(genename_to_ENTREZID) <- c("GENE.ID","LOC","Descrp","org.code","ENTREZ")
genename_to_ENTREZID <- subset(genename_to_ENTREZID, select = c("GENE.ID","ENTREZ"))
full_annotation_matrix <- merge(annotation_matrix, genename_to_ENTREZID, by = "GENE.ID")

#Conversion to ENTREZ
MIC.m.NM.ENTREZ <-  mic_nm_full$GSEA.Mic.m.NM
names(MIC.m.NM.ENTREZ) <- mic_nm_full$GENEPAGE
MIC.m.NM.convert <- as.data.frame(MIC.m.NM.ENTREZ)
MIC.m.NM.convert$GENEPAGE.ID <- rownames(MIC.m.NM.convert)
MIC.m.NM.convert2 <- merge(MIC.m.NM.convert, full_annotation_matrix, by = "GENEPAGE.ID")
MIC.m.NM.ENTREZ <- subset(MIC.m.NM.convert2, select  = c("MIC.m.NM.ENTREZ","ENTREZ"))
MIC.m.NM.ENTREZ <- unique(MIC.m.NM.ENTREZ)
MIC.m.NM.pathways <- MIC.m.NM.ENTREZ$MIC.m.NM.ENTREZ
names(MIC.m.NM.pathways) <- MIC.m.NM.ENTREZ$ENTREZ
MIC.m.NM.pathways <- sort(MIC.m.NM.pathways, decreasing = TRUE)


####Save CSV for app####
write.csv(MIC.m.NM.pathways, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/MIC NM/mic.m.nm.proteomics.csv")