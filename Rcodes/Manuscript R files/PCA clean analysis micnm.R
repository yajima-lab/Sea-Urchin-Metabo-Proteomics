library(pca3d)
library(factoextra)
library(PCAtools)
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
library(pathview)
library(ggiraph)
library(ggfortify)
library(cluster)
library(tidyverse)

mic_nm_file3 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Mic_NM_file3.txt")
mic_nm_file4 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Mic_NM_file4.txt")
mic_nm_file5 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Mic_NM_file5.txt")
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


##########PCA with choosing most variable amt of proteins###########
mic_nm_samples <- as.data.frame(t(mic_nm_full[,c("MIC.1","NM.1","MIC.2","NM.2","MIC.3","NM.3")]))


rv <- rowVars(as.matrix(t(mic_nm_samples)))
hist(rv)
select = order(rv, decreasing=TRUE)[seq_len(4683)]
#most_variable_proteins <- timecourse_samples[ ,select]
pca_mv = prcomp(mic_nm_samples[ ,select], center = TRUE, scale = TRUE)
x <- summary(pca_mv)
per_var_exp <- as.data.frame(x$importance)
#PCA plot
pca.mv_plot <-  as.data.frame(as.matrix(pca_mv$x))
pca.mv_plot$stages <- rownames(pca.mv_plot)
pca.mv_plot$cell <- gsub("\\..*","",pca.mv_plot$stages)

#For app idea (PC1 and PC2 can be changed on x and y axis etc...)
library(ggforce)
ggplot(pca.mv_plot,aes(x=PC1,y=PC2,col=cell, label = stages))+
  geom_mark_ellipse(aes(fill = cell, label = cell)) +
  geom_point()+ 
  geom_text()+
  xlab(paste0("PC1 - percent variance explained: ",per_var_exp[2,c("PC1")])) +
  ylab(paste0("PC2 - percent variance explained: ",per_var_exp[2,c("PC2")])) 
