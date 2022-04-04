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
#####Prot Timecourse####
timecourse_file6 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/timecourse_file6.txt")
colnames(timecourse_file6)
colnames(timecourse_file6) <- c("Analyte.","Protein.ID","Entry.Name","Protein.Description",        
                                "Gene","Contaminants","Protein.Group","Sub.Group",                   
                                "Number.of..Used.Peptide.Ions", "Egg", "2cell","16cell",
                                "Morula","Blastula","Gastrula","Prism","Pluteus", "Numb.Egg",
                                "Numb.2cell","Numb.16cell","Numb.Morula","Numb.Blastula",
                                "Numb.Gastrula","Numb.Prism","Numb.Pluteus","Median.Egg",
                                "Median.2cell","Median.16cell","Median.Morula","Median.Blastula",
                                "Median.Gastrula","Median.Prism","Median.Pluteus","Avg.Egg",
                                "Avg.2cell","Avg.16cell","Avg.Morula","Avg.Blastula",
                                "Avg.Gastrula","Avg.Prism","Avg.Pluteus","Std.Egg",
                                "Std.2cell","Std.16cell","Std.Morula","Std.Blastula",
                                "Std.Gastrula","Std.Prism","Std.Pluteus","CV.Egg",
                                "CV.2cell","CV.16cell","CV.Morula","CV.Blastula",
                                "CV.Gastrula","CV.Prism","CV.Pluteus","Multiple","Eggv2cell","X.1.vs..3",                   
                                "X.1.vs..4","X.1.vs..5", "X.1.vs..6", "X.1.vs..7",                   
                                "X.1.vs..8" ,"2cellv16cell","X.2.vs..4","X.2.vs..5",                   
                                "X.2.vs..6","X.2.vs..7", "X.2.vs..8","16cellvMorula",                 
                                "X.3.vs..5","X.3.vs..6","X.3.vs..7","X.3.vs..8",                 
                                "MorulavBlastula","X.4.vs..6","X.4.vs..7","X.4.vs..8",                  
                                "BlastulavGastrula","X.5.vs..7", "X.5.vs..8","GastrulavPrism",                   
                                "X.6.vs..8","PrismvPluteus")




#Functions
#strsplit the time columns
Egg.df = do.call("rbind", strsplit(timecourse_file6$Egg, ";"))
Egg.df = data.frame(apply(Egg.df, 2, as.numeric))
names(Egg.df) = paste("Egg.", 1:3, sep = "")
head(Egg.df)


twocell.df = do.call("rbind", strsplit(timecourse_file6$"2cell", ";"))
twocell.df = data.frame(apply(twocell.df, 2, as.numeric))
names(twocell.df) = paste("2cell.", 1:3, sep = "")
head(twocell.df)

cell16.df = do.call("rbind", strsplit(timecourse_file6$"16cell", ";"))
cell16.df = data.frame(apply(cell16.df, 2, as.numeric))
names(cell16.df) = paste("16cell.", 1:3, sep = "")
head(cell16.df)

Morula.df = do.call("rbind", strsplit(timecourse_file6$Morula, ";"))
Morula.df = data.frame(apply(Morula.df, 2, as.numeric))
names(Morula.df) = paste("Morula.", 1:3, sep = "")
head(Morula.df)

Blastula.df = do.call("rbind", strsplit(timecourse_file6$Blastula, ";"))
Blastula.df = data.frame(apply(Blastula.df, 2, as.numeric))
names(Blastula.df) = paste("Blastula.", 1:3, sep = "")
head(Blastula.df)

Gastrula.df = do.call("rbind", strsplit(timecourse_file6$Gastrula, ";"))
Gastrula.df = data.frame(apply(Gastrula.df, 2, as.numeric))
names(Gastrula.df) = paste("Gastrula.", 1:3, sep = "")
head(Gastrula.df)

Prism.df = do.call("rbind", strsplit(timecourse_file6$Prism, ";"))
Prism.df = data.frame(apply(Prism.df, 2, as.numeric))
names(Prism.df) = paste("Prism.", 1:3, sep = "")
head(Prism.df)


Pluteus.df = do.call("rbind", strsplit(timecourse_file6$Pluteus, ";"))
Pluteus.df = data.frame(apply(Pluteus.df, 2, as.numeric))
names(Pluteus.df) = paste("Pluteus.", 1:3, sep = "")
head(Pluteus.df)

#cbind since same dimensions are used
timecourse_full <- cbind(timecourse_file6, Egg.df, twocell.df, cell16.df, Morula.df, Blastula.df, Gastrula.df, Prism.df, Pluteus.df)

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

timecourse_full <- completeFun(timecourse_full, c("Egg.1", "Egg.2", "Egg.3",
                                                  "2cell.1", "2cell.2", "2cell.3",
                                                  "16cell.1", "16cell.2", "16cell.3",
                                                  "Morula.1", "Morula.2", "Morula.3",
                                                  "Blastula.1", "Blastula.2", "Blastula.3",
                                                  "Gastrula.1", "Gastrula.2", "Gastrula.3",
                                                  "Prism.1", "Prism.2", "Prism.3",
                                                  "Pluteus.1",  "Pluteus.2",  "Pluteus.3"))

########PCA ##########
#my_pca_function <- function(indat = timecourse_ful) 
#{
 # timecourse_samples <- t(indat[,c("Egg.1", "Egg.2", "Egg.3",
  #                                           "2cell.1", "2cell.2", "2cell.3",
   #                                          "16cell.1", "16cell.2", "16cell.3",
    #                                         "Morula.1", "Morula.2", "Morula.3",
     #                                        "Blastula.1", "Blastula.2", "Blastula.3",
      #                                       "Gastrula.1", "Gastrula.2", "Gastrula.3",
       #                                      "Prism.1", "Prism.2", "Prism.3",
        #                                     "Pluteus.1",  "Pluteus.2",  "Pluteus.3")])
  
#}
timecourse_samples <- t(timecourse_full[,c("Egg.1", "Egg.2", "Egg.3",
                                           "2cell.1", "2cell.2", "2cell.3",
                                           "16cell.1", "16cell.2", "16cell.3",
                                           "Morula.1", "Morula.2", "Morula.3",
                                           "Blastula.1", "Blastula.2", "Blastula.3",
                                           "Gastrula.1", "Gastrula.2", "Gastrula.3",
                                           "Prism.1", "Prism.2", "Prism.3",
                                           "Pluteus.1",  "Pluteus.2",  "Pluteus.3", "Protein.ID")])
timecourse_samples <- as.data.frame(timecourse_samples)
colnames(timecourse_samples) <- timecourse_samples[25,]
timecourse_samples <- timecourse_samples[-25,] 

sapply(1:length(timecourse_samples), function(i) {
  timecourse_samples[, i] <<- as.numeric(as.character(timecourse_samples[, i]))
});


######Finding top variable proteins ()#######
#For convenience timecourse_samples was written as row x column (24 x lots of proteins)
#To get most variable proteins (find important rows from the transposed matrix)
rv <- rowVars(as.matrix(t(timecourse_samples)))
hist(rv)
select = order(rv, decreasing=TRUE)[seq_len(250)]
#most_variable_proteins <- timecourse_samples[ ,select]
pca_mv = prcomp(timecourse_samples[ ,select], center = TRUE, scale = TRUE)
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




##########EVerything below here uses all proteins###########
#######using everything (all proteins)######
group_timecourse <- factor(gsub("\\..*","",rownames(timecourse_samples)))
pca_timecourse <- prcomp(timecourse_samples, center = TRUE, scale = TRUE)
#nice PCA plot with control
PCA <- as.data.frame(as.matrix(pca_timecourse$x))
PCA$stages <- rownames(PCA)
ggplot(PCA,aes(x=PC1,y=PC2,col=stages))+geom_point()


PCA_proteins_loadings <- as.data.frame(as.matrix(pca_timecourse$rotation))
PCA_proteins_loadings$Protein.ID <- rownames(PCA_proteins_loadings)
PCA_proteins_loadings$Protein.ID<- gsub("\\..*","",PCA_proteins_loadings$Protein.ID)














#####Looking at Distribution of loadings####


hist(PCA_proteins_loadings$PC1)
hist(PCA_proteins_loadings$PC2)
hist(PCA_proteins_loadings$PC3)







######Analysis of PCA loadings########
#Find proteins with high loading component score to PC1
#Find proteins with high loading component score to PC2

#Annotate the data
#Conversion of IDs
NCBItoGENEID <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/NCBItoGENEPAGE.txt", header = F)
GENEIDtoGENEPAGE <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/GENEIDtoGENEPAGE.txt", header = F)
NCBItoGENEID <- subset(NCBItoGENEID, select = c("V2","V3"))
GENEIDtoGENEPAGE <- subset(GENEIDtoGENEPAGE, select = c("V1","V3"))

timecourse_full$Protein.ID<- gsub("\\..*","",timecourse_full$Protein.ID)


colnames(NCBItoGENEID) <- c("Protein.ID","GENEID")
colnames(GENEIDtoGENEPAGE) <- c("GENEPAGE","GENEID")
timecourse_full <- merge(timecourse_full, NCBItoGENEID, by = "Protein.ID" )
timecourse_full <- merge(timecourse_full, GENEIDtoGENEPAGE, by = "GENEID")
timecourse_full <- timecourse_full[!duplicated(timecourse_full$GENEPAGE), ]


###Merge the two files with loading and full timecourse###
timecourse_withloading <- merge(timecourse_full, PCA_proteins_loadings, by = "Protein.ID")

###Pathway input files + Small Conversions###
library(KEGGREST)
library(clusterProfiler)
library(enrichplot)
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


####Make a function
ORA_analysis <- function(genelist,FUNCt2g, FUNCver2_t2g, GOt2g, GOt2n,KEGGt2g, KEGGt2n){
  
  #You can change the pvalueCutoff, and pAdjustMethod
  funcver1 = enricher(genelist, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = FUNCt2g)
  funcver2 = enricher(genelist, pvalueCutoff = 1,  TERM2GENE = FUNCver2_t2g)
  GO =  enricher(genelist, pvalueCutoff = 1,  TERM2GENE = GOt2g, TERM2NAME = GOt2n)
  KEGG =  enricher(genelist, pvalueCutoff = 1,  TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n)
  list_ORA <- list(funcver1, funcver2, GO, KEGG)
  return(list_ORA)
}


#Timecourse_withloading

PC1_order <- timecourse_withloading[order(timecourse_withloading$PC1, decreasing = TRUE),]
plot(PC1_order)
top_genesPC1 <- PC1_order[1:50,]$GENEPAGE
pc1 <- ORA_analysis(top_genesPC1,FUNCt2g, FUNCver2_t2g, GOt2g, GOt2n,KEGGt2g, KEGGt2n)
dotplot(pc1[[1]])+ ggtitle("FUNC")
dotplot(pc1[[2]])+ ggtitle("FUNC2")
dotplot(pc1[[3]])+ ggtitle("GO")
dotplot(pc1[[4]])+ ggtitle("KEGG")



PC2_order <- timecourse_withloading[order(timecourse_withloading$PC2, decreasing = TRUE),]
top_genesPC2 <- PC2_order[1:50,]$GENEPAGE
pc2 <- ORA_analysis(top_genesPC2,FUNCt2g, FUNCver2_t2g, GOt2g, GOt2n,KEGGt2g, KEGGt2n)
dotplot(pc2[[1]])+ ggtitle("FUNC")
dotplot(pc2[[2]])+ ggtitle("FUNC2")
dotplot(pc2[[3]])+ ggtitle("GO")
dotplot(pc2[[4]])+ ggtitle("KEGG")


#GSEA?

pc1 <- timecourse_withloading$PC1
names(pc1) <- timecourse_withloading$GENEPAGE
pc1[pc1 == Inf] <- 0
pc1 <- sort(pc1, decreasing = TRUE)



pc2 <- timecourse_withloading$PC2
names(pc2) <- timecourse_withloading$GENEPAGE
pc2[pc2 == Inf] <- 0
pc2 <- sort(pc2, decreasing = TRUE)


pc1.KEGG <- GSEA(pc1, pvalueCutoff = 0.5,  TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, verbose=FALSE)
ridgeplot(pc1.KEGG)

pc2.KEGG <- GSEA(pc2, pvalueCutoff = 0.5,  TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, verbose=FALSE)
ridgeplot(pc2.KEGG)





fviz_eig(pca_timecourse) + ggtitle("Timecourse scree plot")
autoplot(pca_timecourse, data = timecourse_samples, 
         label = TRUE, label.size = 3)
autoplot(pca_timecourse, data = timecourse_samples, 
         label = TRUE, label.size = 3, loadings = TRUE, loadings.colour = 'blue', 
         loadings.label = TRUE, loadings.label.size = 3)
pca3d(pca_timecourse, group = group_timecourse, show.labels = T)


#Add PC components. 