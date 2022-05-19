#####Libraries######
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
library(marray)
#for PCA plotting
library(ggforce)
library(scales)
#heatmap creation
library(ComplexHeatmap)
####Annotations File Maker + Preprocessing##########
#directory for annotations
#C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Full Annotations/

#More Input files
###Pathway input files###
FUNCt2g <- read.csv(file = "Data/Anno/FUNCt2g.csv", row.names = 1)
FUNCt2n <- read.csv(file = "Data/Anno/FUNCt2n.csv", row.names = 1)
GOt2g <- read.csv(file = "Data/Anno/GOt2g.csv", row.names = 1)
GOt2n <- read.csv(file = "Data/Anno/GOt2n.csv", row.names = 1)
KEGGt2g <- read.csv(file = "Data/Anno/KEGGt2g.csv", row.names = 1)
KEGGt2n <- read.csv(file = "Data/Anno/KEGGt2n.csv", row.names = 1)
FUNCver2_t2g <- read.csv(file = "Data/Anno/FUNCver2_t2g.csv", row.names = 1)


#directory for data (Should be similar to app)
timecourse_file6 <- read.delim(file = "Data/prot_data/timecourse_file6.txt")
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


#strsplit the time columns and re-add these columns back into the main file
timepoints <- c("Egg", "2cell", "16cell", "Morula", "Blastula", "Gastrula", "Prism", "Pluteus")
for (i in timepoints){
  time <- as.data.frame(apply(do.call("rbind", strsplit(timecourse_file6[[i]], ";")), 2, as.numeric))
  names(time) <- paste0(i,".",1:3,sep = "")
  timecourse_file6 <- cbind(timecourse_file6, time) #cbind since same dimensions are used
}

timecourse_full <- timecourse_file6



########Area to choose method (restrict to only samples that give values)####
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

#####EXTRACT CSV FILE FROM HERE For app
#write.csv(timecourse_full, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/TIMECOURSE/gsea.timecourse.proteomics.csv")



#########PCA Results (Supplemental Figure)#########
timecourse_samples <- t(timecourse_full[,c("Egg.1", "Egg.2", "Egg.3",
                                           "2cell.1", "2cell.2", "2cell.3",
                                           "16cell.1", "16cell.2", "16cell.3",
                                           "Morula.1", "Morula.2", "Morula.3",
                                           "Blastula.1", "Blastula.2", "Blastula.3",
                                           "Gastrula.1", "Gastrula.2", "Gastrula.3",
                                           "Prism.1", "Prism.2", "Prism.3",
                                           "Pluteus.1",  "Pluteus.2",  "Pluteus.3", "Protein.ID")])
timecourse_samples <- as.data.frame(timecourse_samples)
colnames(timecourse_samples) <- timecourse_samples[rownames(timecourse_samples) == "Protein.ID",]
timecourse_samples <- timecourse_samples[1:24,] 


sapply(1:length(timecourse_samples), function(i) {
  timecourse_samples[, i] <<- as.numeric(as.character(timecourse_samples[, i]))
});

######PCA on every protein found in the dataset#####
pca_tot <- prcomp(timecourse_samples, center = TRUE, scale = TRUE)
#PCA plot
pca.plot <-  as.data.frame(as.matrix(pca_tot$x))
pca.plot$stages <- rownames(pca.plot)
pca.plot$cell_stages <- gsub("\\..*","",pca.plot$stages)
x <- summary(pca_tot)
per_var_exp <- as.data.frame(x$importance)

#Order legend correctly (specify factor levels)
pca.plot$cell_stages <- factor(pca.plot$cell_stages, levels = c("Egg",
                                                                "2cell",
                                                                "16cell",
                                                                "Morula",
                                                                "Blastula",
                                                                "Gastrula",
                                                                "Prism",
                                                                "Pluteus"))

library(scales)
timecourse_PCA <- ggplot(pca.plot,aes(x=PC1,y=PC2,col=cell_stages, label = stages))+
  geom_mark_ellipse(aes(x=PC1,y=PC2,fill = cell_stages), show.legend = FALSE, inherit.aes = FALSE) +
  geom_point()+ 
  geom_text(aes(label = stages), color = "black", hjust=1.1, vjust=1.1, size = 4.5) +
  xlab(paste0("PC1 (",percent(per_var_exp[2,c("PC1")], accuracy = 0.1), ")")) +
  ylab(paste0("PC2 (",percent(per_var_exp[2,c("PC2")], accuracy = 0.1), ")")) +
  ggtitle("PCA of Proteomics Timecourse") +
  theme_light() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 14)
        )  
timecourse_PCA
ggsave(filename = "timecourse_PCA.png", plot = timecourse_PCA,
       path = "Output/timecourse_supplementary_figures",
       width = 14,
       height = 11,
       units = "in",
       dpi = 600)

###########Heatmap of significant proteins (Supplemental Figure: AVERAGE version)############
timecourse_ht_sign_avg <- subset(timecourse_full, select = c("Protein.ID","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                                   "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus", "Multiple"))
sapply(2:10, function(i) {
  timecourse_ht_sign_avg[, i] <<- as.numeric(as.character(timecourse_ht_sign_avg[, i]))
});

timecourse_ht_sign_avg <- subset(timecourse_ht_sign_avg, timecourse_ht_sign_avg$Multiple < 0.05)
timecourse_ht_sign_avg$Multiple <- NULL
timecourse_ht_sign_avg$Protein.ID <- NULL



#remove everything before the period (take out Avg signpost)
timecourse_avg_names <- gsub("^.*\\.","",colnames(timecourse_ht_sign_avg))

level = factor(timecourse_avg_names, levels = c("Egg",
                                            "2cell",
                                            "16cell",
                                            "Morula",
                                            "Blastula",
                                            "Gastrula",
                                            "Prism",
                                            "Pluteus"))

#IDS in the color pallete of ggplot2
#library(scales)
#show_col(hue_pal()(8))
colors_gg <- list('timepoint' = c("Egg" = "#F8766D",
                                  "2cell" = "#CD9600",
                                  "16cell" = "#7CAE00",
                                  "Morula" = "#00BE67",
                                  "Blastula" = "#00BFC4",
                                  "Gastrula" = "#00A9FF",
                                  "Prism" = "#C77CFF",
                                  "Pluteus" = "#FFB1CC"))

ha = HeatmapAnnotation(timepoint = level, annotation_name_side = "left", col = colors_gg)


#####No Z-score version
ht_list = Heatmap(as.matrix(timecourse_ht_sign_avg), #from complexHeatmap github
                  name = "scale",
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list)

#####Z-score version
png("Output/timecourse_supplementary_figures/timecourse_heatmap.png",width=6.25,height=8.25,units="in",res=1200)
ht_list_avg_timecourse = Heatmap(t(scale(t(as.matrix(timecourse_ht_sign_avg)))), #from complexHeatmap github
                  name = "scale",
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list_avg_timecourse)
dev.off()

#######Mfuzz clustering (Supplemental Figure)########
######With significant genes, z-score standardize
timecourse_significant_mfuzz <- subset(timecourse_full, select = c("Protein.ID","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                                   "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus", "Multiple"))
sapply(2:10, function(i) {
  timecourse_significant_mfuzz[, i] <<- as.numeric(as.character(timecourse_significant_mfuzz[, i]))
});

timecourse_significant_mfuzz <- subset(timecourse_significant_mfuzz, timecourse_significant_mfuzz$Multiple < 0.05)
timecourse_significant_mfuzz$Multiple <- NULL


####PATHWAY MAPS Annotation
annotation_matrix <- read.delim(file = "Data/Anno/annotation_matrix.txt")
genename_to_ENTREZID <- read.delim(file = "Data/Anno/GeneExternalRef.txt", header = F)
colnames(genename_to_ENTREZID) <- c("GENE.ID","LOC","Descrp","org.code","ENTREZ")
genename_to_ENTREZID <- subset(genename_to_ENTREZID, select = c("GENE.ID","ENTREZ","Descrp"))
full_annotation_matrix <- merge(annotation_matrix, genename_to_ENTREZID, by = "GENE.ID")

timecourse_significant_mfuzz$Protein.ID<- gsub("\\..*","",timecourse_significant_mfuzz$Protein.ID)
names(timecourse_significant_mfuzz)[names(timecourse_significant_mfuzz) == "Protein.ID"] <- "NCBI.ID"


timecourse_significant_mfuzz <- merge(timecourse_significant_mfuzz, full_annotation_matrix[,c("NCBI.ID","GENEPAGE.ID")], by = "NCBI.ID")


#Take average of duplicated data
timecourse_significant_mfuzz$NCBI.ID <- NULL
timecourse_significant_mfuzz <- aggregate(.~GENEPAGE.ID, data = timecourse_significant_mfuzz, FUN = mean)

#Simplify dataframe
rownames(timecourse_significant_mfuzz) <- timecourse_significant_mfuzz$GENEPAGE.ID
timecourse_significant_mfuzz$GENEPAGE.ID <- NULL



sapply(1:8, function(i) {
  timecourse_significant_mfuzz[, i] <<- as.numeric(as.character(timecourse_significant_mfuzz[, i]))
});

timepoint <- c(0,2,4.5,10,24,48,72,96)
# bind that to the dataframe
timecourse_significant_mfuzz <- rbind(timepoint, timecourse_significant_mfuzz)
rownames(timecourse_significant_mfuzz)[1]<-"time"


#data <- table2eset(timecourse_signficant_mfuzz)
tmp <- tempfile()
write.table(timecourse_significant_mfuzz,file=tmp, sep='\t', quote = F, col.names=NA)
data <- table2eset(file=tmp)

data <- standardise(data)
m1 <- mestimate(data)
#Dmin(data, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)
clust = 6
c <- mfuzz(data,c=clust,m=m1, iter.max = 500, verbose = T)


#Errors for Convergence either: 1.2458 and 1.2183 (Tests done below)
errors <- list()
for(i in 1:100){
  c <- mfuzz(data,c=clust,m=m1, iter.max = 500)
  errors[[i]] <- c$withinerror
}
#Plot of errors after running 100 times
hist(unlist(errors[1:100])) 


#Simple Test to check if we are proceeding with calculations only related to global minimum
if(c$withinerror > 1.2184){
  stop("You detected a local minimum, retry the calculation")
}

#Two important plots. 1.) Clustering 2.) Membership score legend
png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_clustering.png",width=12.25,height=8.25,units="in",res=1200)
mfuzz.plot2(data,cl=c, mfrow=c(3,2),
            time.labels= #c(0,2,4.5,10,24,48,72,96)
              c("Egg","2cell","16cell","Morula","Blastula","Gastrula","Prism","Pluteus"),
            cex.main=1.8,cex.lab=1.6, cex.axis = 1.4 ,x11 =  FALSE)
dev.off()

mfuzzColorBar()
mfrow = c(1,1)
mfuzzColorBar(main="Membership",cex.main=1.4, cex.axis = 3, cex.lab = 2)


cor(t(c[[1]]))

#extracts membership values 
acore <- acore(data,c,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


mfuzz <- list()
for (i in 1:clust) {
  names <- paste0("Cluster",i)
  mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= 0.7,]$NAME
}


#Previous p-value cutoff 0.2
KEGG_clust <- compareCluster(geneCluster = mfuzz,  fun = "enricher", TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, pvalueCutoff = 0.2)
png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_clustering_KEGG.png",width=12.25,height=8.25,units="in",res=1200)
dotplot(KEGG_clust, showCategory = 10) + ggtitle("KEGG Pathway Enrichment in each Cluster")
dev.off()

GO_clust <- compareCluster(geneCluster = mfuzz,  fun = "enricher", TERM2GENE = GOt2g, TERM2NAME = GOt2n, pvalueCutoff = 0.2)
#method to remove NA
#mutate result with complete cases argument
GO_clust@compareClusterResult <- GO_clust@compareClusterResult[complete.cases(GO_clust@compareClusterResult), ]
png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_clustering_GO.png",width=12.25,height=8.25,units="in",res=1200)
dotplot(GO_clust, showCategory = 10) + ggtitle("GO Pathway Enrichment in each Cluster")
dev.off()


FUNC_clust <- compareCluster(geneCluster = mfuzz,  fun = "enricher", TERM2GENE = FUNCt2g, pvalueCutoff = 0.2)
png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_clustering_FUNC.png",width=12.25,height=8.25,units="in",res=1200)
dotplot(FUNC_clust, showCategory = 10) + ggtitle("FUNC Pathway Enrichment in each Cluster")
dev.off()


#####csv files
KEGG_clusters <- KEGG_clust@compareClusterResult
GO_clusters <- GO_clust@compareClusterResult
FUNC_clusters <- FUNC_clust@compareClusterResult

#######Mfuzz clustering only metabolic enzymes (Main Figure)#######
timecourse_significant_mfuzz_metab <- subset(timecourse_full, select = c("Protein.ID","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                                   "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus", "Multiple"))
annotation_matrix <- read.delim(file = "Data/Anno/annotation_matrix.txt")
##data should have columns {Protein.ID, Avg.Egg, Avg.2cell, Avg.16cell, Avg.Morula, Avg.Blastula, Avg.Gastrula, Avg.Prism, Avg.Pluteus, Multiple}
#filtering = "spu01100" #metabolism cat.
#filtering = grepl("spu01100", data_full$KEGG.Pathway.ID)
mfuzz_filtering <- function(data, filter1, filter2, filter3, pvalcutoff){
  data$Protein.ID <- gsub("\\..*","",data$Protein.ID)
  #Add annotations based on protein ID
  annotation_matrix <- read.delim(file = "Data/Anno/annotation_matrix.txt")
  ##Merge the timecourse + protein.to.entrez
  data_full <- merge(data, annotation_matrix, by.x = "Protein.ID", by.y = "NCBI.ID")
  
  ##Filtering by pathway
  #old version: grepl(filter1, data_full[[filter2]]) 
  tc_data_full <- subset(data_full, str_detect(data_full[[filter2]], filter1, negate = filter3))
  
  ##Char. to numeric
  sapply(2:10, function(i) {
    tc_data_full[, i] <<- as.numeric(as.character(tc_data_full[, i]))
  });
  
  ##Filtering by p-value (significance)
  tc_data_full <- subset(tc_data_full, tc_data_full$Multiple < pvalcutoff)
  tc_data_full$Multiple <- NULL
  
  ##Set up for mfuzz plots
  tc_data_full <- subset(tc_data_full, select = c("Protein.ID","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                  "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus", "GENEPAGE.ID"))
  
  #Take average of duplicated data
  tc_data_full$Protein.ID <- NULL
  tc_data_full <- aggregate(.~GENEPAGE.ID, data = tc_data_full, FUN = mean)
  
  #Simplify dataframe
  rownames(tc_data_full) <- tc_data_full$GENEPAGE.ID
  tc_data_full$GENEPAGE.ID <- NULL
  sapply(1:8, function(i) {
    tc_data_full[, i] <<- as.numeric(as.character(tc_data_full[, i]))
  });
  
  
  timepoint <- c(0,2,4.5,10,24,48,72,96)
  # bind that to the dataframe
  tc_data_full <- rbind(timepoint, tc_data_full)
  rownames(tc_data_full)[1]<-"time"
  return(tc_data_full)
}
mfuzz_plot <- function(mfuzz_data,clustnum,memscore){
  tmp <- tempfile()
  write.table(mfuzz_data,file=tmp, sep='\t', quote = F, col.names=NA)
  data <- table2eset(file=tmp)
  
  data <- standardise(data)
  m1 <- mestimate(data)
  #Dmin(data, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)
  clust = clustnum
  c <- mfuzz(data,c=clust,m=m1, iter.max = 500, verbose = T)
  
  #Errors for Convergence either: 1.2458 and 1.2183 (Tests done below)
  errors <- list()
  for(i in 1:100){
    c <- mfuzz(data,c=clust,m=m1, iter.max = 500)
    errors[[i]] <- c$withinerror
  }
  #Plot of errors after running 100 times
  hist(unlist(errors[1:100])) 
  
  #Simple Test to check if we are proceeding with calculations only related to global minimum
  if(c$withinerror > 1.170){
    stop("You detected a local minimum, retry the calculation")
  }
  
  Dmin(data, m=m1, crange=seq(2,22,1), repeats=10, visu=TRUE)
  
  png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_enzyme_only.png",width=12.25,height=8.25,units="in",res=1200)
  mfuzz.plot2(data,cl=c, mfrow=c(3,2),
              time.labels= c("Egg","2cell","16cell","Morula","Blastula","Gastrula","Prism","Pluteus"),
              cex.main=1.8,cex.lab=1.6, cex.axis = 1.4 ,x11 =  FALSE)
  dev.off()
  
  mfuzzColorBar(main="Membership",cex.main=1.4, cex.axis = 3, cex.lab = 2)
  
  #extracts membership values 
  acore <- acore(data,c,min.acore=0)
  acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))
  
  
  mfuzz <- list()
  for (i in 1:clust) {
    names <- paste0("Cluster",i)
    mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= memscore,]$NAME
  }
  
  return(mfuzz)
}
mfuzz_clusters <- function(mfuzz_clusts, pvalcutoff, title,showcatparam){
  pathway_plots <- NULL
  KEGG_clust <- compareCluster(geneCluster = mfuzz_clusts,  fun = "enricher", TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, pvalueCutoff = pvalcutoff)
  pathway_plots[["KEGG"]] <- dotplot(KEGG_clust, showCategory = showcatparam) + ggtitle(paste0("KEGG Pathway Enrichment in each Cluster:", title))
  
  GO_clust <- compareCluster(geneCluster = mfuzz_clusts,  fun = "enricher", TERM2GENE = GOt2g, TERM2NAME = GOt2n, pvalueCutoff = pvalcutoff)
  pathway_plots[["GO"]]<- dotplot(GO_clust, showCategory = showcatparam) + ggtitle(paste0("GO Pathway Enrichment in each Cluster:",title))
  
  FUNC_clust <- compareCluster(geneCluster = mfuzz_clusts,  fun = "enricher", TERM2GENE = FUNCt2g, pvalueCutoff = pvalcutoff)
  pathway_plots[["FUNC"]] <- dotplot(FUNC_clust, showCategory = showcatparam) + ggtitle(paste0("FUNC Pathway Enrichment in each Cluster:",title))
  
  return(pathway_plots)
}

#filter1 = "spu01100"
#filter2 = "KEGG.Pathway.ID"
#filter3 = FALSE
#pvalcutoff = 0.05
#data = timecourse_significant_mfuzz_metab
tc_sign_mfuzz_metab_fin <- mfuzz_filtering(timecourse_significant_mfuzz_metab, "spu01100", "KEGG.Pathway.ID", FALSE, 0.05)


kegg_enzyme_lists <- mfuzz_plot(tc_sign_mfuzz_metab_fin, 4, 0)
enzyme_paths_no_manual <- mfuzz_clusters(kegg_enzyme_lists, 0.2, "KEGG ENZYMES", 10)
png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_enzyme_only_KEGG.png",width=12.25,height=8.25,units="in",res=1200)
plot(enzyme_paths_no_manual$KEGG)
dev.off()

png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_enzyme_only_GO.png",width=12.25,height=8.25,units="in",res=1200)
plot(enzyme_paths_no_manual$GO)
dev.off()

png("Output/timecourse_supplementary_figures/timecourse_fuzzyc_enzyme_only_FUNC.png",width=12.25,height=8.25,units="in",res=1200)
plot(enzyme_paths_no_manual$FUNC)
dev.off()




#######Manual edit: Since Cluster #'s change I manually editted the ordering to produce ClusterProfiler Plots (Main FIgure)#######
# #Have to edit below each time
# names(kegg_enzyme_lists) <- c("Blastula", "Pluteus", "Gastrula","Morula")
# #Edit the list to get "Morula -> Blastula -> Gastrula -> Pluteus"
# vector_order = c(4,1,3,2)
# kegg_enzyme_lists <- kegg_enzyme_lists[vector_order]
# ##Manual edit (end)
# 
# enzyme_paths <- mfuzz_clusters(kegg_enzyme_lists,0.2,"KEGG ENZYMES", 10)
# plot(enzyme_paths$KEGG)
# plot(enzyme_paths$GO)
# plot(enzyme_paths$FUNC)
# 
# 
# tc_sign_mfuzz_metab_fin <- mfuzz_filtering(timecourse_significant_mfuzz_metab, "spu01100", "KEGG.Pathway.ID", FALSE, 0.1)
# kegg_enzyme_lists <- mfuzz_plot(tc_sign_mfuzz_metab_fin, 6, 0)
# enzyme_paths <- mfuzz_clusters(kegg_enzyme_lists,0.1,"KEGG ENZYMES", 10)
# plot(enzyme_paths$KEGG)
# plot(enzyme_paths$GO)
# plot(enzyme_paths$FUNC)