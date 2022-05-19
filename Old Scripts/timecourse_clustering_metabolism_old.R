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
library(ComplexHeatmap)
library(dplyr)
library(HybridMTest)
#Set seed.

######Function Creation#########
#Make function for cleaning the files
####Function #1
metabolite_clean1 <- function(metabolite_list){
  metabolite_list <- gsub(".TMS","",metabolite_list)
  metabolite_list <- gsub(".2TMS","",metabolite_list)
  metabolite_list <- gsub(".3TMS","",metabolite_list)
  metabolite_list <- gsub(".4TMS","",metabolite_list)
  metabolite_list <- gsub(".5TMS","",metabolite_list)
  metabolite_list <- gsub(".6TMS","",metabolite_list)
  metabolite_list <- gsub("X1.","1-",metabolite_list)
  metabolite_list <- gsub("X2.","2-",metabolite_list)
  metabolite_list <- gsub("X3.","3-",metabolite_list)
  metabolite_list <- gsub("X4.","4-",metabolite_list)
  metabolite_list <- gsub("X5.","5-",metabolite_list)
  metabolite_list <- gsub(".meto","",metabolite_list)
  metabolite_list <- gsub("..1","",metabolite_list)
  metabolite_list <- gsub("..2","",metabolite_list)
  metabolite_list <- gsub(".d3.","",metabolite_list)
  metabolite_list <- gsub("."," ", metabolite_list, fixed=TRUE)
  metabolite_list <- trimws(metabolite_list)
  return(metabolite_list)}

###Clean the KEGG IDs from Generic sugar to alpha type...
metabolite_alpha <- function(pathview_set){
  #for glucose-6-phosphate (alpha version)
  rownames(pathview_set)[rownames(pathview_set) == "C00092"] <- "C00668"
  
  #for beta-D-glucose -> alpha-D-glucose
  rownames(pathview_set)[rownames(pathview_set) == "C00221"] <- "C00267"
  return(pathview_set)
}

#Confirmed want Paired T-test for these comparisons
#Need to make global metabolite name matcher:
####New PreProcess####
######Timecourse######
#Using the normalized, quantile data {the heatmapped version}
metab_timecourse <- read.csv(file = "Data/metab_data/timeseries_normalized_quantile.csv")
metab_timecourse <- as.data.frame(t(metab_timecourse))
colnames(metab_timecourse) <- metab_timecourse[1,]
metab_timecourse<- metab_timecourse[-1,]
metab_timecourse<- metab_timecourse[-1,]
metab_timecourse$Query <- metabolite_clean1(rownames(metab_timecourse))
write.csv(metab_timecourse$Query, file = "Output/metab_global/Files/metabname_timecourse.csv")
KEGGMOLID_5 <- read.csv(file = "Data/metab_conv/name_map_timecourse.csv")

metab_timecourse_full <- merge(metab_timecourse, KEGGMOLID_5, by = "Query")
metab_timecourse_full<- subset(metab_timecourse_full, select = c("Egg_A","Egg_B","Egg_C",
                                                                 "2 cell_A","2 cell_B","2 cell_C",
                                                                 "16 cell_A","16 cell_B","16 cell_C",
                                                                 "Morula_A","Morula_B","Morula_C",
                                                                 "Blastula_A","Blastula_B","Blastula_C",
                                                                 "Gastrula_A","Gastrula_B","Gastrula_C",
                                                                 "Prism_A","Prism_B","Prism_C",
                                                                 "Pluteus_A","Pluteus_B","Pluteus_C",
                                                                 "KEGG"))

sapply(1:24, function(i) {
  metab_timecourse_full[, i] <<- as.numeric(as.character(metab_timecourse_full[, i]))
});

###ANOVA Test
anova_res <- row.oneway.anova(metab_timecourse_full[,1:24],c(1,1,1,2,2,2,3,3,3,
                                                   4,4,4,5,5,5,6,6,6,
                                                   7,7,7,8,8,8))
anova_res$p.adjust <- p.adjust(anova_res$pval, method = "BH", n = length(anova_res$pval))
metab_timecourse_full$p.adjust <- anova_res$p.adjust






####Aggregate by average 
#Take mean of duplicate ID's
metab_timecourse_full <- aggregate(.~KEGG, data = metab_timecourse_full, FUN = mean, na.action = na.pass)




###Get average of timepoints
col <- colnames(metab_timecourse_full)[2:25]
metab_timecourse_full$Avg.Egg <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Egg_A", "Egg_B", "Egg_C")]))
metab_timecourse_full$Avg.2cell <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("2 cell_A", "2 cell_B", "2 cell_C")]))
metab_timecourse_full$Avg.16cell <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("16 cell_A", "16 cell_B", "16 cell_C")]))
metab_timecourse_full$Avg.Morula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Morula_A", "Morula_B", "Morula_C")]))
metab_timecourse_full$Avg.Blastula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Blastula_A", "Blastula_B", "Blastula_C")]))
metab_timecourse_full$Avg.Gastrula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Gastrula_A", "Gastrula_B", "Gastrula_C")]))
metab_timecourse_full$Avg.Prism <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Prism_A", "Prism_B", "Prism_C")]))
metab_timecourse_full$Avg.Pluteus <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Pluteus_A", "Pluteus_B", "Pluteus_C")]))




#Get timecourse for mfuzz clust:
metab_timecourse_full_pathview <- subset(metab_timecourse_full, select = c("KEGG","Avg.Egg","Avg.2cell","Avg.16cell","Avg.Morula","Avg.Blastula",
                                                                           "Avg.Gastrula","Avg.Prism","Avg.Pluteus","p.adjust"))
rownames(metab_timecourse_full_pathview) <- metab_timecourse_full_pathview$KEGG
metab_timecourse_full_pathview$KEGG <- NULL

###Modified ID changes
metab_timecourse_full_pathview <- metabolite_alpha(metab_timecourse_full_pathview)


####NEw MFuzz code####
metab_timecourse_mfuzz_sign <- metab_timecourse_full_pathview
metab_timecourse_mfuzz_sign <- metab_timecourse_mfuzz_sign[order(metab_timecourse_mfuzz_sign$p.adjust),]
metab_timecourse_mfuzz_sign <- metab_timecourse_mfuzz_sign[1:50,] #top 50 metabolites
metab_timecourse_mfuzz_sign$p.adjust <- NULL



timepoint <- c(0,2,4.5,10,24,48,72,96)
# bind that to the dataframe
metab_timecourse_mfuzz_sign <- rbind(timepoint,metab_timecourse_mfuzz_sign)
rownames(metab_timecourse_mfuzz_sign)[1]<-"time"



#######Mfuzz - zscore
#data <- table2eset(timecourse_signficant_mfuzz)
tmp <- tempfile()
write.table(metab_timecourse_mfuzz_sign,file=tmp, sep='\t', quote = F, col.names=NA)
data <- table2eset(file=tmp)


#Standardize to mean = 0 = zscore
data <- standardise(data)
m1 <- mestimate(data)
Dmin(data, m=m1, crange=seq(2,22,1), repeats=10, visu=TRUE)

#"Try setting the fuzziness parameter m to 1, this will give you the equivalent of k-means clustering"
#m1 = 1
clust = 5
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
if(c$withinerror > 0.744){
  stop("You detected a local minimum, retry the calculation")
}



#Using verbose True to find global minimum.
png("Output/metab_timecourse_figures/timecourse_fuzzyc_metabolism.png",width=12.25,height=8.25,units="in",res=1200)
mfuzz.plot2(data,cl=c,mfrow=c(3,2),time.labels= c("Egg","2cell","16cell","Morula","Blastula","Gastrula","Prism","Pluteus"),
           cex.main=1.8,cex.lab=1.6, cex.axis = 1.4,
           new.window=FALSE,x11 =  FALSE)
dev.off()

cor(t(c[[1]]))



###Membership score + metabolite names:
#extracts membership values 
acore <- acore(data,c,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


mfuzz <- list()
for (i in 1:clust) {
  names <- paste0("Cluster",i)
  mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= 0,]$NAME
  #mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i,]$NAME
}

KEGGt2n <- read.csv(file = "Data/Anno/KEGGt2n.csv", row.names = 1)
KEGGt2m <- read.csv(file = "Data/Anno/KEGG_ENZYME_t2m.csv", row.names = 1)
KEGGt2m <- KEGGt2m[ , c("PathwayID","KEGG")]




####Change ordering manually, mfuzz
#names(mfuzz) <- c("Egg", "Blastula", "Pluteus","Morula.peak","Morula.mat")
#ordering <- c(1,4,5,2,3)
#mfuzz <- mfuzz[ordering] 
###End of manual edit.


KEGG_clust <- compareCluster(geneCluster = mfuzz,  fun = "enricher", TERM2GENE = KEGGt2m, TERM2NAME = KEGGt2n, pvalueCutoff = 0.2)
png("Output/metab_timecourse_figures/timecourse_fuzzyc_metabolism_KEGG.png",width=12.25,height=8.25,units="in",res=1200)
dotplot(KEGG_clust, showCategory = 10) + ggtitle("KEGG Pathway Enrichment in each Cluster (Metabolites)")
dev.off()


########TRASH OLD BELOW#########

####Preprocess####
#Obtain different metab_timecourse with quantile normalization instead
metab_timecourse <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/Timecourse/timeseries_normalized_quantile.csv")
metab_timecourse <- as.data.frame(t(metab_timecourse))
colnames(metab_timecourse) <- metab_timecourse[1,]
metab_timecourse<- metab_timecourse[-1,]
metab_timecourse<- metab_timecourse[-1,]

# metabolite_list <- rownames(metab_timecourse)
# metabolite_list <- gsub(".TMS","",metabolite_list)
# metabolite_list <- gsub(".2TMS","",metabolite_list)
# metabolite_list <- gsub(".3TMS","",metabolite_list)
# metabolite_list <- gsub(".4TMS","",metabolite_list)
# metabolite_list <- gsub(".5TMS","",metabolite_list)
# metabolite_list <- gsub(".6TMS","",metabolite_list)
# metabolite_list <- gsub("X1.","1-",metabolite_list)
# metabolite_list <- gsub("X2.","2-",metabolite_list)
# metabolite_list <- gsub("X3.","3-",metabolite_list)
# metabolite_list <- gsub("X4.","4-",metabolite_list)
# metabolite_list <- gsub("X5.","5-",metabolite_list)
# metabolite_list <- gsub("X7.","7-",metabolite_list)
# metabolite_list <- gsub(".meto","",metabolite_list)
# metabolite_list <- gsub("..1","",metabolite_list)
# metabolite_list <- gsub("..2","",metabolite_list)
# metabolite_list <- gsub(".d3.","",metabolite_list)
# metabolite_list <- gsub("."," ", metabolite_list, fixed=TRUE)
# metabolite_list <- trimws(metabolite_list)
# write.csv(metabolite_list, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/Timecourse/metabnametimecourse.csv")
# 
# 
# 
# KEGGMOLID <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/Timecourse/name_map_timecourse.csv")
# 
# metab_timecourse$Query <- metabolite_list
metab_timecourse$Names <- rownames(metab_timecourse)

sapply(1:24, function(i) {
  metab_timecourse[, i] <<- as.numeric(as.character(metab_timecourse[, i]))
});
#BiocManager::install("HybridMTest")
library(HybridMTest)
lols <- row.oneway.anova(metab_timecourse[,1:24],c(1,1,1,2,2,2,3,3,3,
                                                   4,4,4,5,5,5,6,6,6,
                                                   7,7,7,8,8,8))
lols$p.adjust <- p.adjust(lols$pval, method = "BH", n = length(lols$pval))
metab_timecourse$p.adjust <- lols$p.adjust


# 
# #metab_timecourse_full <- merge(metab_timecourse, KEGGMOLID, by = "Query")
# metab_timecourse_full <- subset(metab_timecourse_full, select = c("Egg_A","Egg_B","Egg_C","2 cell_A","2 cell_B",  
#                                                                   "2 cell_C","16 cell_A","16 cell_B","16 cell_C",
#                                                                   "Morula_A","Morula_B","Morula_C","Blastula_A","Blastula_B",
#                                                                   "Blastula_C","Gastrula_A","Gastrula_B","Gastrula_C",
#                                                                   "Prism_A","Prism_B","Prism_C","Pluteus_A",  
#                                                                   "Pluteus_B","Pluteus_C","KEGG", "Query","Names"))
# 








#Use everything from metab_timecourse (normalized file in google drive)
# sapply(1:24, function(i) {
#   metab_timecourse[, i] <<- as.numeric(as.character(metab_timecourse[, i]))
# });


col <- colnames(metab_timecourse)[1:24]
metab_timecourse$Avg.Egg <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Egg_A", "Egg_B", "Egg_C")]))
metab_timecourse$Avg.2cell <- apply(metab_timecourse[,col],1,function(x) mean(x[c("2 cell_A", "2 cell_B", "2 cell_C")]))
metab_timecourse$Avg.16cell <- apply(metab_timecourse[,col],1,function(x) mean(x[c("16 cell_A", "16 cell_B", "16 cell_C")]))
metab_timecourse$Avg.Morula <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Morula_A", "Morula_B", "Morula_C")]))
metab_timecourse$Avg.Blastula <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Blastula_A", "Blastula_B", "Blastula_C")]))
metab_timecourse$Avg.Gastrula <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Gastrula_A", "Gastrula_B", "Gastrula_C")]))
metab_timecourse$Avg.Prism <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Prism_A", "Prism_B", "Prism_C")]))
metab_timecourse$Avg.Pluteus <- apply(metab_timecourse[,col],1,function(x) mean(x[c("Pluteus_A", "Pluteus_B", "Pluteus_C")]))

######Top 50 metabolites by ANOVA#####
metab_timecourse_mfuzz <- subset(metab_timecourse, select = c("Names","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                              "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus",
                                                              "p.adjust"))
metabolite_list <- metab_timecourse_mfuzz$Names
metabolite_list <- gsub(".TMS","",metabolite_list)
metabolite_list <- gsub(".2TMS","",metabolite_list)
metabolite_list <- gsub(".3TMS","",metabolite_list)
metabolite_list <- gsub(".4TMS","",metabolite_list)
metabolite_list <- gsub(".5TMS","",metabolite_list)
metabolite_list <- gsub(".6TMS","",metabolite_list)
metabolite_list <- gsub("X1.","1-",metabolite_list)
metabolite_list <- gsub("X2.","2-",metabolite_list)
metabolite_list <- gsub("X3.","3-",metabolite_list)
metabolite_list <- gsub("X4.","4-",metabolite_list)
metabolite_list <- gsub("X5.","5-",metabolite_list)
metabolite_list <- gsub("X7.","7-",metabolite_list)
metabolite_list <- gsub(".meto","",metabolite_list)
metabolite_list <- gsub("..1","",metabolite_list)
metabolite_list <- gsub("..2","",metabolite_list)
metabolite_list <- gsub(".d3.","",metabolite_list)
metabolite_list <- gsub("."," ", metabolite_list, fixed=TRUE)
metabolite_list <- trimws(metabolite_list)

metab_timecourse_mfuzz$Alt.Names <- metabolite_list
write.csv(metabolite_list, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/timecourse_.csv")





















metab_timecourse_mfuzz_sign <- metab_timecourse_mfuzz[order(metab_timecourse_mfuzz$p.adjust),]
metab_timecourse_mfuzz_sign <- metab_timecourse_mfuzz_sign[1:50,] #top 50 metabolites
rownames(metab_timecourse_mfuzz_sign) <- metab_timecourse_mfuzz_sign$Names
metab_timecourse_mfuzz_sign$Names <- NULL
metab_timecourse_mfuzz_sign$p.adjust <- NULL

sapply(1:8, function(i) {
  metab_timecourse_mfuzz_sign[, i] <<- as.numeric(as.character(metab_timecourse_mfuzz_sign[, i]))
});

timepoint <- c(0,2,4.5,10,24,48,72,96)
# bind that to the dataframe
metab_timecourse_mfuzz_sign <- rbind(timepoint,metab_timecourse_mfuzz_sign)
rownames(metab_timecourse_mfuzz_sign)[1]<-"time"



#######Mfuzz - zscore
#data <- table2eset(timecourse_signficant_mfuzz)
tmp <- tempfile()
write.table(metab_timecourse_mfuzz_sign,file=tmp, sep='\t', quote = F, col.names=NA)
data <- table2eset(file=tmp)


#Standardize to mean = 0 = zscore
data <- standardise(data)
m1 <- mestimate(data)
Dmin(data, m=m1, crange=seq(2,22,1), repeats=10, visu=TRUE)

#"Try setting the fuzziness parameter m to 1, this will give you the equivalent of k-means clustering"
#m1 = 1
clust = 5
c <- mfuzz(data,c=clust,m=m1)
mfuzz.plot(data,cl=c,mfrow=c(2,2),time.labels= c("Egg","2cell","16cell","Morula","Blastula","Gastrula","Prism","Pluteus"),new.window=FALSE)
cor(t(c[[1]]))



###Membership score + metabolite names:
#extracts membership values 
acore <- acore(data,c,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


mfuzz <- list()
for (i in 1:clust) {
  names <- paste0("Cluster",i)
  mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= 0,]$NAME
  #mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i,]$NAME
}


library(plyr)
clusters <- plyr::ldply(mfuzz, cbind)

metabolite_list <- clusters$"1"
metabolite_list <- gsub(".TMS","",metabolite_list)
metabolite_list <- gsub(".2TMS","",metabolite_list)
metabolite_list <- gsub(".3TMS","",metabolite_list)
metabolite_list <- gsub(".4TMS","",metabolite_list)
metabolite_list <- gsub(".5TMS","",metabolite_list)
metabolite_list <- gsub(".6TMS","",metabolite_list)
metabolite_list <- gsub("X1.","1-",metabolite_list)
metabolite_list <- gsub("X2.","2-",metabolite_list)
metabolite_list <- gsub("X3.","3-",metabolite_list)
metabolite_list <- gsub("X4.","4-",metabolite_list)
metabolite_list <- gsub("X5.","5-",metabolite_list)
metabolite_list <- gsub("X7.","7-",metabolite_list)
metabolite_list <- gsub(".meto","",metabolite_list)
metabolite_list <- gsub("..1","",metabolite_list)
metabolite_list <- gsub("..2","",metabolite_list)
metabolite_list <- gsub(".d3.","",metabolite_list)
metabolite_list <- gsub("."," ", metabolite_list, fixed=TRUE)
metabolite_list <- trimws(metabolite_list)

clusters$name <- metabolite_list
write.csv(clusters, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/cluster_top50_metabtimecourse.csv")



#Input new csv file with KEGG IDs and run ClusterProfiler.
cluster_name_conv <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/name_map_clustertop50_metab_new.csv", header = T)
conversions <- subset(cluster_name_conv, select = c("Query","KEGG"))
clusters_fin <- merge(clusters, conversions, by.x = "name", by.y = "Query")
clusters_fin_path <- unique(subset(clusters_fin, select = c(".id","KEGG")))


#KEGG annotations again.
KEGGt2n <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGGanno_t2n.csv", row.names = 1)
KEGGt2n$PathwayID <- gsub("path:","",KEGGt2n$PathwayID)
KEGGt2m <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Anno/KEGG_ENZYME_t2m.csv", row.names = 1)
KEGGt2m <- KEGGt2m[ , c("PathwayID","KEGG")]





KEGG_clust <- compareCluster(geneCluster = clusters_fin_path,  fun = "enricher", TERM2GENE = KEGGt2m, TERM2NAME = KEGGt2n, pvalueCutoff = 0.2)
dotplot(KEGG_clust, showCategory = 10) + ggtitle("KEGG Pathway Enrichment in each Cluster (Metabolites)")




######All metabolites by ANOVA########

metab_timecourse_mfuzz <- subset(metab_timecourse, select = c("Names","Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                              "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus"
                                                              
))

#Simplify dataframe
rownames(metab_timecourse_mfuzz) <- metab_timecourse_mfuzz$Names
metab_timecourse_mfuzz$Names <- NULL


sapply(1:8, function(i) {
  metab_timecourse_mfuzz[, i] <<- as.numeric(as.character(metab_timecourse_mfuzz[, i]))
});

timepoint <- c(0,2,4.5,10,24,48,72,96)
# bind that to the dataframe
metab_timecourse_mfuzz <- rbind(timepoint,metab_timecourse_mfuzz)
rownames(metab_timecourse_mfuzz)[1]<-"time"


length(which(metab_timecourse$p.adjust < 0.05))





#######Mfuzz - zscore#########
#data <- table2eset(timecourse_signficant_mfuzz)
tmp <- tempfile()
write.table(metab_timecourse_mfuzz,file=tmp, sep='\t', quote = F, col.names=NA)
data <- table2eset(file=tmp)


#Standardize to mean = 0 = zscore
data <- standardise(data)
m1 <- mestimate(data)
Dmin(data, m=m1, crange=seq(2,22,1), repeats=10, visu=TRUE)

#"Try setting the fuzziness parameter m to 1, this will give you the equivalent of k-means clustering"
#m1 = 1
clust = 6
c <- mfuzz(data,c=clust,m=m1)
mfuzz.plot(data,cl=c,mfrow=c(2,2),time.labels=c(0,2,4.5,10,24,48,72,96),new.window=FALSE)
cor(t(c[[1]]))


#extracts membership values 
acore <- acore(data,c,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


mfuzz <- list()
for (i in 1:clust) {
  names <- paste0("Cluster",i)
  mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= 0.6,]$NAME
  #mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i,]$NAME
}


library(plyr)
clusters <- plyr::ldply(mfuzz, cbind)

metabolite_list <- clusters$"1"
metabolite_list <- gsub(".TMS","",metabolite_list)
metabolite_list <- gsub(".2TMS","",metabolite_list)
metabolite_list <- gsub(".3TMS","",metabolite_list)
metabolite_list <- gsub(".4TMS","",metabolite_list)
metabolite_list <- gsub(".5TMS","",metabolite_list)
metabolite_list <- gsub(".6TMS","",metabolite_list)
metabolite_list <- gsub("X1.","1-",metabolite_list)
metabolite_list <- gsub("X2.","2-",metabolite_list)
metabolite_list <- gsub("X3.","3-",metabolite_list)
metabolite_list <- gsub("X4.","4-",metabolite_list)
metabolite_list <- gsub("X5.","5-",metabolite_list)
metabolite_list <- gsub("X7.","7-",metabolite_list)
metabolite_list <- gsub(".meto","",metabolite_list)
metabolite_list <- gsub("..1","",metabolite_list)
metabolite_list <- gsub("..2","",metabolite_list)
metabolite_list <- gsub(".d3.","",metabolite_list)
metabolite_list <- gsub("."," ", metabolite_list, fixed=TRUE)
metabolite_list <- trimws(metabolite_list)

clusters$name <- metabolite_list
write.csv(clusters, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/Timecourse/timecoursemetabolism_clusters.csv")








####Mfuzz - standardize to egg#####

#data <- table2eset(timecourse_signficant_mfuzz)
tmp <- tempfile()
write.table(metab_timecourse_mfuzz,file=tmp, sep='\t', quote = F, col.names=NA)
data <- table2eset(file=tmp)


#Standardize to mean = 0 = zscore
data <- standardise2(data, timepoint = 1)
m1 <- mestimate(data)
Dmin(data, m=m1, crange=seq(2,22,1), repeats=3, visu=TRUE)

#"Try setting the fuzziness parameter m to 1, this will give you the equivalent of k-means clustering"
#m1 = 1
clust = 9
c <- mfuzz(data,c=clust,m=m1)
mfuzz.plot(data,cl=c,mfrow=c(3,3),time.labels=c(0,2,4.5,10,24,48,72,96),new.window=FALSE)
cor(t(c[[1]]))


#extracts membership values 
acore <- acore(data,c,min.acore=0)
acore_list <- do.call(rbind, lapply(seq_along(acore), function(i){ data.frame(CLUSTER=i, acore[[i]])}))


mfuzz <- list()
for (i in 1:clust) {
  names <- paste0("Cluster",i)
  mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i & acore_list$MEM.SHIP >= 0.6,]$NAME
  #mfuzz[[names]] <- acore_list[acore_list$CLUSTER == i,]$NAME
}
































##############Heatmap stuff - for all metabolites########
metab_timecourse_heatmap <- metab_timecourse[,1:25]
#rownames(metab_timecourse_heatmap) <- metab_timecourse_heatmap$Query
metab_timecourse_heatmap$Query <- NULL
rownames(metab_timecourse_heatmap) <- gsub("X1.","1-",rownames(metab_timecourse_heatmap))
rownames(metab_timecourse_heatmap) <- gsub("X2.","2-",rownames(metab_timecourse_heatmap))
rownames(metab_timecourse_heatmap) <- gsub("X3.","3-",rownames(metab_timecourse_heatmap))
rownames(metab_timecourse_heatmap) <- gsub("X4.","4-",rownames(metab_timecourse_heatmap))
rownames(metab_timecourse_heatmap) <- gsub("X5.","5-",rownames(metab_timecourse_heatmap))
rownames(metab_timecourse_heatmap) <- gsub("X7.","7-",rownames(metab_timecourse_heatmap))


timecourse_names <- gsub("\\_.*","",colnames(metab_timecourse_heatmap))
ha = HeatmapAnnotation(timepoint = timecourse_names , annotation_name_side = "left")


#Normalized - no z-score
timecourse_heatmap_metab_all <- as.matrix(metab_timecourse_heatmap)
ht_list = Heatmap(timecourse_heatmap_metab_all,
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list)


#Normalized - no z-score - no column clustering
timecourse_heatmap_metab_all <- as.matrix(metab_timecourse_heatmap)
ht_list = Heatmap(timecourse_heatmap_metab_all,
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE)
draw(ht_list)





#Normalized - with z-score
timecourse_heatmap_metab_all <- as.matrix(metab_timecourse_heatmap)
ht_list = Heatmap(as.matrix(t(scale(t(timecourse_heatmap_metab_all )))),
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list)


#Normalized - with z-score - no column clustering
timecourse_heatmap_metab_all <- as.matrix(metab_timecourse_heatmap)
ht_list = Heatmap(as.matrix(t(scale(t(timecourse_heatmap_metab_all )))),
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE)
draw(ht_list)




######Average of everything

metab_timecourse_heatmap_avg <- subset(metab_timecourse, select = c("Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                                                 "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus"))

rownames(metab_timecourse_heatmap_avg) <- gsub("X1.","1-",rownames(metab_timecourse_heatmap_avg))
rownames(metab_timecourse_heatmap_avg) <- gsub("X2.","2-",rownames(metab_timecourse_heatmap_avg))
rownames(metab_timecourse_heatmap_avg) <- gsub("X3.","3-",rownames(metab_timecourse_heatmap_avg))
rownames(metab_timecourse_heatmap_avg) <- gsub("X4.","4-",rownames(metab_timecourse_heatmap_avg))
rownames(metab_timecourse_heatmap_avg) <- gsub("X5.","5-",rownames(metab_timecourse_heatmap_avg))
rownames(metab_timecourse_heatmap_avg) <- gsub("X7.","7-",rownames(metab_timecourse_heatmap_avg))





timecourse_names <- gsub("\\_.*","",colnames(metab_timecourse_heatmap_avg))
hanno = HeatmapAnnotation(timepoint = timecourse_names , annotation_name_side = "left")

#Average
timecourse_heatmap_metab_avg <- as.matrix(metab_timecourse_heatmap_avg)
ht_list = Heatmap(as.matrix(t(scale(t(timecourse_heatmap_metab_avg)))),
                  top_annotation = hanno,
                  row_title = NULL,
                  show_row_names =  TRUE,
                  show_row_dend = FALSE,
                  row_names_side = "left",
                  row_names_max_width = unit(12, "cm"),
                  row_names_gp = gpar(fontsize = 8),
                  cluster_columns = FALSE)
draw(ht_list)















##############Heatmap stuff - for significant metabolites
#Signficiant set
#metab_timecourse_sign <- subset(metab_timecourse_full, metab_timecourse_full$p.adjust < 0.05)
metab_timecourse_significant <- top_n(metab_timecourse, 50, -log10(p.adjust))

metab_timecourse_heatmap_sign <- subset(metab_timecourse_significant, select = c("Avg.Egg","Avg.2cell","Avg.16cell", "Avg.Morula",
                                                                     "Avg.Blastula","Avg.Gastrula","Avg.Prism","Avg.Pluteus"))

rownames(metab_timecourse_heatmap_sign) <- gsub("X1.","1-",rownames(metab_timecourse_heatmap_sign))
rownames(metab_timecourse_heatmap_sign) <- gsub("X2.","2-",rownames(metab_timecourse_heatmap_sign))
rownames(metab_timecourse_heatmap_sign) <- gsub("X3.","3-",rownames(metab_timecourse_heatmap_sign))
rownames(metab_timecourse_heatmap_sign) <- gsub("X4.","4-",rownames(metab_timecourse_heatmap_sign))
rownames(metab_timecourse_heatmap_sign) <- gsub("X5.","5-",rownames(metab_timecourse_heatmap_sign))
rownames(metab_timecourse_heatmap_sign) <- gsub("X7.","7-",rownames(metab_timecourse_heatmap_sign))


timecourse_names <- gsub("\\_.*","",colnames(metab_timecourse_heatmap_sign))
ha = HeatmapAnnotation(timepoint = timecourse_names , annotation_name_side = "left")


#Normalized - no z-score
timecourse_heatmap_metab_sign <- as.matrix(metab_timecourse_heatmap_sign)
ht_list = Heatmap(timecourse_heatmap_metab_sign,
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE)
draw(ht_list)


#Normalized - z-score
timecourse_heatmap_metab_sign <- as.matrix(metab_timecourse_heatmap_sign)
ht_list = Heatmap(as.matrix(t(scale(t(timecourse_heatmap_metab_sign)))),
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names =  TRUE,
                  show_row_dend = FALSE,
                  row_names_side = "left",
                  row_names_max_width = unit(12, "cm"),
                  row_names_gp = gpar(fontsize = 8),
                   cluster_columns = FALSE)
draw(ht_list)




