#######Libraries#######
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


####MFuzz code####
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
