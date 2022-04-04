library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(KEGGREST)
library(pathview)


####Common Metabolites that aren't converted easily: {manually edited}
#3-Aminoglutaric Acid -> isoglutamte -> C05574
#Ribose 5 phosphate -> C00117

####Things to focus on
#add the correct label of glucose 6 phosphate (appropriate label..)

###All p-values are PAIRED (b/c we want to understand differences across each sample...)


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


####Function dataprocessing
metab_stat_WE <- function(data, KEGGconv){
  data_full <- merge(data, KEGGconv, by = "Query")
  
  data_full <- subset(data_full, select = c("Query","2DG-1","2DG-2","2DG-3", 
                                            "Cer-1","Cer-2","Cer-3", 
                                            "DMSO-1","DMSO-2","DMSO-3",
                                            "NaN3-1","NaN3-2","NaN3-3","KEGG"))
  
  #Paired T-test
  col <- c("2DG-1","2DG-2","2DG-3", 
           "Cer-1","Cer-2","Cer-3", 
           "DMSO-1","DMSO-2","DMSO-3",
           "NaN3-1","NaN3-2","NaN3-3"
  ) 
  
  #log(x/y) = log(x) - log(y)
  sapply(2:13, function(i) {
    data_full[, i] <<- as.numeric(as.character(data_full[, i]))
  });
  
  data_full$ttest.2DG.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("2DG-1","2DG-2","2DG-3")], x[c("DMSO-1","DMSO-2","DMSO-3")], paired = TRUE)$p.value})
  data_full$ttest.Cer.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("Cer-1","Cer-2","Cer-3")], x[c("DMSO-1","DMSO-2","DMSO-3")], paired = TRUE)$p.value})
  data_full$ttest.NaN3.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("NaN3-1","NaN3-2","NaN3-3")], x[c("DMSO-1","DMSO-2","DMSO-3")], paired = TRUE)$p.value})
  
  data_full$logFC.2DG.m.DMSO  <- apply( data_full[,col], 1, function(x) log2(mean(x[c("2DG-1","2DG-2","2DG-3")])) - log2(mean(x[c("DMSO-1","DMSO-2","DMSO-3")])))
  data_full$logFC.Cer.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("Cer-1","Cer-2","Cer-3")])) - log2(mean(x[c("DMSO-1","DMSO-2","DMSO-3")])))
  data_full$logFC.NaN3.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("NaN3-1","NaN3-2","NaN3-3")])) - log2(mean(x[c("DMSO-1","DMSO-2","DMSO-3")])))
  
  data_full$GSEA.2DG.m.DMSO <- apply( data_full[,c("ttest.2DG.m.DMSO", "logFC.2DG.m.DMSO")],1,function(x) -log10(x[c("ttest.2DG.m.DMSO")])/sign(x[c("logFC.2DG.m.DMSO")]))
  data_full$GSEA.Cer.m.DMSO <- apply( data_full[,c("ttest.Cer.m.DMSO", "logFC.Cer.m.DMSO")],1,function(x) -log10(x[c("ttest.Cer.m.DMSO")])/sign(x[c("logFC.Cer.m.DMSO")]))
  data_full$GSEA.NaN3.m.DMSO <- apply( data_full[,c("ttest.NaN3.m.DMSO", "logFC.NaN3.m.DMSO")],1,function(x) -log10(x[c("ttest.NaN3.m.DMSO")])/sign(x[c("logFC.NaN3.m.DMSO")]))
  
  
  data_full_pathview <- subset(data_full, select = c("GSEA.2DG.m.DMSO","GSEA.Cer.m.DMSO","GSEA.NaN3.m.DMSO","KEGG"))
  
  
  ####Aggregate by average 
  #Take mean of duplicate ID's
  data_full_pathview <- aggregate(.~KEGG, data = data_full_pathview, FUN = mean,na.action = na.pass)
  rownames(data_full_pathview) <- data_full_pathview$KEGG
  data_full_pathview$KEGG <- NULL
  return(data_full_pathview)
}
metab_stat_micinh <- function(data, KEGGconv){
  data_full <- merge(data, KEGGconv, by = "Query")
  
  data_full <- subset(data_full, select = c("Query",
                                            "2DG-Mic-1","2DG-Mic-2","2DG-Mic-3",
                                            "Cer-Mic-1","Cer-Mic-2","Cer-Mic-3",
                                            "DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3",
                                            "NaN3-Mic-1","NaN3-Mic-2","NaN3-Mic-3",
                                            "KEGG"))
  
  #Paired T-test
  col <- c("2DG-Mic-1","2DG-Mic-2","2DG-Mic-3",
           "Cer-Mic-1","Cer-Mic-2","Cer-Mic-3",
           "DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3",
           "NaN3-Mic-1","NaN3-Mic-2","NaN3-Mic-3") 
  
  #log(x/y) = log(x) - log(y)
  sapply(2:13, function(i) {
    data_full[, i] <<- as.numeric(as.character(data_full[, i]))
  });
  
  data_full$ttest.2DG.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("2DG-Mic-1","2DG-Mic-2","2DG-Mic-3")], x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")], paired = TRUE)$p.value})
  data_full$ttest.Cer.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("Cer-Mic-1","Cer-Mic-2","Cer-Mic-3")], x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")], paired = TRUE)$p.value})
  data_full$ttest.NaN3.m.DMSO  <- apply( data_full[,col], 1, function(x){t.test(x[c("NaN3-Mic-1","NaN3-Mic-2","NaN3-Mic-3")], x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")], paired = TRUE)$p.value})
  
  data_full$logFC.2DG.m.DMSO  <- apply( data_full[,col], 1, function(x) log2(mean(x[c("2DG-Mic-1","2DG-Mic-2","2DG-Mic-3")])) - log2(mean(x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")])))
  data_full$logFC.Cer.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("Cer-Mic-1","Cer-Mic-2","Cer-Mic-3")])) - log2(mean(x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")])))
  data_full$logFC.NaN3.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("NaN3-Mic-1","NaN3-Mic-2","NaN3-Mic-3")])) - log2(mean(x[c("DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3")])))
  
  data_full$GSEA.2DG.m.DMSO <- apply( data_full[,c("ttest.2DG.m.DMSO", "logFC.2DG.m.DMSO")],1,function(x) -log10(x[c("ttest.2DG.m.DMSO")])/sign(x[c("logFC.2DG.m.DMSO")]))
  data_full$GSEA.Cer.m.DMSO <- apply( data_full[,c("ttest.Cer.m.DMSO", "logFC.Cer.m.DMSO")],1,function(x) -log10(x[c("ttest.Cer.m.DMSO")])/sign(x[c("logFC.Cer.m.DMSO")]))
  data_full$GSEA.NaN3.m.DMSO <- apply( data_full[,c("ttest.NaN3.m.DMSO", "logFC.NaN3.m.DMSO")],1,function(x) -log10(x[c("ttest.NaN3.m.DMSO")])/sign(x[c("logFC.NaN3.m.DMSO")]))
  
  
  data_full_pathview <- subset(data_full, select = c("GSEA.2DG.m.DMSO","GSEA.Cer.m.DMSO","GSEA.NaN3.m.DMSO","KEGG"))
  
  
  ####Aggregate by average 
  #Take mean of duplicate ID's
  data_full_pathview <- aggregate(.~KEGG, data = data_full_pathview, FUN = mean,na.action = na.pass)
  rownames(data_full_pathview) <- data_full_pathview$KEGG
  data_full_pathview$KEGG <- NULL
  
  return(data_full_pathview)
  
  #check above...
}
metab_stat_micnmr_inh <- function(data, KEGGconv){
  data_full <- merge(data, KEGGconv, by = "Query")
  data_full <- subset(data_full, select = c("Query",
                                            "2DG-Mic/NM-1","2DG-Mic/NM-2","2DG-Mic/NM-3",
                                            "Cer-Mic/NM-1","Cer-Mic/NM-2","Cer-Mic/NM-3",
                                            "DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3",
                                            "NaN3-Mic/NM-1","NaN3-Mic/NM-2","NaN3-Mic/NM-3","KEGG"))
  
  #Paired T-test
  col <- c("2DG-Mic/NM-1","2DG-Mic/NM-2","2DG-Mic/NM-3",
           "Cer-Mic/NM-1","Cer-Mic/NM-2","Cer-Mic/NM-3",
           "DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3",
           "NaN3-Mic/NM-1","NaN3-Mic/NM-2","NaN3-Mic/NM-3") 
  
  #log(x/y) = log(x) - log(y)
  sapply(2:13, function(i) {
    data_full[, i] <<- as.numeric(as.character(data_full[, i]))
  });
  
  data_full$ttest.2DG.m.DMSO  <- apply(data_full[,col], 1, function(x){tryCatch(
    expr = {t.test(x[c("2DG-Mic/NM-1","2DG-Mic/NM-2","2DG-Mic/NM-3")], x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")], paired = TRUE)$p.value},
    error = function(e){
      return(NA)
    })})
  data_full$ttest.Cer.m.DMSO  <- apply(data_full[,col], 1, function(x){tryCatch(
    expr = {t.test(x[c("Cer-Mic/NM-1","Cer-Mic/NM-2","Cer-Mic/NM-3")], x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")], paired = TRUE)$p.value},
    error = function(e) {
      return(NA)})})
  
  data_full$ttest.NaN3.m.DMSO  <- apply(data_full[,col], 1, function(x){tryCatch(
    expr = {t.test(x[c("NaN3-Mic/NM-1","NaN3-Mic/NM-2","NaN3-Mic/NM-3")], x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")], paired = TRUE)$p.value},
    error = function(e){
      return(NA)})})
  
  data_full$logFC.2DG.m.DMSO  <- apply( data_full[,col], 1, function(x) log2(mean(x[c("2DG-Mic/NM-1","2DG-Mic/NM-2","2DG-Mic/NM-3")])) - log2(mean(x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")])))
  data_full$logFC.Cer.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("Cer-Mic/NM-1","Cer-Mic/NM-2","Cer-Mic/NM-3")])) - log2(mean(x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")])))
  data_full$logFC.NaN3.m.DMSO <- apply( data_full[,col], 1, function(x) log2(mean(x[c("NaN3-Mic/NM-1","NaN3-Mic/NM-2","NaN3-Mic/NM-3")])) - log2(mean(x[c("DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3")])))
  
  data_full$GSEA.2DG.m.DMSO <- apply( data_full[,c("ttest.2DG.m.DMSO", "logFC.2DG.m.DMSO")],1,function(x) -log10(x[c("ttest.2DG.m.DMSO")])/sign(x[c("logFC.2DG.m.DMSO")]))
  data_full$GSEA.Cer.m.DMSO <- apply( data_full[,c("ttest.Cer.m.DMSO", "logFC.Cer.m.DMSO")],1,function(x) -log10(x[c("ttest.Cer.m.DMSO")])/sign(x[c("logFC.Cer.m.DMSO")]))
  data_full$GSEA.NaN3.m.DMSO <- apply( data_full[,c("ttest.NaN3.m.DMSO", "logFC.NaN3.m.DMSO")],1,function(x) -log10(x[c("ttest.NaN3.m.DMSO")])/sign(x[c("logFC.NaN3.m.DMSO")]))
  
  
  data_full_pathview <- subset(data_full, select = c("GSEA.2DG.m.DMSO","GSEA.Cer.m.DMSO","GSEA.NaN3.m.DMSO","KEGG"))
  
  ####Aggregate by average 
  #Take mean of duplicate ID's
  
  data_full_pathview <- aggregate(.~KEGG, data = data_full_pathview, FUN = mean,na.action = na.pass)
  rownames(data_full_pathview) <- data_full_pathview$KEGG
  data_full_pathview$KEGG <- NULL
  
  return(data_full_pathview)
}
#metab_stat_micnm 
#No stat for timecourse just used the normalized values to form heatmap...




#Confirmed want Paired T-test for these comparisons
#Need to make global metabolite name matcher:
#######16 cell WE inh#########
metab_16cell_inh <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/inhibitors 2021/inh_16cell_dataprocessed.csv")
metab_16cell_inh <- as.data.frame(t(metab_16cell_inh))
colnames(metab_16cell_inh) = metab_16cell_inh[1,]
metab_16cell_inh <- metab_16cell_inh[-1,]
metab_16cell_inh <- metab_16cell_inh[-1,]
metab_16cell_inh$Query <- metabolite_clean1(rownames(metab_16cell_inh))


write.csv(metab_16cell_inh$Query, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/metabnameinh_16cell.csv")
KEGGMOLID <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_16cell.csv")


metab_16cell_inh_pathview <- metab_stat_WE(metab_16cell_inh, KEGGMOLID)

###Add more speicfic KEGG labels
metab_16cell_inh_pathview <- metabolite_alpha(metab_16cell_inh_pathview)

#########2 day WE inh##########
metab_2day_inh <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/inhibitors 2021/inh_2day_dataprocessed.csv")
metab_2day_inh <- as.data.frame(t(metab_2day_inh))
colnames(metab_2day_inh) = metab_2day_inh[1,]
metab_2day_inh <- metab_2day_inh[-1,]
metab_2day_inh <- metab_2day_inh[-1,]
metab_2day_inh$Query <- metabolite_clean1(rownames(metab_2day_inh))

write.csv(metab_2day_inh$Query, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/metabnameinh_2day.csv")
KEGGMOLID_2 <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_2day.csv")

metab_2day_inh_pathview <- metab_stat_WE(metab_2day_inh, KEGGMOLID_2)

###Add more specific KEGG labels
metab_2day_inh_pathview <- metabolite_alpha(metab_2day_inh_pathview)

########MIC/NM set (Mic vs DMSO)#########
metab_micnm_inh <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/micnm inhibitors 2021/micinh_data_processed.csv")
metab_micnm_inh <- as.data.frame(t(metab_micnm_inh))
colnames(metab_micnm_inh) = metab_micnm_inh[1,]
metab_micnm_inh <- metab_micnm_inh[-1,]
metab_micnm_inh <- metab_micnm_inh[-1,]
metab_micnm_inh$Query <- metabolite_clean1(rownames(metab_micnm_inh))

write.csv(metab_micnm_inh$Query, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/metabname_micnm_inh.csv")
KEGGMOLID_3 <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_micnminh.csv")

metab_micvDMSO_inh_pathview <- metab_stat_micinh(metab_micnm_inh,KEGGMOLID_3)
###Add more specific KEGG labels

metab_micvDMSO_inh_pathview<- metabolite_alpha(metab_micvDMSO_inh_pathview)



#####ratio Mic/NM vs DMSO#####
sapply(1:24, function(i) {
  metab_micnm_inh[, i] <<- as.numeric(as.character(metab_micnm_inh[, i]))
});
length(metab_micnm_inh$Query) == unique(length(metab_micnm_inh$Query))

metab_micnm_inh_mic <- subset(metab_micnm_inh, select = c("2DG-Mic-1","2DG-Mic-2","2DG-Mic-3",
                                                          "Cer-Mic-1","Cer-Mic-2","Cer-Mic-3",
                                                          "DMSO-Mic-1","DMSO-Mic-2","DMSO-Mic-3",
                                                          "NaN3-Mic-1","NaN3-Mic-2","NaN3-Mic-3","Query"))

metab_micnm_inh_nm <- subset(metab_micnm_inh, select = c("2DG-NM-1","2DG-NM-2","2DG-NM-3",
                                                         "Cer-NM-1","Cer-NM-2","Cer-NM-3",
                                                         "DMSO-NM-1","DMSO-NM-2","DMSO-NM-3",
                                                         "NaN3-NM-1","NaN3-NM-2","NaN3-Nm-3","Query"))

metab_micnm_inh_ratio <- metab_micnm_inh_mic[,1:12]/metab_micnm_inh_nm[,1:12]
colnames(metab_micnm_inh_ratio) <- c("2DG-Mic/NM-1","2DG-Mic/NM-2","2DG-Mic/NM-3",
                                     "Cer-Mic/NM-1","Cer-Mic/NM-2","Cer-Mic/NM-3",
                                     "DMSO-Mic/NM-1","DMSO-Mic/NM-2","DMSO-Mic/NM-3",
                                     "NaN3-Mic/NM-1","NaN3-Mic/NM-2","NaN3-Mic/NM-3")
metab_micnm_inh_ratio$Query <- metabolite_clean1(rownames(metab_micnm_inh_ratio))

metab_micnmratio_pathview <- metab_stat_micnmr_inh(metab_micnm_inh_ratio, KEGGMOLID_3)
#Used as a check
#metab_micnm_inh$'2DG-ratio-1' <- metab_micnm_inh$`2DG-Mic-1`/metab_micnm_inh$`2DG-NM-1`


###Add more specific KEGG labels
metab_micnmratio_pathview <- metabolite_alpha(metab_micnmratio_pathview)

######MIC/NM 2019 sample######
metab_micnm <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/data_processed_2019micnm.csv")
metab_micnm <- as.data.frame(t(metab_micnm))
colnames(metab_micnm) = metab_micnm[1,]
metab_micnm <- metab_micnm[-1,]
metab_micnm$Query <- metabolite_clean1(rownames(metab_micnm))

write.csv(metab_micnm$Query, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/metabname_micnm2019.csv")

KEGGMOLID_4 <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_micnm.csv")

metab_micnm_full <- merge(metab_micnm, KEGGMOLID_4, by = "Query")
metab_micnm_full <- subset(metab_micnm_full, select = c("MicromereC","MicromereA","MicromereB",
                                                        "Non-MicromereC","Non-MicromereA","Non-MicromereB",
                                                        "KEGG"
))

#Paired T-test
col <- c("MicromereC","MicromereA","MicromereB",
         "Non-MicromereC","Non-MicromereA","Non-MicromereB") 
sapply(1:6, function(i) {
  metab_micnm_full[, i] <<- as.numeric(as.character(metab_micnm_full[, i]))
});


metab_micnm_full$ttest <- apply(metab_micnm_full[,col], 1, function(x){t.test(x[c("MicromereC","MicromereA","MicromereB")], x[c("Non-MicromereC","Non-MicromereA","Non-MicromereB")], paired = TRUE)$p.value})
metab_micnm_full$logFC.Mic.m.NM <- apply(metab_micnm_full[,col], 1, function(x) log2(mean(x[c("MicromereC","MicromereA","MicromereB")])) - log2(mean(x[c("Non-MicromereC","Non-MicromereA","Non-MicromereB")])))
metab_micnm_full$GSEA.Mic.m.NM <- apply(metab_micnm_full[,c("ttest", "logFC.Mic.m.NM")],1,function(x) -log10(x[c("ttest")])/sign(x[c("logFC.Mic.m.NM")]))

metab_micnm_full_pathview <- subset(metab_micnm_full, select = c("KEGG","GSEA.Mic.m.NM"))
metab_micnm_full_pathview <- aggregate(.~KEGG, data = metab_micnm_full_pathview, FUN = mean,na.action = na.pass)
rownames(metab_micnm_full_pathview) <- metab_micnm_full_pathview$KEGG
metab_micnm_full_pathview$KEGG <- NULL

###Add more specific KEGG labels
#for glucose-6-phosphate
metab_micnm_full_pathview <- metabolite_alpha(metab_micnm_full_pathview)



######Timecourse######
#Using the normalized, quantile data {the heatmapped version}
metab_timecourse <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/Timecourse/timeseries_normalized_quantile.csv")
metab_timecourse <- as.data.frame(t(metab_timecourse))
colnames(metab_timecourse) <- metab_timecourse[1,]
metab_timecourse<- metab_timecourse[-1,]
metab_timecourse<- metab_timecourse[-1,]
metab_timecourse$Query <- metabolite_clean1(rownames(metab_timecourse))
write.csv(metab_timecourse$Query, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/metabname_timecourse.csv")
KEGGMOLID_5 <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_timecourse.csv")

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

####Aggregate by average 
#Take mean of duplicate ID's
metab_timecourse_full <- aggregate(.~KEGG, data = metab_timecourse_full, FUN = mean, na.action = na.pass)


col <- colnames(metab_timecourse_full)[2:25]
metab_timecourse_full$Avg.Egg <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Egg_A", "Egg_B", "Egg_C")]))
metab_timecourse_full$Avg.2cell <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("2 cell_A", "2 cell_B", "2 cell_C")]))
metab_timecourse_full$Avg.16cell <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("16 cell_A", "16 cell_B", "16 cell_C")]))
metab_timecourse_full$Avg.Morula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Morula_A", "Morula_B", "Morula_C")]))
metab_timecourse_full$Avg.Blastula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Blastula_A", "Blastula_B", "Blastula_C")]))
metab_timecourse_full$Avg.Gastrula <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Gastrula_A", "Gastrula_B", "Gastrula_C")]))
metab_timecourse_full$Avg.Prism <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Prism_A", "Prism_B", "Prism_C")]))
metab_timecourse_full$Avg.Pluteus <- apply(metab_timecourse_full[,col],1,function(x) mean(x[c("Pluteus_A", "Pluteus_B", "Pluteus_C")]))

metab_timecourse_full_pathview <- subset(metab_timecourse_full, select = c("KEGG","Avg.Egg","Avg.2cell","Avg.16cell","Avg.Morula","Avg.Blastula",
                                                                          "Avg.Gastrula","Avg.Prism","Avg.Pluteus"))
rownames(metab_timecourse_full_pathview) <- metab_timecourse_full_pathview$KEGG
metab_timecourse_full_pathview$KEGG <- NULL

###Modified ID changes
metab_timecourse_full_pathview <- metabolite_alpha(metab_timecourse_full_pathview)

######Timecourse Egg scaled########
egg_scaling <- function(vector){
  norm <- (vector - vector[1])/(max(vector) - min(vector))
  return(norm)
} 

#Test
#test1 <- c(7,1,10,5,5,9,12,0.1)
#egg_scaling(test1)

metab_timecourse_full_eggscale_pathview <- as.data.frame(t(apply(metab_timecourse_full_pathview, 1, egg_scaling)))



#####Consolidating the ID conversion files together to check#####
ALL_ID_conversions <- do.call("rbind", list(KEGGMOLID, KEGGMOLID_2, KEGGMOLID_3, KEGGMOLID_4, KEGGMOLID_5))
ALL_ID_conversions <- unique(ALL_ID_conversions)
write.csv(ALL_ID_conversions, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metab Conv/name_map_alldatasets_manu.csv")
######Global Pathways#####
#Items to Prepare:
#May Replace with new color (figure dat S1 - S10) in the folder!





#mic/nm 2019 samples
global_metab_micnm_2019 <- pathview(cpd.data = metab_micnm_full_pathview,
                       pathway.id = 'spu01100',
                       species    = 'spu',
                       out.suffix = "micnm2019",
                       low = list(cpd = "blue"),
                       mid = list(cpd = "gray"),
                       high = list(cpd = "red"))

#DG <- subset(metab_micnmratio_pathview, select = c("GSEA.2DG.m.DMSO"))

#mic/nm ratio inh
global_metab_micnmratio_inh_2DG <- pathview(cpd.data =  subset(metab_micnmratio_pathview, select = c("GSEA.2DG.m.DMSO")),
                                            pathway.id = 'spu01100',
                                            species    = 'spu',
                                            out.suffix = "micnminh2DG",
                                            low = list(cpd = "blue"),
                                            mid = list(cpd = "gray"),
                                            high = list(cpd = "red"))

global_metab_micnmratio_inh_Cer <- pathview(cpd.data =  subset(metab_micnmratio_pathview, select = c("GSEA.Cer.m.DMSO")),
                                            pathway.id = 'spu01100',
                                            species    = 'spu',
                                            out.suffix = "micnminhCer",
                                            low = list(cpd = "blue"),
                                            mid = list(cpd = "gray"),
                                            high = list(cpd = "red"))

global_metab_micnmratio_inh_NaN3 <- pathview(cpd.data =  subset(metab_micnmratio_pathview, select = c("GSEA.NaN3.m.DMSO")),
                                            pathway.id = 'spu01100',
                                            species    = 'spu',
                                            out.suffix = "micnminhNaN3",
                                            low = list(cpd = "blue"),
                                            mid = list(cpd = "gray"),
                                            high = list(cpd = "red"))

#WE 2day inh
global_metab_2day_inh_2DG <- pathview(cpd.data =  subset(metab_2day_inh_pathview, select = c("GSEA.2DG.m.DMSO")),
                                     pathway.id = 'spu01100',
                                     species    = 'spu',
                                     out.suffix = "2dayinh2DG",
                                     low = list(cpd = "blue"),
                                     mid = list(cpd = "gray"),
                                     high = list(cpd = "red"))

global_metab_2day_inh_Cer <- pathview(cpd.data =  subset(metab_2day_inh_pathview, select = c("GSEA.Cer.m.DMSO")),
                                      pathway.id = 'spu01100',
                                      species    = 'spu',
                                      out.suffix = "2dayinhCer",
                                      low = list(cpd = "blue"),
                                      mid = list(cpd = "gray"),
                                      high = list(cpd = "red"))

global_metab_2day_inh_NaN3 <- pathview(cpd.data =  subset(metab_2day_inh_pathview, select = c("GSEA.NaN3.m.DMSO")),
                                      pathway.id = 'spu01100',
                                      species    = 'spu',
                                      out.suffix = "2dayinhNaN3",
                                      low = list(cpd = "blue"),
                                      mid = list(cpd = "gray"),
                                      high = list(cpd = "red"))

#WE 16 cell inh
global_metab_16cell_inh_2DG <- pathview(cpd.data =  subset(metab_16cell_inh_pathview, select = c("GSEA.2DG.m.DMSO")),
                                      pathway.id = 'spu01100',
                                      species    = 'spu',
                                      out.suffix = "16cellinh2DG",
                                      low = list(cpd = "blue"),
                                      mid = list(cpd = "gray"),
                                      high = list(cpd = "red"))

global_metab_16cell_inh_Cer <- pathview(cpd.data =  subset(metab_16cell_inh_pathview, select = c("GSEA.Cer.m.DMSO")),
                                      pathway.id = 'spu01100',
                                      species    = 'spu',
                                      out.suffix = "16cellinhCer",
                                      low = list(cpd = "blue"),
                                      mid = list(cpd = "gray"),
                                      high = list(cpd = "red"))

global_metab_16cell_inh_NaN3 <- pathview(cpd.data =  subset(metab_16cell_inh_pathview, select = c("GSEA.NaN3.m.DMSO")),
                                       pathway.id = 'spu01100',
                                       species    = 'spu',
                                       out.suffix = "16cellinhNaN3",
                                       low = list(cpd = "blue"),
                                       mid = list(cpd = "gray"),
                                       high = list(cpd = "red"))


###Timecourse Egg scaled norm.
global_metab_timecourse_2cell_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.2cell")),
                                        pathway.id = 'spu01100',
                                        species    = 'spu',
                                        out.suffix = "2cellscaleegg",
                                        low = list(cpd = "blue"),
                                        mid = list(cpd = "gray"),
                                        high = list(cpd = "red"))

global_metab_timecourse_16cell_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.16cell")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "16cellscaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))

global_metab_timecourse_morula_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.Morula")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "morulascaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))

global_metab_timecourse_blastula_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.Blastula")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "blastulascaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))

global_metab_timecourse_gastrula_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.Gastrula")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "gastrulascaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))

global_metab_timecourse_prism_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.Prism")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "prismscaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))

global_metab_timecourse_pluteus_eggscale <- pathview(cpd.data =  subset(metab_timecourse_full_eggscale_pathview, select = c("Avg.Pluteus")),
                                                   pathway.id = 'spu01100',
                                                   species    = 'spu',
                                                   out.suffix = "pluteusscaleegg",
                                                   low = list(cpd = "blue"),
                                                   mid = list(cpd = "gray"),
                                                   high = list(cpd = "red"))













###Timecourse old:
global_metab_timecourse_EGG <- pathview(cpd.data =  subset(metab_timecourse_full_pathview, select = c("Avg.Egg")),
                                      pathway.id = 'spu01100',
                                      species    = 'spu',
                                      low = list(cpd = "blue"),
                                      mid = list(cpd = "gray"),
                                      high = list(cpd = "red"))

global_metab_timecourse_2cell <- pathview(cpd.data =  subset(metab_timecourse_full_pathview, select = c("Avg.2cell")),
                                        pathway.id = 'spu01100',
                                        species    = 'spu',
                                        low = list(cpd = "blue"),
                                        mid = list(cpd = "gray"),
                                        high = list(cpd = "red"))




###Test to see if the cloud version is working
  


#######Specific Pathways######
glycolysis <- pathview(cpd.data =  subset(metab_micnmratio_pathview, select = c("GSEA.2DG.m.DMSO")),
                       pathway.id = 'spu00010',
                       species    = 'spu',
                       low = list(cpd = "blue"),
                       mid = list(cpd = "gray"),
                       high = list(cpd = "red"))

glycolysis <- pathview(cpd.data = metab_micnm_full_pathview,
                       pathway.id = 'spu00010',
                       species    = 'spu',
                       out.suffix = "try1",
                       low = list(cpd = "blue"),
                       mid = list(cpd = "gray"),
                       high = list(cpd = "green"))




######Pathways######
glycolysis <- pathview(cpd.data = metab_16cell_inh_pathview,
                       pathway.id = 'spu00010',
                       species    = 'spu',
                       low = list(cpd = "blue"),
                       mid = list(cpd = "gray"),
                       high = list(cpd = "red"))

glycolysis <- pathview(cpd.data = metab_16cell_inh_pathview,
                       pathway.id = 'spu00010',
                       species    = 'spu',
                       low = list(cpd = "blue"),
                       mid = list(cpd = "gray"),
                       high = list(cpd = "green"))


######Making CSV Files for Apps######
#MIC vs NM
write.csv(metab_micnm_full_pathview, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/MIC NM/mic.m.nm.metabolomics.csv")
#Mic vs NM inh
write.csv(metab_micnmratio_pathview, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/MIC NM/mic.m.nm.inh.metabolomics.csv")

#Timecourse
write.csv(metab_timecourse_full_pathview, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Data/timecourse_metabolomics_global.csv")

