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
#BiocManager::install("marray")
library(marray)

#for PCA plotting
library(ggforce)
library(scales)

#heatmap creation
library(ComplexHeatmap)

#####Notes#####
#Used unpaired t-test for analysis to stay consistent.
####Preprocessing####
####Metaboanalyst 2019 samples####
metab_micnm <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/data_processed_2019micnm.csv")
metab_micnm <- as.data.frame(t(metab_micnm))
colnames(metab_micnm) = metab_micnm[1,]
metab_micnm <- metab_micnm[-1,]



metabolite_list <- rownames(metab_micnm)
write.csv(metabolite_list, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/metabnamenofilter.csv")



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

metabolite_list <- as.data.frame(metabolite_list)
write.csv(metabolite_list$metabolite_list, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/metabname.csv")
KEGGMOLID <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Metabolomics files/name_map.csv")

metab_micnm$Query <- metabolite_list$metabolite_list 
metab_micnm_full <- merge(metab_micnm, KEGGMOLID, by = "Query")
metab_micnm_full <- subset(metab_micnm_full, select = c("Query","MicromereC","MicromereA","MicromereB",
                                                        "Non-MicromereC","Non-MicromereA","Non-MicromereB",
                                                        "KEGG"
))

#T-test (not paired)

col <- c("MicromereC","MicromereA","MicromereB",
         "Non-MicromereC","Non-MicromereA","Non-MicromereB") 
#log(x/y) = log(x) - log(y)
metab_micnm_full$MicromereC <- as.numeric(metab_micnm_full$MicromereC)
metab_micnm_full$MicromereA <- as.numeric(metab_micnm_full$MicromereA)
metab_micnm_full$MicromereB <- as.numeric(metab_micnm_full$MicromereB)
metab_micnm_full$`Non-MicromereC` <- as.numeric(metab_micnm_full$`Non-MicromereC`)
metab_micnm_full$`Non-MicromereA` <- as.numeric(metab_micnm_full$`Non-MicromereA`)
metab_micnm_full$`Non-MicromereB` <- as.numeric(metab_micnm_full$`Non-MicromereB`)



#chose unpaired here!!!
metab_micnm_full$ttest <- apply(metab_micnm_full[,col], 1, function(x){t.test(x[c("MicromereC","MicromereA","MicromereB")], x[c("Non-MicromereC","Non-MicromereA","Non-MicromereB")], paired = FALSE)$p.value})
metab_micnm_full$logFC.Mic.m.NM <- apply(metab_micnm_full[,col], 1, function(x) log2(mean(x[c("MicromereC","MicromereA","MicromereB")])) - log2(mean(x[c("Non-MicromereC","Non-MicromereA","Non-MicromereB")])))
metab_micnm_full$GSEA.Mic.m.NM <- apply(metab_micnm_full[,c("ttest", "logFC.Mic.m.NM")],1,function(x) -log10(x[c("ttest")])/sign(x[c("logFC.Mic.m.NM")]))

#Check no names
metab_names_check <- metab_micnm_full[is.na(metab_micnm_full$KEGG),]$Query

####For app (CSV file)#####

#columns: KEGG, 2019  (exclude 2020)
metab_micnm_pathview <- metab_micnm_full[,c("KEGG","GSEA.Mic.m.NM")]

#Take mean of duplicate ID's
metab_micnm_pathview <- aggregate(.~KEGG, data = metab_micnm_pathview, FUN = mean)

#change column names
colnames(metab_micnm_pathview) <- c("KEGG", "2019")

#write csv directly into app folder
write.csv(metab_micnm_pathview, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/MIC NM/mic.m.nm.metabolomics.csv", row.names = FALSE)






#######Make plots######
#Calculate logFC
metab_micnm_full$logFC.Mic.m.NM.1 <- apply(metab_micnm_full[,col], 1, function(x) log2(x[c("MicromereA")]) - log2(x[c("Non-MicromereA")]))
metab_micnm_full$logFC.Mic.m.NM.2 <- apply(metab_micnm_full[,col], 1, function(x) log2(x[c("MicromereB")]) - log2(x[c("Non-MicromereB")]))
metab_micnm_full$logFC.Mic.m.NM.3 <- apply(metab_micnm_full[,col], 1, function(x) log2(x[c("MicromereC")]) - log2(x[c("Non-MicromereC")]))










#Of avg log2 FC
metab_micnm_fc <- metab_micnm_full[,c("Query", "logFC.Mic.m.NM")]
ggplot(metab_micnm_fc, aes(x = reorder(Query, -logFC.Mic.m.NM), logFC.Mic.m.NM)) +
  geom_segment( aes(x=reorder(Query, -logFC.Mic.m.NM), xend=Query, y=0, yend=logFC.Mic.m.NM), color="grey") +
  geom_point( color="blue", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Metabolite Name") +
  ylab("Avg.log2FC")


# Create data
data <- data.frame(
  x=LETTERS[1:26],
  y=abs(rnorm(26))
)

# Change baseline
ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=1, yend=y), color="grey") +
  geom_point( color="orange", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Value of Y")



######Preprocessing MicNM Inhibitors#######
metab_micnm_inh <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metabolomics Processed Data/Mic NM inhibitors ratio norm/mic inh processed.csv")
metab_micnm_inh <- as.data.frame(metab_micnm_inh)
colnames(metab_micnm_inh)


colnames(metab_micnm_inh) = c("Label","2DG.1","2DG.2","2DG.3","Azide.1","Azide.2","Azide.3","Cer.1","Cer.2","Cer.3","DMSO.1","DMSO.2","DMSO.3")
metab_micnm_inh<- metab_micnm_inh[-1,]

rownames(metab_micnm_inh) <- metab_micnm_inh$Label
metabolite_list_inhmicnm <- rownames(metab_micnm_inh)


metabolite_list_inhmicnm <- gsub("-TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-2TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-3TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-4TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-5TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-6TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-8TMS","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("X1-","1-",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("X2-","2-",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("X3-","3-",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("X4-","4-",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("X5-","5-",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("-meto","",metabolite_list_inhmicnm)

metabolite_list_inhmicnm <- gsub("\\(1","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("\\(2","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- gsub("\\)","",metabolite_list_inhmicnm)
metabolite_list_inhmicnm <- trimws(metabolite_list_inhmicnm)

metabolite_list_inhmicnm <- as.data.frame(metabolite_list_inhmicnm)
write.csv(metabolite_list_inhmicnm$metabolite_list_inhmicnm, file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metabolomics Processed Data/metabname_micnminh.csv")


#Need to write down programs from here...

KEGGMOLID_micinh <- read.csv(file ="C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Metabolomics Processed Data/namemap_micnminh.csv")

metab_micnm_inh$Query <- metabolite_list_inhmicnm$metabolite_list_inhmicnm
metab_micnm_inh_full <- merge(metab_micnm_inh, KEGGMOLID_micinh, by = "Query")
metab_micnm_inh_full <- subset(metab_micnm_inh_full, select = c("Query","2DG.1","2DG.2","2DG.3",
                                                                "Azide.1","Azide.2","Azide.3",
                                                                "Cer.1","Cer.2","Cer.3",
                                                                "DMSO.1","DMSO.2","DMSO.3","KEGG"
))

#T-test (not paired)

col <- c("2DG.1","2DG.2","2DG.3",
         "Azide.1","Azide.2","Azide.3",
         "Cer.1","Cer.2","Cer.3",
         "DMSO.1","DMSO.2","DMSO.3") 
#log(x/y) = log(x) - log(y)
sapply(2:13, function(i) {
  metab_micnm_inh_full[, i] <<- as.numeric(as.character(metab_micnm_inh_full[, i]))
});


#chose paired here!!!
#choose the direction stufff.
metab_micnm_inh_full$ttest.Cer.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x){t.test(x[c("Cer.1","Cer.2","Cer.3")], x[c("DMSO.1","DMSO.2","DMSO.3")], paired = FALSE)$p.value})
metab_micnm_inh_full$ttest.NaN3.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x) {tryCatch(
  expr = {t.test(x[c("Azide.1","Azide.2","Azide.3")], x[c("DMSO.1","DMSO.2","DMSO.3")], paired = FALSE)$p.value},
  error = function(e) {
    return(NA)
  }
  )})
metab_micnm_inh_full$ttest.2DG.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x) {tryCatch(
  expr = {t.test(x[c("2DG.1","2DG.2","2DG.3")], x[c("DMSO.1","DMSO.2","DMSO.3")], paired = FALSE)$p.value},
  error = function(e) {
    return(NA)
  }
)})
#^above if std dev~ 0 then return error and set result to NA.

metab_micnm_inh_full$logFC.Cer.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x) log2(mean(x[c("Cer.1","Cer.2","Cer.3")])) - log2(mean(x[c("DMSO.1","DMSO.2","DMSO.3")])))
metab_micnm_inh_full$logFC.NaN3.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x) log2(mean(x[c("Azide.1","Azide.2","Azide.3")])) - log2(mean(x[c("DMSO.1","DMSO.2","DMSO.3")])))
metab_micnm_inh_full$logFC.2DG.m.DMSO <- apply(metab_micnm_inh_full[,col], 1, function(x) log2(mean(x[c("2DG.1","2DG.2","2DG.3")])) - log2(mean(x[c("DMSO.1","DMSO.2","DMSO.3")])))

metab_micnm_inh_full$GSEA.Cer.m.DMSO <- apply(metab_micnm_inh_full[,c("ttest.Cer.m.DMSO", "logFC.Cer.m.DMSO")],1,function(x) -log10(x[c("ttest.Cer.m.DMSO")])/sign(x[c("logFC.Cer.m.DMSO")]))
metab_micnm_inh_full$GSEA.NaN3.m.DMSO <- apply(metab_micnm_inh_full[,c("ttest.NaN3.m.DMSO", "logFC.NaN3.m.DMSO")],1,function(x) -log10(x[c("ttest.NaN3.m.DMSO")])/sign(x[c("logFC.NaN3.m.DMSO")]))
metab_micnm_inh_full$GSEA.2DG.m.DMSO <- apply(metab_micnm_inh_full[,c("ttest.2DG.m.DMSO", "logFC.2DG.m.DMSO")],1,function(x) -log10(x[c("ttest.2DG.m.DMSO")])/sign(x[c("logFC.2DG.m.DMSO")]))


#####Mic NM inhibitor data for MIC/NM app#######

#columns: KEGG, Cer, NaN3, 2DG
metab_micnm_inh_pathview <- metab_micnm_inh_full[,c("KEGG","GSEA.2DG.m.DMSO","GSEA.Cer.m.DMSO","GSEA.NaN3.m.DMSO")]

#Take mean of duplicate ID's
metab_micnm_inh_pathview <- aggregate(.~KEGG, data = metab_micnm_inh_pathview, FUN = mean)

#change column names
colnames(metab_micnm_inh_pathview) <- c("KEGG", "2DG", "Cer", "NaN3")

#write csv directly into app folder
write.csv(metab_micnm_inh_pathview, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/MIC NM/mic.m.nm.inh.metabolomics.csv", row.names = FALSE)



