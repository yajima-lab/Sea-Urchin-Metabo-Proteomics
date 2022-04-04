library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(KEGGREST)
library(pathview)
library(ggrepel)
library(ggfortify)
library(cluster)
library(tidyverse)

vasa_MO_file6 <- read.delim(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Vasa_MO_file6.txt")

#strsplit the Cont and Vasa Columns
ContMO1day.df = do.call("rbind", strsplit(vasa_MO_file6$Value.ContMO1day, ";"))
ContMO1day.df = data.frame(apply(ContMO1day.df, 2, as.numeric))
names(ContMO1day.df) = paste("ContMO1day.", 1:3, sep = "")
head(ContMO1day.df)

VasaMO1day.df = do.call("rbind", strsplit(vasa_MO_file6$Value.VasaMO1day, ";"))
VasaMO1day.df = data.frame(apply(VasaMO1day.df, 2, as.numeric))
names(VasaMO1day.df) = paste("VasaMO1day.", 1:3, sep = "")
head(VasaMO1day.df)


ContMO16cell.df = do.call("rbind", strsplit(vasa_MO_file6$Value.ContMO16cell, ";"))
ContMO16cell.df = data.frame(apply(ContMO16cell.df, 2, as.numeric))
names(ContMO16cell.df) = paste("ContMO16cell.", 1:3, sep = "")
head(ContMO16cell.df)

VasaMO16cell.df = do.call("rbind", strsplit(vasa_MO_file6$Value.VasaMO16cell, ";"))
VasaMO16cell.df = data.frame(apply(VasaMO16cell.df, 2, as.numeric))
names(VasaMO16cell.df) = paste("VasaMO16cell.", 1:3, sep = "")
head(VasaMO16cell.df)



###Combine and take into account only proteins with values for all 3 samples
VasaMO_full <- cbind(vasa_MO_file6, VasaMO16cell.df, ContMO16cell.df, VasaMO1day.df, ContMO1day.df)

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

VasaMO_full <- completeFun(VasaMO_full, c("VasaMO1day.1","ContMO1day.1","VasaMO1day.2","ContMO1day.2","VasaMO1day.3","ContMO1day.3",
                                          "VasaMO16cell.1","ContMO16cell.1","VasaMO16cell.2","ContMO16cell.2","VasaMO16cell.3","ContMO16cell.3"))



###Get Log2FC for Mic-NM
col <- c("VasaMO1day.1","ContMO1day.1","VasaMO1day.2","ContMO1day.2","VasaMO1day.3","ContMO1day.3",
         "VasaMO16cell.1","ContMO16cell.1","VasaMO16cell.2","ContMO16cell.2","VasaMO16cell.3","ContMO16cell.3") 

#log(x/y) = log(x) - log(y)
VasaMO_full$logFC.Vasa1day.m.Cont <- apply(VasaMO_full[,col], 1, function(x) log2(mean(x[c("VasaMO1day.1","VasaMO1day.2","VasaMO1day.3")])) - log2(mean(x[c("ContMO1day.1","ContMO1day.2","ContMO1day.3")])))
VasaMO_full$logFC.Vasa16cell.m.Cont <- apply(VasaMO_full[,col], 1, function(x) log2(mean(x[c("VasaMO16cell.1","VasaMO16cell.2","VasaMO16cell.3")])) - log2(mean(x[c("ContMO16cell.1","ContMO16cell.2","ContMO16cell.3")])))

#paired t-test
VasaMO_full$ttest.Vasa1day.m.Cont <- apply(VasaMO_full[,col], 1, function(x){t.test(x[c("VasaMO1day.1","VasaMO1day.2","VasaMO1day.3")], x[c("ContMO1day.1","ContMO1day.2","ContMO1day.3")], paired = TRUE)$p.value})
VasaMO_full$ttest.Vasa16cell.m.Cont <- apply(VasaMO_full[,col], 1, function(x){t.test(x[c("VasaMO16cell.1","VasaMO16cell.2","VasaMO16cell.3")], x[c("ContMO16cell.1","ContMO16cell.2","ContMO16cell.3")], paired = TRUE)$p.value})


#GSEA
VasaMO_full$GSEA.Vasa1day.m.Cont <- apply(VasaMO_full[,c("logFC.Vasa1day.m.Cont", "ttest.Vasa1day.m.Cont")],1,function(x) -log10(x[c("ttest.Vasa1day.m.Cont")])/sign(x[c("logFC.Vasa1day.m.Cont")]))
VasaMO_full$GSEA.Vasa16cell.m.Cont <- apply(VasaMO_full[,c("logFC.Vasa16cell.m.Cont", "ttest.Vasa16cell.m.Cont")],1,function(x) -log10(x[c("ttest.Vasa16cell.m.Cont")])/sign(x[c("logFC.Vasa16cell.m.Cont")]))


#######PCA (based on most variable genes, etc.) 16 and 1 day#######
vasa_samples <- as.data.frame(t(VasaMO_full[,c("VasaMO1day.1","ContMO1day.1","VasaMO1day.2","ContMO1day.2","VasaMO1day.3","ContMO1day.3",
                                                 "VasaMO16cell.1","ContMO16cell.1","VasaMO16cell.2","ContMO16cell.2","VasaMO16cell.3","ContMO16cell.3")]))


rv <- rowVars(as.matrix(t(vasa_samples)))
hist(rv)
select = order(rv, decreasing=TRUE)[seq_len(500)]
#most_variable_proteins <- timecourse_samples[ ,select]
pca_mv = prcomp(vasa_samples[ ,select], center = TRUE, scale = TRUE)
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


######PCA (based on most variable genes) 1 day####
vasa_samples_1day <- as.data.frame(t(VasaMO_full[,c("VasaMO1day.1","ContMO1day.1","VasaMO1day.2","ContMO1day.2","VasaMO1day.3","ContMO1day.3")]))


rv <- rowVars(as.matrix(t(vasa_samples_1day)))
hist(rv)
select = order(rv, decreasing=TRUE)[seq_len(1000)]
#most_variable_proteins <- timecourse_samples[ ,select]
pca_mv = prcomp(vasa_samples_1day[ ,select], center = TRUE, scale = TRUE)
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


###### PCA, 16 cell #####
vasa_samples <- as.data.frame(t(VasaMO_full[,c("VasaMO1day.1","ContMO1day.1","VasaMO1day.2","ContMO1day.2","VasaMO1day.3","ContMO1day.3",
                                               "VasaMO16cell.1","ContMO16cell.1","VasaMO16cell.2","ContMO16cell.2","VasaMO16cell.3","ContMO16cell.3")]))




rv <- rowVars(as.matrix(t(vasa_samples_1day)))
hist(rv)
select = order(rv, decreasing=TRUE)[seq_len(1000)]
#most_variable_proteins <- timecourse_samples[ ,select]
pca_mv = prcomp(vasa_samples_1day[ ,select], center = TRUE, scale = TRUE)
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
