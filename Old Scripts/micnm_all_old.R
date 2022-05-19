####Libraries####
library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggupset)
library(ggridges)
library(KEGGREST)
library(ggrepel)
#for PCA plotting
library(ggforce)
library(scales)
#heatmap creation
library(ComplexHeatmap)

#######Data Preprocessing#####
#Original File
mic_nm_file6 <- read.delim(file = "Data/prot_data/Mic_NM_file6.txt")

#strsplit the MIC and NM columns and re-add these columns back into the main file
NM.df = do.call("rbind", strsplit(mic_nm_file6$Value.NM, ";"))
NM.df = data.frame(apply(NM.df, 2, as.numeric))
names(NM.df) = paste("NM.", 1:3, sep = "")
head(NM.df)

MIC.df = do.call("rbind", strsplit(mic_nm_file6$Value.MIC, ";"))
MIC.df = data.frame(apply(MIC.df, 2, as.numeric))
names(MIC.df) = paste("MIC.", 1:3, sep = "")
head(MIC.df)

mic_nm_full <- cbind(mic_nm_file6, MIC.df, NM.df)

###Get Log2FC for Mic-NM
col <- c("MIC.1","NM.1","MIC.2","NM.2","MIC.3","NM.3") 
#log(x/y) = log(x) - log(y)
mic_nm_full$logFC.Mic.m.NM <- apply(mic_nm_full[,col], 1, function(x) log2(mean(x[c("MIC.1","MIC.2","MIC.3")])) - log2(mean(x[c("NM.1","NM.2","NM.3")])))
mic_nm_full$logFC.NM.m.Mic <- apply(mic_nm_full[,col], 1, function(x) log2(mean(x[c("NM.1","NM.2","NM.3")])) - log2(mean(x[c("MIC.1","MIC.2","MIC.3")])))


mic_nm_full$pval.NM.vs.MIC.ttest <- as.numeric(mic_nm_full$pval.NM.vs.MIC.ttest)
mic_nm_full$GSEA.Mic.m.NM <- apply(mic_nm_full[,c("pval.NM.vs.MIC.ttest", "logFC.Mic.m.NM")],1,function(x) -log10(x[c("pval.NM.vs.MIC.ttest")])/sign(x[c("logFC.Mic.m.NM")]))
mic_nm_full$GSEA.NM.m.Mic <- apply(mic_nm_full[,c("pval.NM.vs.MIC.ttest","logFC.NM.m.Mic")],1,function(x) -log10(x[c("pval.NM.vs.MIC.ttest")])/sign(x[c("logFC.NM.m.Mic")]))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

mic_nm_full <- completeFun(mic_nm_full, c("GSEA.Mic.m.NM", "GSEA.NM.m.Mic"))


NCBItoGENEID <- read.delim(file = "Data/Anno/NCBItoGENEPAGE.txt", header = F)
GENEIDtoGENEPAGE <- read.delim(file = "Data/Anno/GENEIDtoGENEPAGE.txt", header = F)
NCBItoGENEID <- subset(NCBItoGENEID, select = c("V2","V3"))
GENEIDtoGENEPAGE <- subset(GENEIDtoGENEPAGE, select = c("V1","V3"))

mic_nm_full$Protein.ID<- gsub("\\..*","",mic_nm_full$Protein.ID)

colnames(NCBItoGENEID) <- c("Protein.ID","GENEID")
colnames(GENEIDtoGENEPAGE) <- c("GENEPAGE","GENEID")
mic_nm_full <- merge(mic_nm_full, NCBItoGENEID, by = "Protein.ID" )
mic_nm_full <- merge(mic_nm_full, GENEIDtoGENEPAGE, by = "GENEID")


#######Annotations########
#Pathway annotation Files
KEGGt2g <- read.csv(file = "Data/Anno/KEGGt2g.csv", row.names = 1)
KEGGt2n <- read.csv(file = "Data/Anno/KEGGt2n.csv", row.names = 1)
GOt2g <- read.csv(file = "Data/Anno/GOt2g.csv", row.names = 1)
GOt2n <- read.csv(file = "Data/Anno/GOt2n.csv", row.names = 1)
FUNCt2g <- read.csv(file = "Data/Anno/FUNCt2g.csv", row.names = 1)
FUNCt2n <- read.csv(file = "Data/Anno/FUNCt2n.csv", row.names = 1)
FUNCver2_t2g <- read.csv(file = "Data/Anno/FUNCver2_t2g.csv", row.names = 1)
#######PCA#######
#All proteins
mic_nm_samples <- t(mic_nm_full[,c("MIC.1","MIC.2","MIC.3","NM.1","NM.2","NM.3","Protein.ID")])
mic_nm_samples <- as.data.frame(mic_nm_samples)
colnames(mic_nm_samples) <- mic_nm_samples[rownames(mic_nm_samples) == "Protein.ID",]
mic_nm_samples <- mic_nm_samples[1:6,] 


sapply(1:length(mic_nm_samples), function(i) {
  mic_nm_samples[, i] <<- as.numeric(as.character(mic_nm_samples[, i]))
});

######PCA on every protein found in the dataset
pca_tot <- prcomp(mic_nm_samples, center = TRUE, scale = TRUE)
#PCA plot
pca.plot <-  as.data.frame(as.matrix(pca_tot$x))
pca.plot$stages <- rownames(pca.plot)
pca.plot$spatial_lineage <- gsub("\\..*","",pca.plot$stages)
x <- summary(pca_tot)
per_var_exp <- as.data.frame(x$importance)

#Order legend correctly (specify factor levels)
pca.plot$spatial_lineage <- factor(pca.plot$spatial_lineage, levels = c("MIC",
                                                                "NM"))
levels(pca.plot$spatial_lineage) <- c("Micromere","Non-Micromere")
pca.plot$spatial_lineage

library(scales)
micnm_PCA <- ggplot(pca.plot,aes(x=PC1,y=PC2,col=spatial_lineage, label = stages))+
  geom_mark_ellipse(aes(x=PC1,y=PC2,fill = spatial_lineage), show.legend = FALSE, inherit.aes = FALSE) +
  geom_point()+ 
  geom_text(aes(label = stages), color = "black", hjust=1.1, vjust=1.1, size = 4.5) +
  xlab(paste0("PC1 (",percent(per_var_exp[2,c("PC1")], accuracy = 0.1), ")")) +
  ylab(paste0("PC2 (",percent(per_var_exp[2,c("PC2")], accuracy = 0.1), ")")) +
  ggtitle("PCA of Proteomics Micromere vs Non-Micromere") +
  xlim(-100,70)+
  ylim(-80,60)+
  theme_light() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 14),
        legend.text=element_text(size = 14)
  )
micnm_PCA

#Output Plot:
ggsave(filename = "micnm_PCA.png", plot = micnm_PCA,
       path = "Output/micnm_supplementary_figures",
       width = 14,
       height = 11,
       units = "in",
       dpi = 600)


#######Heatmap#######
#Heatmaps - with most significant genes (ANOVA)
micnm_heatmap = subset(mic_nm_full, select = c("MIC.1","MIC.2","MIC.3",
                                              "NM.1","NM.2","NM.3",
                                              "Protein.ID", "pval.NM.vs.MIC.ttest"))
micnm_heatmap_significant <- subset(micnm_heatmap, as.numeric(micnm_heatmap$pval.NM.vs.MIC.ttest) < 0.05)
rownames(micnm_heatmap_significant) <- micnm_heatmap_significant$Protein.ID 
micnm_heatmap_significant$Protein.ID <- NULL
micnm_heatmap_significant$pval.NM.vs.MIC.ttest <- NULL



micnm_names <- gsub("\\..*","",colnames(micnm_heatmap_significant))

level = factor(micnm_names, levels = c("MIC","NM"))
levels(level) <- c("Micromere","Non-Micromere")
level

#IDS in the color pallete of ggplot2
#library(scales)
#show_col(hue_pal()(2))
colors_gg <- list('lineage' = c("Micromere" = "#F8766D",
                                  "Non-Micromere" = "#00BFC4"))
                                  

ha = HeatmapAnnotation(lineage = level, annotation_name_side = "left", col = colors_gg)

######No Z-score (Preferred)
micnm_heatmap_significant <- as.matrix(micnm_heatmap_significant)
ht_list = Heatmap(micnm_heatmap_significant,
                  name = "scale",
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list)



######Z-score version
png("Output/micnm_supplementary_figures/micnm_heatmap.png",width=6.25,height=4.25,units="in",res=1200)
ht_list_zscore = Heatmap(t(scale(t(micnm_heatmap_significant))), #from complexHeatmap github
                  name = "scale",
                  top_annotation = ha,
                  row_title = NULL,
                  show_row_names = FALSE,
                  cluster_columns = FALSE,
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "complete")
draw(ht_list_zscore)
dev.off()


#######Overrepresentation Analysis: Significant Proteins########
#Analogous to the mfuzz analysis
micnm_significant_prot <- subset(mic_nm_full, as.numeric(mic_nm_full$pval.NM.vs.MIC.ttest) < 0.05)

#ORA analysis, make a list with significant proteins
micnm_ora <- list()
micnm_ora[["Micromere Up-regulated"]] <- micnm_significant_prot[micnm_significant_prot$logFC.Mic.m.NM > 0,]$GENEPAGE
micnm_ora[["Micromere Down-regulated"]] <- micnm_significant_prot[micnm_significant_prot$logFC.Mic.m.NM < 0,]$GENEPAGE

#Previous p-value cutoff 0.2
#KEGG
KEGG_clust <- compareCluster(geneCluster = micnm_ora,  fun = "enricher", TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, pvalueCutoff = 0.2)
png("Output/micnm_supplementary_figures/micnm_kegg_ora.png",width=10.25,height=6.25,units="in",res=1200)
dotplot(KEGG_clust, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")
dev.off()

#GO
GO_clust <- compareCluster(geneCluster = micnm_ora,  fun = "enricher", TERM2GENE = GOt2g, TERM2NAME = GOt2n, pvalueCutoff = 0.2)
#method to remove NA, mutate result with complete cases argument
GO_clust@compareClusterResult <- GO_clust@compareClusterResult[complete.cases(GO_clust@compareClusterResult), ]
png("Output/micnm_supplementary_figures/micnm_go_ora.png",width=10.25,height=6.25,units="in",res=1200)
dotplot(GO_clust, showCategory = 10) + ggtitle("GO Pathway Enrichment")
dev.off()

#FUNC
FUNC_clust <- compareCluster(geneCluster = micnm_ora,  fun = "enricher", TERM2GENE = FUNCt2g, pvalueCutoff = 0.2)
png("Output/micnm_supplementary_figures/micnm_func_ora.png",width=10.25,height=6.25,units="in",res=1200)
dotplot(FUNC_clust, showCategory = 10) + ggtitle("FUNC Pathway Enrichment")
dev.off()

#######GSEA Analysis#######
#make GSEA genelist
MIC.m.NM <-  mic_nm_full$GSEA.Mic.m.NM
names(MIC.m.NM) <- mic_nm_full$GENEPAGE
MIC.m.NM <- sort(MIC.m.NM, decreasing = TRUE)

##Functional list (generalized)
gsea.MIC = GSEA(MIC.m.NM, pvalueCutoff = 0.1,  TERM2GENE = FUNCver2_t2g, verbose=FALSE)
head(gsea.MIC)
png("Output/micnm_supplementary_figures/micnm_funcver2_gsea.png",width=10.25,height=8.25,units="in",res=1200)
ridgeplot(gsea.MIC) + ggtitle("FUNC: GSEA enrichment in Micromeres")
dev.off()
gsea.FUNCv1allmic <- gsea.MIC@result

##Functional list (specific)
gsea.MIC = GSEA(MIC.m.NM, pvalueCutoff = 0.1,  TERM2GENE = FUNCt2g, verbose=FALSE)
head(gsea.MIC)
png("Output/micnm_supplementary_figures/micnm_func_gsea.png",width=10.25,height=8.25,units="in",res=1200)
ridgeplot(gsea.MIC) + ggtitle("FUNC: GSEA enrichment in Micromeres")
dev.off()
gsea.FUNCv2allmic <- gsea.MIC@result

##GO list
gsea.GO.mic <- GSEA(MIC.m.NM, pvalueCutoff = 0.1,  TERM2GENE = GOt2g, TERM2NAME = GOt2n, verbose=FALSE)
gsea.GO.mic@result<- gsea.GO.mic@result[complete.cases(gsea.GO.mic@result), ]
head(gsea.GO.mic)
png("Output/micnm_supplementary_figures/micnm_go_gsea.png",width=10.25,height=8.25,units="in",res=1200)
ridgeplot(gsea.GO.mic) + ggtitle("GO: GSEA enrichment in Micromeres")
dev.off()
gsea.GO.allmic <- gsea.GO.mic@result

##KEGG list
gsea.KEGG.mic <-  GSEA(MIC.m.NM, pvalueCutoff = 0.1,  TERM2GENE = KEGGt2g, TERM2NAME = KEGGt2n, verbose=FALSE)
head(gsea.KEGG.mic)
png("Output/micnm_supplementary_figures/micnm_kegg_gsea.png",width=10.25,height=8.25,units="in",res=1200)
ridgeplot(gsea.KEGG.mic) + ggtitle("KEGG: GSEA enrichment in Micromeres")
dev.off()
gsea.KEGG.allmic <- gsea.KEGG.mic@result







####Extra Stuff####
#Example of 2D plot to show enrichment (Not used for manuscript)
gseaGO_full <- gsea.GO.mic@result
gseaGO_full$log10pval <- -log10(gseaGO_full$p.adjust)
gseaGO_plot <- ggplot(gseaGO_full, aes(x = enrichmentScore, y = log10pval, color = enrichmentScore, size = setSize)) +
  geom_point()

gseaGO_plot +   geom_text_repel(aes(label = Description), size = 4)

