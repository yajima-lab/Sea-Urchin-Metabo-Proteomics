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