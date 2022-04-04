library(pathview)

#Initialize datasets from "global_metab_maps.R" and "timecourse_all.R"
timecourse_metab <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Data/timecourse_metabolomics_global.csv", row.names = 1)
timecourse_prot <- read.csv(file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/Manuscript R files/Data/timecourse_proteomics_global.csv", row.names = 1)

#Remove ANOVA column
timecourse_prot <- timecourse_prot[,1:8]



######Z-score both datasets {method #1}###########
timecourse_prot_zscore  <- as.data.frame(t(scale(t(timecourse_prot))))
timecourse_metab_zscore  <- as.data.frame(t(scale(t(timecourse_metab))))

#Pathview plots: (by z-score)
pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         out.suffix = "glycolysis.zscore",
         pathway.id = "spu00010",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         out.suffix = "TCA.zscore",
         pathway.id = "spu00020",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         out.suffix = "OXPHOS.zscore",
         pathway.id = "spu00190",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         out.suffix = "FAbio.zscore",
         pathway.id = "spu00061",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         pathway.id = "spu00030",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         pathway.id = "spu00061",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_zscore,
         cpd.data =  timecourse_metab_zscore,
         pathway.id = "spu00290",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

#########By max value {method #2}#########
min_max_scaling <- function(vector){
  norm <- (vector - min(vector))/(max(vector) - min(vector))
  return(norm)
} 

#Test
test1 <- c(0,1,10,5,5,9,12,0.1)
min_max_scaling(test1)
  
timecourse_prot_maxscale <- as.data.frame(t(apply(timecourse_prot, 1, min_max_scaling)))
timecourse_metab_maxscale <- as.data.frame(t(apply(timecourse_metab, 1, min_max_scaling)))

pathview(gene.data  = timecourse_prot_maxscale,
         cpd.data =  timecourse_metab_maxscale,
         out.suffix = "min.max.3",
         pathway.id = "spu00010",
         species    = 'spu',
         both.dirs = list(gene=FALSE, cpd=FALSE),
         low = list(gene = "green", cpd = "blue"),
         mid = list(gene = "gray", cpd = "gray"),
         high = list(gene = "red", cpd = "yellow"),
         limit = list(gene = 1, cpd = 1), #"provide limit values when putting in color format."
         #discrete = list(gene = T, cpd = F),
         match.data = TRUE)

pathview(gene.data  = timecourse_prot_maxscale,
         cpd.data =  timecourse_metab_maxscale,
         out.suffix = "min.max.2",
         pathway.id = "spu00010",
         species    = 'spu',
         limit = list(gene = 1, cpd = 1), #"provide limit values when putting in color format."
         discrete = list(gene = T, cpd = F),
         match.data = TRUE)


#########Normalize by Egg {method #3}#########
egg_scaling <- function(vector){
  norm <- (vector - vector[1])/(max(vector) - min(vector))
  return(norm)
} 

#Test
test1 <- c(7,1,10,5,5,9,12,0.1)
egg_scaling(test1)

timecourse_prot_eggscale <- as.data.frame(t(apply(timecourse_prot, 1, egg_scaling)))
timecourse_metab_eggscale <- as.data.frame(t(apply(timecourse_metab, 1, egg_scaling)))
write.csv(timecourse_prot_eggscale, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/TIMECOURSE/timecourse.prot.eggscale.csv")
write.csv(timecourse_metab_eggscale, file = "C:/Users/shaks/OneDrive/Documents/R/Proteomics Summer/App/TIMECOURSE/timecourse.metab.eggscale.csv")





pathview(gene.data  = timecourse_prot_eggscale,
         cpd.data =  timecourse_metab_eggscale,
         out.suffix = "glycolysis.eggscale",
         pathway.id = "spu00010",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_eggscale,
         cpd.data =  timecourse_metab_eggscale,
         out.suffix = "tca.eggscale",
         pathway.id = "spu00020",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_eggscale,
         cpd.data =  timecourse_metab_eggscale,
         out.suffix = "fa.eggscale",
         pathway.id = "spu00061",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_eggscale,
         cpd.data =  timecourse_metab_eggscale,
         out.suffix = "oxphos.eggscale",
         pathway.id = "spu00190",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)

pathview(gene.data  = timecourse_prot_eggscale,
         cpd.data =  timecourse_metab_eggscale,
         out.suffix = "FAbio.eggscale",
         pathway.id = "spu00061",
         species    = 'spu',
         #limit = list(gene = 2, cpd = 2), #"provide limit values when putting in color format."
         discrete = list(gene = TRUE, cpd = TRUE),
         match.data = FALSE)
