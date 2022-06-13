library(pathview)

#Initialize datasets from "global_metab_maps.R" and "timecourse_all.R"
timecourse_metab <- read.csv(file = "Data/App_data/timecourse_metabolomics_global.csv", row.names = 1)
timecourse_prot <- read.csv(file = "Data/App_data/timecourse_proteomics_global.csv", row.names = 1)

#Remove ANOVA column
timecourse_prot <- timecourse_prot[,1:8]

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
write.csv(timecourse_prot_eggscale, file = "Output/data_files_for_app/timecourse.prot.eggscale.csv")
write.csv(timecourse_metab_eggscale, file = "Output/data_files_for_app/timecourse.metab.eggscale.csv")
