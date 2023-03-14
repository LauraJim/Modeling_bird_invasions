# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####
library(raster)
library(embarcadero)
library(blockCV)
library(modEvA)
library(kuenm)
library(foreach)
library(doParallel)


# GET TARGET SPECIES ####
occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the BART predict function to work with data frames instead of raster layers

##########################################################################################################################################
###CLIMATE 5E#############################################################################################################################
##########################################################################################################################################
nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  #get the threshold for the threshold dependent evaluation metrics
  native_range <- raster(paste0("../../../0_Analyses/0_Nf_modeling/ResultsClimate/",species,"_suitability_map_nat.tif"))
  native_range_occs <- read.csv(paste0("../../../0_Analyses/0_Nf_modeling/occurrences-v3/",species,"_native.csv"))
  native_range_suitabilities <- raster::extract(native_range,native_range_occs[c(2:3)])
  
  pred_nat <- read.csv(paste0("../../../0_Analyses/1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"))[c(2:3,7)]
  pred_nat$fne_suit <- raster::extract(native_range,pred_nat[c(1:2)]) 
  
  thresh <- quantile(native_range_suitabilities, probs = 0.05, na.rm = TRUE)  # pred for minimum 5% presence in the training data

  proc <- kuenm_proc(occ.test = native_range_occs[c(2:3)], model = native_range, threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary
  threshmeas <- threshMeasures(obs = pred_nat$presence, pred = pred_nat$fne_suit, measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
  tempeval <- data.frame(data.frame(t(proc)),data.frame(t(threshmeas)))
  tempeval$preds <- "clim"
  tempeval$species <- species
  tempeval$range <- "native"
  tempeval$threshold <- "5"
  write.csv(tempeval, paste0("../../../0_Analyses/0_Nf_modeling/eval_metrics/NATIVE/BIAS/5E/fne_native_", species, ".csv"), row.names = FALSE)
}
rm(cl)

##########################################################################################################################################
###CLIMATE 2.5E#############################################################################################################################
##########################################################################################################################################
nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  #get the threshold for the threshold dependent evaluation metrics
  native_range <- raster(paste0("../../../0_Analyses/0_Nf_modeling/ResultsClimate/",species,"_suitability_map_nat.tif"))
  native_range_occs <- read.csv(paste0("../../../0_Analyses/0_Nf_modeling/occurrences-v3/",species,"_native.csv"))
  native_range_suitabilities <- raster::extract(native_range,native_range_occs[c(2:3)])
  
  pred_nat <- read.csv(paste0("../../../0_Analyses/1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"))[c(2:3,7)]
  pred_nat$fne_suit <- raster::extract(native_range,pred_nat[c(1:2)]) 
  
  thresh <- quantile(native_range_suitabilities, probs = 0.025, na.rm = TRUE)  # pred for minimum 5% presence in the training data
  
  proc <- kuenm_proc(occ.test = native_range_occs[c(2:3)], model = native_range, threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary
  threshmeas <- threshMeasures(obs = pred_nat$presence, pred = pred_nat$fne_suit, measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
  tempeval <- data.frame(data.frame(t(proc)),data.frame(t(threshmeas)))
  tempeval$preds <- "clim"
  tempeval$species <- species
  tempeval$range <- "native"
  tempeval$threshold <- "2.5"
  write.csv(tempeval, paste0("../../../0_Analyses/0_Nf_modeling/eval_metrics/NATIVE/BIAS/2.5E/fne_native_", species, ".csv"), row.names = FALSE)
}
rm(cl)

##########################################################################################################################################
###CLIMATE BIOTIC HABITAT 5E##############################################################################################################
##########################################################################################################################################
nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  #get the threshold for the threshold dependent evaluation metrics
  native_range <- raster(paste0("../../../0_Analyses/0_Nf_modeling/ResultsClimBioticHabitat/",species,"_suitability_map_nat.tif"))
  native_range_occs <- read.csv(paste0("../../../0_Analyses/0_Nf_modeling/occurrencesClimBioticHabitat/",species,"_native.csv"))
  native_range_suitabilities <- raster::extract(native_range,native_range_occs[c(2:3)])
  
  pred_nat <- read.csv(paste0("../../../0_Analyses/1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/pred_CSVsBIAS_BIOTIC_HABITAT/", species, "_pred_crossval.csv"))[c(2:3,11)]
  pred_nat$fne_suit <- raster::extract(native_range,pred_nat[c(1:2)]) 
  
  thresh <- quantile(native_range_suitabilities, probs = 0.05, na.rm = TRUE)  # pred for minimum 5% presence in the training data
  
  proc <- kuenm_proc(occ.test = native_range_occs[c(2:3)], model = native_range, threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary
  threshmeas <- threshMeasures(obs = pred_nat$presence, pred = pred_nat$fne_suit, measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
  tempeval <- data.frame(data.frame(t(proc)),data.frame(t(threshmeas)))
  tempeval$preds <- "clim_hb"
  tempeval$species <- species
  tempeval$range <- "native"
  tempeval$threshold <- "5"
  write.csv(tempeval, paste0("../../../0_Analyses/0_Nf_modeling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/fne_native_", species, ".csv"), row.names = FALSE)
}
rm(cl)

##########################################################################################################################################
###CLIMATE BIOTIC HABITAT 2.5E##############################################################################################################
##########################################################################################################################################
nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  #get the threshold for the threshold dependent evaluation metrics
  native_range <- raster(paste0("../../../0_Analyses/0_Nf_modeling/ResultsClimBioticHabitat/",species,"_suitability_map_nat.tif"))
  native_range_occs <- read.csv(paste0("../../../0_Analyses/0_Nf_modeling/occurrencesClimBioticHabitat/",species,"_native.csv"))
  native_range_suitabilities <- raster::extract(native_range,native_range_occs[c(2:3)])
  
  pred_nat <- read.csv(paste0("../../../0_Analyses/1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/pred_CSVsBIAS_BIOTIC_HABITAT/", species, "_pred_crossval.csv"))[c(2:3,11)]
  pred_nat$fne_suit <- raster::extract(native_range,pred_nat[c(1:2)]) 
  
  thresh <- quantile(native_range_suitabilities, probs = 0.025, na.rm = TRUE)  # pred for minimum 5% presence in the training data
  
  proc <- kuenm_proc(occ.test = native_range_occs[c(2:3)], model = native_range, threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary
  threshmeas <- threshMeasures(obs = pred_nat$presence, pred = pred_nat$fne_suit, measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
  tempeval <- data.frame(data.frame(t(proc)),data.frame(t(threshmeas)))
  tempeval$preds <- "clim_hb"
  tempeval$species <- species
  tempeval$range <- "native"
  tempeval$threshold <- "2.5"
  write.csv(tempeval, paste0("../../../0_Analyses/0_Nf_modeling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/2.5E/fne_native_", species, ".csv"), row.names = FALSE)
}
rm(cl)
