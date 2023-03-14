# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(fuzzySim)
library(modEvA)
library(maps)
library(ENMeval)


# CREATE OUTPUT FOLDERS ####

dir.create("../pred_rastersBIAS_CLAMPING/native", recursive = TRUE)
dir.create("../pred_rastersBIAS_CLAMPING/EU/invasive", recursive = TRUE)
dir.create("../pred_JPGsBIAS_CLAMPING", recursive = TRUE)


# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# GET MODEL PREDICTION RASTERS FOR EACH SPECIES ####

for (species in species_list) {
  
  #pa_nat$glm_p <- predict(mod_glm, pa_nat, type = "response")
  #pred_bart <- predict_bart_df(mod_bart, pa_nat, quantiles = c(0.025, 0.975))
  #head(pred_bart)
  #pred_bart$uncert <- pred_bart[ , 3] - pred_bart[ , 2]
  #names(pred_bart) <- c("bart_p", "bart_q0025", "bart_q0975", "bart_uncert")
  #pa_nat <- data.frame(pa_nat, pred_bart)
  #head(pa_nat)
  
  #pa_nat_sv <- vect(pa_nat, geom = c("lon", "lat"))
  #par(mfrow = c(3, 1))
  #plot(pa_nat_sv, "presence", cex = 0.1, col = c(NA, "black"))
  #plot(pa_nat_sv, "glm_p", cex = 0.1)
  #plot(pa_nat_sv, "bart_p", cex = 0.1)
  
  gc()
  
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  climgrids_nat$bias <- (climgrids_nat[[1]]*0)
  
  climgrids_inv <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.invasive.tif"))

  
###ADJUST climgrids_inv FOR CLAMPING###
  sp.id <- read.csv("../../0_Nf_modeling/allspecies_samplesizes_v3.csv",header=T)[,c(1,4,6)]
  j<-species
  f.occ.nat <- paste0("../../0_Nf_modeling/accessible-areas-v3/",j,"_native_range.csv")
  sp.occ.nat <- read.csv(f.occ.nat,header=T)[,-1]
  
  #spatRaster issues
  temp <- stack(climgrids_inv)
  temp[[1]]
  temp[[2]]
  names(temp) <- names(sp.occ.nat[c(3:4)])
  
  climgrids_inv <- clamp.vars(temp, sp.occ.nat[c(3:4)], left = NULL, right = NULL, categoricals = NULL)
  climgrids_inv$bias <- (climgrids_inv[[1]]*0)
  
  mod_glm <- readRDS(paste0("../model_objectsBIAS/", species, "_mod_glm.rds"))
  mod_bart <- readRDS(paste0("../model_objectsBIAS/", species, "_mod_bart.rds"))
  
  message("\n", species, "\npredicting probability on raster maps...")
  var_names <- names(mod_glm$coefficients)[-1]
  names(climgrids_nat) <- names(climgrids_inv) <- var_names
  
  message(" - GLM...")
  glm_p_nat <- predict(climgrids_nat, mod_glm, type = "response")
  glm_p_inv <- predict(climgrids_inv, mod_glm, type = "response")
  message(" - BART...")
  bart_p_nat <- predict(mod_bart, raster::stack(climgrids_nat))
  bart_p_nat <- rast(bart_p_nat)
  gc()
  bart_p_inv <- predict(mod_bart, raster::stack(climgrids_inv))
  bart_p_inv <- rast(bart_p_inv)
  
  message("converting to favourability...")
  glm_f_nat <- fuzzySim::Fav(pred = glm_p_nat, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_nat <- fuzzySim::Fav(pred = bart_p_nat, sample.preval = modEvA::prevalence(model = mod_bart))
  glm_f_inv <- fuzzySim::Fav(pred = glm_p_inv, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_inv <- fuzzySim::Fav(pred = bart_p_inv, sample.preval = modEvA::prevalence(model = mod_bart))
  
  message("saving prediction raster maps...")
  writeRaster(glm_p_nat, paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_glm.tif"),overwrite=TRUE)
  writeRaster(bart_p_nat, paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_bart.tif"),overwrite=TRUE)
  writeRaster(glm_f_nat, paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_f_glm.tif"),overwrite=TRUE)
  writeRaster(bart_f_nat, paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_f_bart.tif"),overwrite=TRUE)
  
  writeRaster(glm_p_inv, paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_glm.tif"),overwrite=TRUE)
  writeRaster(bart_p_inv, paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_bart.tif"),overwrite=TRUE)
  writeRaster(glm_f_inv, paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_f_glm.tif"),overwrite=TRUE)
  writeRaster(bart_f_inv, paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_f_bart.tif"),overwrite=TRUE)
  
  gc()
}  # end predicting


# GET PREDICTION RASTER MAP FILE NAMES ####

pred_files_nat <- list.files("../pred_rastersBIAS_CLAMPING/native", full.names = TRUE)
pred_files_nat <- pred_files_nat[grep("_f_", pred_files_nat)]
pred_files_nat

pred_files_inv <- list.files("../pred_rastersBIAS_CLAMPING/EU/invasive", full.names = TRUE)
pred_files_inv <- pred_files_inv[grep("_f_", pred_files_inv)]
pred_files_inv


# IMPORT EACH RASTER AND SAVE AS JPG ####

for (f in c(pred_files_nat, pred_files_inv)) {
  rst <- rast(f)
  filename <- basename(tools::file_path_sans_ext(f))
  species <- sapply(strsplit(filename, "_"), `[`, 1)
  jpeg(paste0("../pred_JPGsBIAS_CLAMPING/", filename, ".jpg"), width = 800, height = 430)
  maps::map("world", col = "grey", mar = c(1, 1, 2, 6))
  plot(rst, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(basename(tools::file_path_sans_ext(f)))
  dev.off()
}
