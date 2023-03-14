# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(fuzzySim)
library(modEvA)
library(maps)
library(foreach)
library(doParallel)

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)


# CREATE OUTPUT FOLDERS ####

dir.create("../pred_rastersBIAS_BIOTIC_HABITAT/native", recursive = TRUE)
dir.create("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive", recursive = TRUE)
dir.create("../pred_JPGsBIAS_BIOTIC_HABITAT", recursive = TRUE)


# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# GET MODEL PREDICTION RASTERS FOR EACH SPECIES ####

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','maps'))  %dopar% {
  
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
  
  #import biotic grids
  bioticgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_bioticgrids/", species, "/", species, "_bioticGridsPCA.native.tif"))
  bioticgrids_inv <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_bioticgrids/", species, "/", species, "_bioticGridsPCA.invasive.tif"))
  
  
  #import habitat grids
  habitatgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_habitatgrids/", species, "/", species, "_habitatGridsPCA.native.tif"))
  habitatgrids_inv <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_habitatgrids/", species, "/", species, "_habitatGridsPCA.invasive.tif"))
  
  
  #import climgrids:
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  climgrids_inv <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.invasive.tif"))
  
  #native range grids
  mask_nat <- raster(climgrids_nat[[1]]) + raster(climgrids_nat[[2]]) + raster(bioticgrids_nat[[1]]) + raster(bioticgrids_nat[[2]]) + raster(habitatgrids_nat[[1]]) + raster(habitatgrids_nat[[2]]) 
  mask_nat <- mask_nat*0
  PC1 <- mask_nat+raster(climgrids_nat[[1]])
  PC2 <- mask_nat+raster(climgrids_nat[[2]])
  PCbiotic1 <- mask_nat+raster(bioticgrids_nat[[1]])
  PCbiotic2 <- mask_nat+raster(bioticgrids_nat[[2]])
  PChabitat1 <- mask_nat+raster(habitatgrids_nat[[1]])
  PChabitat2 <- mask_nat+raster(habitatgrids_nat[[2]])
    predictorgrids_nat <- stack(PC1,PC2,PCbiotic1,PCbiotic2,PChabitat1,PChabitat2)
      predictorgrids_nat$bias <- predictorgrids_nat[[1]]*0
      plot(predictorgrids_nat)
  
  #invasive range grids
  mask_inv <- raster(climgrids_inv[[1]]) + raster(climgrids_inv[[2]]) + raster(bioticgrids_inv[[1]]) + raster(bioticgrids_inv[[2]]) + raster(habitatgrids_inv[[1]]) + raster(habitatgrids_inv[[2]]) 
  mask_inv <- mask_inv*0
  PC1 <- mask_inv+raster(climgrids_inv[[1]])
  PC2 <- mask_inv+raster(climgrids_inv[[2]])
  PCbiotic1 <- mask_inv+raster(bioticgrids_inv[[1]])
  PCbiotic2 <- mask_inv+raster(bioticgrids_inv[[2]])
  PChabitat1 <- mask_inv+raster(habitatgrids_inv[[1]])
  PChabitat2 <- mask_inv+raster(habitatgrids_inv[[2]])
    predictorgrids_inv <- stack(PC1,PC2,PCbiotic1,PCbiotic2,PChabitat1,PChabitat2)
      predictorgrids_inv$bias <- predictorgrids_inv[[1]]*0
      plot(predictorgrids_inv)
  
  mod_glm <- readRDS(paste0("../model_objectsBIAS_BIOTIC_HABITAT/", species, "_mod_glm.rds"))
  mod_bart <- readRDS(paste0("../model_objectsBIAS_BIOTIC_HABITAT/", species, "_mod_bart.rds"))
  
  message("\n", species, "\npredicting probability on raster maps...")
  var_names <- names(mod_glm$coefficients)[-1]
  names(predictorgrids_nat) <- names(predictorgrids_inv) <- var_names
  
  message(" - GLM...")
  glm_p_nat <- predict(predictorgrids_nat, mod_glm, type = "response")
  glm_p_inv <- predict(predictorgrids_inv, mod_glm, type = "response")
  message(" - BART...")
  bart_p_nat <- predict(mod_bart, raster::stack(predictorgrids_nat))
  bart_p_nat <- rast(bart_p_nat)
  gc()
  bart_p_inv <- predict(mod_bart, raster::stack(predictorgrids_inv))
  bart_p_inv <- rast(bart_p_inv)
  
  message("converting to favourability...")
  glm_f_nat <- fuzzySim::Fav(pred = glm_p_nat, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_nat <- fuzzySim::Fav(pred = bart_p_nat, sample.preval = modEvA::prevalence(model = mod_bart))
  glm_f_inv <- fuzzySim::Fav(pred = glm_p_inv, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_inv <- fuzzySim::Fav(pred = bart_p_inv, sample.preval = modEvA::prevalence(model = mod_bart))
  
  message("saving prediction raster maps...")
  writeRaster(glm_p_nat, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/native/", species, "_nat_p_glm.tif"),overwrite=TRUE)
  writeRaster(bart_p_nat, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/native/", species, "_nat_p_bart.tif"),overwrite=TRUE)
  writeRaster(glm_f_nat, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/native/", species, "_nat_f_glm.tif"),overwrite=TRUE)
  writeRaster(bart_f_nat, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/native/", species, "_nat_f_bart.tif"),overwrite=TRUE)
  
  writeRaster(glm_p_inv, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive/", species, "_inv_p_glm.tif"),overwrite=TRUE)
  writeRaster(bart_p_inv, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive/", species, "_inv_p_bart.tif"),overwrite=TRUE)
  writeRaster(glm_f_inv, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive/", species, "_inv_f_glm.tif"),overwrite=TRUE)
  writeRaster(bart_f_inv, paste0("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive/", species, "_inv_f_bart.tif"),overwrite=TRUE)
  
  gc()
}  # end predicting


# GET PREDICTION RASTER MAP FILE NAMES ####

pred_files_nat <- list.files("../pred_rastersBIAS_BIOTIC_HABITAT/native", full.names = TRUE)
pred_files_nat <- pred_files_nat[grep("_f_", pred_files_nat)]
pred_files_nat

pred_files_inv <- list.files("../pred_rastersBIAS_BIOTIC_HABITAT/EU/invasive", full.names = TRUE)
pred_files_inv <- pred_files_inv[grep("_f_", pred_files_inv)]
pred_files_inv


# IMPORT EACH RASTER AND SAVE AS JPG ####

for (f in c(pred_files_nat, pred_files_inv)) {
  rst <- rast(f)
  filename <- basename(tools::file_path_sans_ext(f))
  species <- sapply(strsplit(filename, "_"), `[`, 1)
  jpeg(paste0("../pred_JPGsBIAS_BIOTIC_HABITAT/", filename, ".jpg"), width = 800, height = 430)
  maps::map("world", col = "grey", mar = c(1, 1, 2, 6))
  plot(rst, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(basename(tools::file_path_sans_ext(f)))
  dev.off()
}
