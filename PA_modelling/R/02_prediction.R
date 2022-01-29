# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(fuzzySim)
library(modEvA)


# CREATE OUTPUT FOLDERS ####

dir.create("../pred_rasters/native", recursive = TRUE)
dir.create("../pred_rasters/invasive", recursive = TRUE)


# GET TARGET SPECIES ####

occ_files <- list.files("../../Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

# "psikra" crashes my RStudio when predicting on raster, so leave it out for now:
species_list <- species_list[-grep("psikra", species_list)]
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
  
  climgrids_nat <- rast(paste0("../../occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  climgrids_inv <- rast(paste0("../../occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.invasive.tif"))
  mod_glm <- readRDS(paste0("../model_objects/", species, "_mod_glm.rds"))
  mod_bart <- readRDS(paste0("../model_objects/", species, "_mod_bart.rds"))
  
  message("\n", species, "\npredicting probability on raster maps...")
  var_names <- names(mod_glm$coefficients)[-1]
  names(climgrids_nat) <- names(climgrids_inv) <- var_names
  
  message(" - GLM...")
  glm_p_nat <- predict(climgrids_nat, mod_glm, type = "response")
  glm_p_inv <- predict(climgrids_inv, mod_glm, type = "response")
  message(" - BART...")
  bart_p_nat <- predict(mod_bart, raster::stack(climgrids_nat))
  bart_p_nat <- rast(bart_p_nat)
  bart_p_inv <- predict(mod_bart, raster::stack(climgrids_inv))
  bart_p_inv <- rast(bart_p_inv)
  
  message("converting to favourability...")
  glm_f_nat <- fuzzySim::Fav(pred = glm_p_nat, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_nat <- fuzzySim::Fav(pred = bart_p_nat, sample.preval = modEvA::prevalence(model = mod_bart))
  glm_f_inv <- fuzzySim::Fav(pred = glm_p_inv, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_inv <- fuzzySim::Fav(pred = bart_p_inv, sample.preval = modEvA::prevalence(model = mod_bart))
  
  message("saving prediction raster maps...")
  writeRaster(glm_p_nat, paste0("../pred_rasters/native/", species, "_nat_p_glm.tif"))
  writeRaster(bart_p_nat, paste0("../pred_rasters/native/", species, "_nat_p_bart.tif"))
  writeRaster(glm_f_nat, paste0("../pred_rasters/native/", species, "_nat_f_glm.tif"))
  writeRaster(bart_f_nat, paste0("../pred_rasters/native/", species, "_nat_f_bart.tif"))
  
  writeRaster(glm_p_inv, paste0("../pred_rasters/invasive/", species, "_inv_p_glm.tif"))
  writeRaster(bart_p_inv, paste0("../pred_rasters/invasive/", species, "_inv_p_bart.tif"))
  writeRaster(glm_f_inv, paste0("../pred_rasters/invasive/", species, "_inv_f_glm.tif"))
  writeRaster(bart_f_inv, paste0("../pred_rasters/invasive/", species, "_inv_f_bart.tif"))
  
  gc()
}  # end predicting

