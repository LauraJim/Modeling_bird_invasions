# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(spatialEco)
library(foreach)
library(doParallel)

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

# CREATE OUTPUT FOLDER ####

dir.create("../model_objectsBIAS_BIOTIC_HABITAT/", recursive = TRUE)


# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# COMPUTE PA (PRESENCE/ABSENCE) MODELS ####

foreach (species = species_list,.packages=c('terra','embarcadero','spatialEco'))  %dopar% {
  message("\n", species, "\npreparing data...")
  
  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  head(occ_nat)
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  head(acc_nat)
  
  
  #import biasgrid
  biasgrid_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_biasgrids/", species, "/", species, "_bias.tif"))

  #import biotic grids
  bioticgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_bioticgrids/", species, "/", species, "_bioticGridsPCA.native.tif"))

  #import habitat grids
  
  habitatgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_habitatgrids/", species, "/", species, "_habitatGridsPCA.native.tif"))

  #exctract bias, biotic and habitat values
  occ_nat$bias <- terra::extract(biasgrid_nat, occ_nat[c(2:3)])[2]$layer
  occ_nat$PCbiotic1 <- terra::extract(bioticgrids_nat, occ_nat[c(2:3)])$layer.1
  occ_nat$PCbiotic2 <- terra::extract(bioticgrids_nat, occ_nat[c(2:3)])$layer.2
  occ_nat$PChabitat1 <- terra::extract(habitatgrids_nat, occ_nat[c(2:3)])$layer.1
  occ_nat$PChabitat2 <- terra::extract(habitatgrids_nat, occ_nat[c(2:3)])$layer.2
  head(occ_nat)
  
  acc_nat$bias <- terra::extract(biasgrid_nat, acc_nat[c(2:3)])[2]$layer
  acc_nat$PCbiotic1 <- terra::extract(bioticgrids_nat, acc_nat[c(2:3)])$layer.1
  acc_nat$PCbiotic2 <- terra::extract(bioticgrids_nat, acc_nat[c(2:3)])$layer.2
  acc_nat$PChabitat1 <- terra::extract(habitatgrids_nat, acc_nat[c(2:3)])$layer.1
  acc_nat$PChabitat2 <- terra::extract(habitatgrids_nat, acc_nat[c(2:3)])$layer.2
  head(acc_nat)
  
  # import climgrids:
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  mask <- raster(climgrids_nat[[1]]) + raster(climgrids_nat[[2]]) + raster(bioticgrids_nat[[1]]) + raster(bioticgrids_nat[[2]]) + raster(habitatgrids_nat[[1]]) + raster(habitatgrids_nat[[2]]) + raster(biasgrid_nat)
  mask <- mask*0
  
  PC1 <- mask+raster(climgrids_nat[[1]])
  PC2 <- mask+raster(climgrids_nat[[2]])
  bias <- mask+raster(biasgrid_nat)
  PCbiotic1 <- mask+raster(bioticgrids_nat[[1]])
  PCbiotic2 <- mask+raster(bioticgrids_nat[[2]])
  PChabitat1 <- mask+raster(habitatgrids_nat[[1]])
  PChabitat2 <- mask+raster(habitatgrids_nat[[2]])
    
  predictorgrids_nat <- stack(PC1,PC2,bias,PCbiotic1,PCbiotic2,PChabitat1,PChabitat2)
  plot(predictorgrids_nat)
  
  #plot(climgrids_nat)
  
  # get presence (occupied) and absence (accessible not occupied) pixels:
  occ_nat_cells <- terra::cellFromXY(object = predictorgrids_nat, xy = as.matrix(occ_nat[ , c("lon", "lat")]))
  acc_nat_cells <- terra::cellFromXY(object = predictorgrids_nat, xy = as.matrix(acc_nat[ , c("lon", "lat")]))
  # length(occ_nat_cells)
  # length(acc_nat_cells)
  abs_nat <- acc_nat[!(acc_nat_cells %in% occ_nat_cells), ]
  occ_nat$presence <- 1
  abs_nat$presence <- 0
  pa_nat <- rbind(occ_nat, abs_nat)
  # head(pa_nat)
  # sum(pa_nat$presence == 1)
  # sum(pa_nat$presence == 0)
  
  
  # compute presence-absence models:
  
  message("computing models...")
  #names(pa_nat)
  var_names <- c(names(pa_nat)[grep("PC", names(pa_nat))],"bias")
  var_names
  
  form_glm <- as.formula(paste("presence ~", paste(var_names, collapse = "+")))
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  #summary(mod_glm)
  saveRDS(mod_glm, paste0("../model_objectsBIAS_BIOTIC_HABITAT/", species, "_mod_glm.rds"))
  
  set.seed(grep(species, species_list))
  mod_bart <- bart(x.train = pa_nat[ , var_names], y.train = pa_nat[ , "presence"], keeptrees = TRUE, verbose = FALSE)
  #summary(mod_bart)
  invisible(mod_bart$fit$state)
  saveRDS(mod_bart, paste0("../model_objectsBIAS_BIOTIC_HABITAT/", species, "_mod_bart.rds"))
  
  gc()
}  # end for spc
