# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(fuzzySim)
library(modEvA)
library(maps)


# DATA ####

occ_files <- list.files("../../Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

#species_list <- c("syrrev", "vidmac", "psikra")
species_list <- species_list[c(1:17, 19:20, 18)]

dir.create("../pred_rasters/native", recursive = TRUE)
dir.create("../pred_rasters/invaded", recursive = TRUE)
dir.create("../pred_JPGs", recursive = TRUE)
dir.create("../model_objects/", recursive = TRUE)


# compute PA models ####

mods_glm <- vector("list", length(species_list))
mods_bart <- vector("list", length(species_list))
names(mods_glm) <- names(mods_bart) <- species_list

for (species in species_list) {
  message(species)
  
  # import occurrences and accessible (background) points:
  # from native areas:
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  head(occ_nat)
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  head(acc_nat)
  
  # import climgrids:
  climgrids_nat <- rast(paste0("../../occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  #plot(climgrids_nat)
  
  # get presence (occupied) and absence (accessible not occupied) pixels:
  occ_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(occ_nat[ , c("lon", "lat")]))
  acc_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(acc_nat[ , c("lon", "lat")]))
  # length(occ_nat_cells)
  # length(acc_nat_cells)
  abs_nat <- acc_nat[!(acc_nat_cells %in% occ_nat_cells), ]
  occ_nat$presence <- 1
  abs_nat$presence <- 0
  pa_nat <- rbind(occ_nat, abs_nat)
  # head(pa_nat)
  # sum(pa_nat$presence)
  # sum(pa_nat$presence == 0)
  
  #names(pa_nat)
  #var_cols
  
  # compute presence-absence models:
  var_cols <- names(pa_nat)[grep("PC", names(pa_nat))]

  form_glm <- as.formula(paste("presence ~", paste(var_cols, collapse = "+")))
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  #summary(mod_glm)
  saveRDS(mod_glm, paste0("../model_objects/", species, "_mod_glm"))
  
  mod_bart <- bart(x.train = pa_nat[ , var_cols], y.train = pa_nat[ , "presence"], keeptrees = TRUE)
  #summary(mod_bart)
  invisible(mod_bart$fit$state)
  saveRDS(mod_bart, paste0("../model_objects/", species, "_mod_bart"))
  
  
  mods_glm[[species]] <- mod_glm
  mods_bart[[species]] <- mod_bart
  
  # get presence-absence model predictions:
  
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
  message("predicting...")
  message(" -GLM")
  names(climgrids_nat) <- var_cols
  glm_p <- predict(climgrids_nat, mod_glm, type = "response")
  message(" -BART")
  bart_p <- predict(mod_bart, stack(climgrids_nat))
  bart_p <- rast(bart_p)
  
  # convert to favourability:
  
  message("computing favourability...")
  glm_f <- Fav(pred = glm_p, sample.preval = prevalence(model = mod_glm))
  bart_f <- Fav(pred = bart_p, sample.preval = prevalence(model = mod_bart))
  
  message("saving JPEGs...")
  jpeg(paste0("../pred_JPGs/", species, "_glm_f.jpg"), width = 600, height = 480)
  maps::map("world", col = "grey")
  plot(glm_f, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(paste(species, "GLM fav"))
  dev.off()
  
  jpeg(paste0("../pred_JPGs/", species, "_bart_f.jpg"), width = 600, height = 480)
  maps::map("world", col = "grey")
  plot(bart_f, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(paste(species, "BART fav"))
  dev.off()

  # plot(glm_f, main = paste(species, "GLM fav"))
  # plot(bart_f, main = paste(species, "BART fav"))
  
  # plot(pa_nat_sv, "presence", cex = 0.1, col = c(NA, "black"))
  # plot(glm_p, col = hcl.colors(100))
  # plot(bart_p, col = hcl.colors(100))
  
  message("saving rasters...")
  writeRaster(glm_p, paste0("../pred_rasters/native/", species, "_glm_p.tif"))
  writeRaster(bart_p, paste0("../pred_rasters/native/", species, "_bart_p.tif"))
  writeRaster(glm_f, paste0("../pred_rasters/native/", species, "_glm_f.tif"))
  writeRaster(bart_f, paste0("../pred_rasters/native/", species, "_bart_f.tif"))
  
  gc()
}  # end for spc
