# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(modEvA)
library(kuenm)
library(raster)
#library(ecospat)


# GET TARGET SPECIES ####

occ_files <- list.files("../../Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

#species <- species_list[1]


# CREATE INVASIVE EVALUATION DATASETS ####

dir.create("../inv_region_tables/")
dir.create("../MESS_JPGs")

for (species in species_list) {
  
  message("\n", species)
  
  # import and join invasive occurrence data:
  occ_inv <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_invasive.csv"))
  acc_inv <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_invaded_region.csv"))
  occ_inv$presence <- 1
  acc_inv$presence <- 0
  dat_inv <- rbind(occ_inv, acc_inv)
  
  # add model predictions:
  dat_inv <- vect(dat_inv, keepgeom = TRUE)
  glm_pred_inv_rst <- rast(paste0("../pred_rasters/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rasters/invasive/", species, "_inv_p_bart.tif"))
  dat_inv$glm_pred_inv <- terra::extract(glm_pred_inv_rst, dat_inv)[,2]
  dat_inv$bart_pred_inv <- terra::extract(bart_pred_inv_rst, dat_inv)[,2]
  # head(dat_inv)
  
  # add invaded region info:
  range_inv <- vect(paste0("../../ranges/", species, "/for_", species, "_europe_inv.shp"))
  # plot(dat_inv, "presence", cex = 0.3, col = c("grey", "black"))
  # plot(range_inv, add = TRUE)
  dat_inv_invaded_range <- dat_inv[range_inv, ]
  dat_inv$inv_range <- 0
  dat_inv$inv_range[dat_inv$X %in% dat_inv_invaded_range$X] <- 1
  # plot(dat_inv, "inv_range")
  # head(dat_inv)
  
  # add MESS info:
  dat_inv <- as.data.frame(dat_inv)
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  dat_nat <- rbind(occ_nat, acc_nat)
  # head(dat_nat)
  # dat_inv$ExDet <- ecospat.climan(dat_nat[ , c("PC1", "PC2")], dat_inv[ , c("PC1", "PC2")])
  # dat_inv$analog <- ifelse(dat_inv$ExDet >= 0 & dat_inv$ExDet <= 1, 1, 0)
  jpeg(paste0("../MESS_JPGs/", species, "_MESS.jpg"))
  dat_inv$mess <- MESS(dat_nat[ , c("PC1", "PC2")], dat_inv[ , c("PC1", "PC2")])$TOTAL
  dat_inv$analog <- ifelse(dat_inv$mess >= 0, 1, 0)
  # head(dat_inv)
  plot(vect(dat_inv, geom = c("lon", "lat"), crs = "epsg:4326"), "mess", type = "continuous", main = species)
  points(subset(dat_inv, analog == 1, select = c("lon", "lat")), pch = 19, cex = 0.01, col = "black")
  points(subset(dat_inv, presence == 1, select = c("lon", "lat")), pch = 20, cex = 0.8, col = "red")
  plot(range_inv, border = "red", add = TRUE)
  legend("topleft", legend = c("presence", "analog climate"), col = c("red", "black"), pch = c(20, 20), cex = 0.8, bty = "n")
  dev.off()
  
  # export table to disk:
  write.csv(dat_inv, paste0("../inv_region_tables/", species, "_inv_region_table.csv"), row.names = FALSE)
  
}; gc()


# EVALUATE NATIVE MODELS ON INVADED REGIONS ####

# try out MOP
# source("MOP_Owens_et_al_Supp_Mat.R")  # I made some edits to fix errors that arose initially
# mop_list <- MOP_NB(m1 = dat_nat[,-1], m2 = dat_inv[-1], c1 = c("PC1", "PC2"), c2 = c("PC1", "PC2"), decil = as.character(10), p1 = 0.5, p2 = 0.5, Xcol = "lon", Ycol = "lat", MxMESS = "Y")
# 
# lapply(mop_list, head)
# but it requires subjective case-specific choices for 'decil', 'p1', 'p2'


# create empty tables to receive the evaluation results:

eval_metrics_eur_full <- eval_metrics_eur_mess <- eval_metrics_inv_full <- eval_metrics_inv_mess <- as.data.frame(matrix(NA, nrow = length(species_list), ncol = 5 * 2))
rownames(eval_metrics_eur_full) <- rownames(eval_metrics_eur_mess) <- rownames(eval_metrics_inv_full) <- rownames(eval_metrics_inv_mess) <- species_list
colnames(eval_metrics_eur_full) <- colnames(eval_metrics_eur_mess) <- colnames(eval_metrics_inv_full) <- colnames(eval_metrics_inv_mess) <- c("AUCratio_glm", "AUCratio_bart", "pROCpval_glm", "pROCpval_bart", "omis_glm", "omis_bart", "sens_glm", "sens_bart", "spec_glm", "spec_bart")
# eval_metrics_eur_full <- as.data.frame(eval_metrics_eur_full)
# eval_metrics_eur_mess <- as.data.frame(eval_metrics_mess)
head(eval_metrics_eur_full)
head(eval_metrics_inv_mess)


# fill the tables with the results:

for (species in species_list) {
  
  message("\n", species)
  
  # import native occ+pred and get 2.5% minimum training presence threshold:
  message("computing 2.5% training omission threshold...")
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  glm_pred_nat_rst <- rast(paste0("../pred_rasters/native/", species, "_nat_p_glm.tif"))
  bart_pred_nat_rst <- rast(paste0("../pred_rasters/native/", species, "_nat_p_bart.tif"))
  glm_pred_occnat <- terra::extract(glm_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  bart_pred_occnat <- terra::extract(bart_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  glm_min_2.5_thresh <- quantile(glm_pred_occnat, probs = 0.025, na.rm = TRUE)
  bart_min_2.5_thresh <- quantile(bart_pred_occnat, probs = 0.025, na.rm = TRUE)
  
  # import invasive occurrence data and pred rasters:
  dat_inv <- read.csv(paste0("../inv_region_tables/", species, "_inv_region_table.csv"))
  glm_pred_inv_rst <- rast(paste0("../pred_rasters/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rasters/invasive/", species, "_inv_p_bart.tif"))
  
  message("evaluating on full European background...")
  bg <- dat_inv
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = glm_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = bart_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  glm_threshmeasures
  bart_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  
  message("evaluating on invaded European background...")
  bg <- dat_inv[dat_inv$inv_range == 1, ]
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = glm_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = bart_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  # glm_threshmeasures
  # bart_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures

  message("evaluating on MESS European background...")
  bg <- dat_inv[dat_inv$analog == 1, ]
  eval_metrics_eur_mess[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = glm_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_mess[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = bart_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  glm_threshmeasures
  bart_threshmeasures
  eval_metrics_eur_mess[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_mess[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  
  message("evaluating on MESS invaded background...")
  bg <- dat_inv[dat_inv$inv_range == 1 & dat_inv$analog == 1, ]
  eval_metrics_inv_mess[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = glm_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_mess[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = bart_min_2.5_thresh * 100, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  glm_threshmeasures
  bart_threshmeasures
  eval_metrics_inv_mess[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_mess[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
}; gc()


# see and plot eval metrics:

head(eval_metrics_eur_full)
head(eval_metrics_eur_mess)
head(eval_metrics_inv_full)
head(eval_metrics_inv_mess)

for (species in species_list) {
  barplot(as.matrix(eval_metrics_eur_full[species, ]), las = 2, main = paste(species, "Eur full"))
  barplot(as.matrix(eval_metrics_eur_mess[species, ]), las = 2, main = paste(species, "Eur MESS"))
  barplot(as.matrix(eval_metrics_inv_full[species, ]), las = 2, main = paste(species, "inv full"))
  barplot(as.matrix(eval_metrics_inv_mess[species, ]), las = 2, main = paste(species, "inv MESS"))
}


# export eval metrics:

dir.create("../eval_metrics")
write.csv(data.frame(species = rownames(eval_metrics_eur_full), eval_metrics_eur_full), "../eval_metrics/eval_metrics_eur_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_eur_mess), eval_metrics_eur_mess), "../eval_metrics/eval_metrics_eur_mess.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_full), eval_metrics_inv_full), "../eval_metrics/eval_metrics_inv_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_mess), eval_metrics_inv_mess), "../eval_metrics/eval_metrics_inv_mess.csv", row.names = FALSE)

