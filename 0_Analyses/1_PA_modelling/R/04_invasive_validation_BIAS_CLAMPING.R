# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(modEvA)
library(kuenm)
library(raster)
library(enmSdm)
library(foreach)
library(doParallel)

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)


# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

#species <- species_list[1]


# CREATE INVASIVE EVALUATION DATASETS ####

dir.create("../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/inv_region_tables/")
dir.create("../eval_metrics/INVASIVE/BIAS/MOP_JPGs")

foreach (species = species_list,.packages=c('terra','modEvA','kuenm','raster','enmSdm'))  %dopar% {
  
  message("\n", species)
  
  #import and join invasive occurrence data:
  occ_inv <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_invasive.csv"))
  acc_inv <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_invaded_region.csv"))
  occ_inv$presence <- 1
  acc_inv$presence <- 0
  dat_inv <- rbind(occ_inv, acc_inv)
  dat_inv_lon <- dat_inv$lon
  dat_inv_lat <- dat_inv$lat
  
  # add model predictions:
  dat_inv <- vect(dat_inv)
  glm_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_bart.tif"))
  fne_pred_inv_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_inv.tif"))
  dat_inv$glm_pred_inv <- terra::extract(glm_pred_inv_rst, dat_inv)[,2]
  dat_inv$bart_pred_inv <- terra::extract(bart_pred_inv_rst, dat_inv)[,2]
  dat_inv$fne_pred_inv <- terra::extract(fne_pred_inv_rst, dat_inv)[,2]
  # head(dat_inv)
  
  # add invaded region info:
  range_inv <- vect(paste0("../../../1_Data/ranges/", species, "/for_", species, "_europe_inv.shp"))
  # plot(dat_inv, "presence", cex = 0.3, col = c("grey", "black"))
  # plot(range_inv, add = TRUE)
  dat_inv_invaded_range <- dat_inv[range_inv, ]
  dat_inv$inv_range <- 0
  dat_inv$inv_range[dat_inv$X %in% dat_inv_invaded_range$X] <- 1
  # plot(dat_inv, "inv_range")
  # head(dat_inv)
  
  # add MOP info:
  dat_inv <- as.data.frame(dat_inv)
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  dat_nat <- rbind(occ_nat, acc_nat)
  #head(dat_nat)
  set1 <- dat_nat[ , c("PC1", "PC2")]
  set2 <- dat_inv[ , c("PC1", "PC2")]
  mop_area <- mop(set1,set2, p=c(0.10), index=TRUE)
  jpeg(paste0("../eval_metrics/INVASIVE/BIAS/MOP_JPGs/", species, "_MOP.jpg"))
  
  #add similar datapoints to df
  native <- data.frame(unlist(mop_area$set1))
  native$similar <- rep("yes",nrow(native))
  colnames(native)[1] <- "ID"
  invasive <- data.frame(unlist(mop_area$set2))
  invasive$similar <- rep("yes",nrow(invasive))
  colnames(invasive)[1] <- "ID"
  
  dat_nat$ID <- seq(1:nrow(dat_nat))
  dat_inv$ID <- seq(1:nrow(dat_inv))
  
  dat_inv <- merge(dat_inv,invasive,by='ID',all=TRUE)
  dat_inv[is.na(dat_inv)] <- "no"
  dat_inv$lon <- dat_inv_lon
  dat_inv$lat <- dat_inv_lat
  
  dat_inv$mop <- ifelse(dat_inv$similar == "yes",1,0)
  dat_inv$analog <- ifelse(dat_inv$similar == "no",1,0)

  # head(dat_inv)
  plot(vect(dat_inv, geom = c("lon", "lat"), crs = "epsg:4326"), "mop", type = "continuous", main = species)
  points(subset(dat_inv, analog == 1, select = c("lon", "lat")), pch = 19, cex = 0.01, col = "black")
  points(subset(dat_inv, presence == 1, select = c("lon", "lat")), pch = 20, cex = 0.8, col = "red")
  plot(range_inv, border = "red", add = TRUE)
  legend("topleft", legend = c("presence", "analog climate"), col = c("red", "black"), pch = c(20, 20), cex = 0.8, bty = "n")
  dev.off()
  
  # export table to disk:
  write.csv(dat_inv, paste0("../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/inv_region_tables/", species, "_inv_region_table.csv"), row.names = FALSE)
  
}; gc()

remove(cl)

# EVALUATE NATIVE MODELS ON INVADED REGIONS ####

# create empty tables to receive the evaluation results:

eval_metrics_eur_full <- eval_metrics_eur_mop <- eval_metrics_inv_full <- eval_metrics_inv_mop <- as.data.frame(matrix(NA, nrow = length(species_list), ncol = 5 * 3))
rownames(eval_metrics_eur_full) <- rownames(eval_metrics_eur_mop) <- rownames(eval_metrics_inv_full) <- rownames(eval_metrics_inv_mop) <- species_list
colnames(eval_metrics_eur_full) <- colnames(eval_metrics_eur_mop) <- colnames(eval_metrics_inv_full) <- colnames(eval_metrics_inv_mop) <- c("AUCratio_glm", "AUCratio_bart","AUCratio_fne", "pROCpval_glm", "pROCpval_bart","pROCpval_fne","omis_glm", "omis_bart","omis_fne","sens_glm", "sens_bart", "sens_fne","spec_glm", "spec_bart","spec_fne")
# eval_metrics_eur_full <- as.data.frame(eval_metrics_eur_full)
# eval_metrics_eur_mop <- as.data.frame(eval_metrics_mop)
head(eval_metrics_eur_full)
head(eval_metrics_inv_mop)


mtp_thresholds <- data.frame(species = species_list, glm_thresh = NA_real_, bart_thresh = NA_real_,fne_thresh = NA_real_)


# fill the tables with the results:

for (species in species_list) {
  
  message("\n", species)
  
  # import native occ+pred and get 5% minimum training presence threshold:
  message("computing 5% training omission threshold...")
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  glm_pred_nat_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_glm.tif"))
  bart_pred_nat_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_bart.tif"))
  fne_pred_nat_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_nat.tif"))
  glm_pred_occnat <- terra::extract(glm_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  bart_pred_occnat <- terra::extract(bart_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  fne_pred_occnat <- terra::extract(fne_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  glm_min_5_thresh <- quantile(glm_pred_occnat, probs = 0.05, na.rm = TRUE)
  bart_min_5_thresh <- quantile(bart_pred_occnat, probs = 0.05, na.rm = TRUE)
  fne_min_5_thresh <- quantile(fne_pred_occnat, probs = 0.05, na.rm = TRUE)
  mtp_thresholds[mtp_thresholds$species == species, "glm_thresh"] <- glm_min_5_thresh
  mtp_thresholds[mtp_thresholds$species == species, "bart_thresh"] <- bart_min_5_thresh
  mtp_thresholds[mtp_thresholds$species == species, "fne_thresh"] <- fne_min_5_thresh
  
  # import invasive occurrence data and pred rasters:
  dat_inv <- read.csv(paste0("../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/inv_region_tables/", species, "_inv_region_table.csv"))
  glm_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_bart.tif"))
  fne_pred_inv_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_inv.tif"))
  
  message("evaluating on full European background...")
  bg <- dat_inv
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
  message("evaluating on invaded European background...")
  bg <- dat_inv[dat_inv$inv_range == 1, ]
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  
   glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  # glm_threshmeasures
  # bart_threshmeasures
  # fne_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  

  message("evaluating on MOP European background...")
  bg <- dat_inv[dat_inv$analog == 1, ]
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
  message("evaluating on MOP invaded background...")
  bg <- dat_inv[dat_inv$inv_range == 1 & dat_inv$analog == 1, ]
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
}; gc()


# see and plot eval metrics:

mtp_thresholds
range(sort(mtp_thresholds$glm_thresh))
range(sort(mtp_thresholds$bart_thresh))

head(eval_metrics_eur_full)
head(eval_metrics_eur_mop)
head(eval_metrics_inv_full)
head(eval_metrics_inv_mop)

for (species in species_list) {
  barplot(as.matrix(eval_metrics_eur_full[species, ]), las = 2, main = paste(species, "Eur full"))
  barplot(as.matrix(eval_metrics_eur_mop[species, ]), las = 2, main = paste(species, "Eur MOP"))
  barplot(as.matrix(eval_metrics_inv_full[species, ]), las = 2, main = paste(species, "inv full"))
  barplot(as.matrix(eval_metrics_inv_mop[species, ]), las = 2, main = paste(species, "inv MOP"))
}


# export mtp_thresholds and eval metrics:

temp <- as.data.frame(matrix(NA, nrow = length(mtp_thresholds$species), ncol = 4))
temp$V1 <- rep("yes",nrow(temp))
temp$V2 <- rep("climate",nrow(temp))
temp$V3 <- rep("clamping",nrow(temp))
temp$V4 <- rep("5",nrow(temp))
colnames(temp) <- c("BIAS","preds","extrapol","threshold")

mtp_thresholds <- cbind(mtp_thresholds,temp)

eval_metrics_eur_mop <- cbind(eval_metrics_eur_mop,temp)
  eval_metrics_eur_mop$mop <- rep("yes_mop",nrow(eval_metrics_eur_mop))

eval_metrics_eur_full <- cbind(eval_metrics_eur_full,temp)
  eval_metrics_eur_full$mop <- rep("no",nrow(eval_metrics_eur_full))

eval_metrics_inv_full <- cbind(eval_metrics_inv_full,temp)
  eval_metrics_inv_full$mop <- rep("no",nrow(eval_metrics_inv_full))

eval_metrics_inv_mop <- cbind(eval_metrics_inv_mop,temp)
  eval_metrics_inv_mop$mop <- rep("yes_mop",nrow(eval_metrics_inv_mop))

write.csv(mtp_thresholds, "../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/mtp_thresholds.csv", row.names = FALSE)

write.csv(data.frame(species = rownames(eval_metrics_eur_full), eval_metrics_eur_full), "../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/eval_metrics_eur_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_eur_mop), eval_metrics_eur_mop), "../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/eval_metrics_eur_mop.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_full), eval_metrics_inv_full), "../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/eval_metrics_inv_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_mop), eval_metrics_inv_mop), "../eval_metrics/INVASIVE/BIAS/CLAMPING/5E/eval_metrics_inv_mop.csv", row.names = FALSE)

##################################################################################################################################################################################################


rm(list = ls())

# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# CREATE INVASIVE EVALUATION DATASETS ####

dir.create("../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/inv_region_tables/")
dir.create("../eval_metrics/INVASIVE/BIAS/MOP_JPGs")

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)


foreach (species = species_list,.packages=c('terra','modEvA','kuenm','raster','enmSdm'))  %dopar% {
  
  
  message("\n", species)
  
  # import and join invasive occurrence data:
  occ_inv <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_invasive.csv"))
  acc_inv <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_invaded_region.csv"))
  occ_inv$presence <- 1
  acc_inv$presence <- 0
  dat_inv <- rbind(occ_inv, acc_inv)
  dat_inv_lon <- dat_inv$lon
  dat_inv_lat <- dat_inv$lat
  
  # add model predictions:
  dat_inv <- vect(dat_inv)
  glm_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_bart.tif"))
  fne_pred_inv_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_inv.tif"))
  dat_inv$glm_pred_inv <- terra::extract(glm_pred_inv_rst, dat_inv)[,2]
  dat_inv$bart_pred_inv <- terra::extract(bart_pred_inv_rst, dat_inv)[,2]
  dat_inv$fne_pred_inv <- terra::extract(fne_pred_inv_rst, dat_inv)[,2]
  # head(dat_inv)
  
  # add invaded region info:
  range_inv <- vect(paste0("../../../1_Data/ranges/", species, "/for_", species, "_europe_inv.shp"))
  # plot(dat_inv, "presence", cex = 0.3, col = c("grey", "black"))
  # plot(range_inv, add = TRUE)
  dat_inv_invaded_range <- dat_inv[range_inv, ]
  dat_inv$inv_range <- 0
  dat_inv$inv_range[dat_inv$X %in% dat_inv_invaded_range$X] <- 1
  # plot(dat_inv, "inv_range")
  # head(dat_inv)
  
  # add MOP info:
  dat_inv <- as.data.frame(dat_inv)
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  dat_nat <- rbind(occ_nat, acc_nat)
  #head(dat_nat)
  set1 <- dat_nat[ , c("PC1", "PC2")]
  set2 <- dat_inv[ , c("PC1", "PC2")]
  mop_area <- mop(set1,set2, p=c(0.10), index=TRUE)
  jpeg(paste0("../eval_metrics/INVASIVE/BIAS/MOP_JPGs/", species, "_MOP.jpg"))
  
  #add similar datapoints to df
  native <- data.frame(unlist(mop_area$set1))
  native$similar <- rep("yes",nrow(native))
  colnames(native)[1] <- "ID"
  invasive <- data.frame(unlist(mop_area$set2))
  invasive$similar <- rep("yes",nrow(invasive))
  colnames(invasive)[1] <- "ID"
  
  dat_nat$ID <- seq(1:nrow(dat_nat))
  dat_inv$ID <- seq(1:nrow(dat_inv))
  
  dat_inv <- merge(dat_inv,invasive,by='ID',all=TRUE)
  dat_inv[is.na(dat_inv)] <- "no"
  dat_inv$lon <- dat_inv_lon
  dat_inv$lat <- dat_inv_lat
  
  dat_inv$mop <- ifelse(dat_inv$similar == "yes",1,0)
  dat_inv$analog <- ifelse(dat_inv$similar == "no",1,0)
  
  # head(dat_inv)
  plot(vect(dat_inv, geom = c("lon", "lat"), crs = "epsg:4326"), "mop", type = "continuous", main = species)
  points(subset(dat_inv, analog == 1, select = c("lon", "lat")), pch = 19, cex = 0.01, col = "black")
  points(subset(dat_inv, presence == 1, select = c("lon", "lat")), pch = 20, cex = 0.8, col = "red")
  plot(range_inv, border = "red", add = TRUE)
  legend("topleft", legend = c("presence", "analog climate"), col = c("red", "black"), pch = c(20, 20), cex = 0.8, bty = "n")
  dev.off()
  
  # export table to disk:
  write.csv(dat_inv, paste0("../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/inv_region_tables/", species, "_inv_region_table.csv"), row.names = FALSE)
  
}; gc()


# EVALUATE NATIVE MODELS ON INVADED REGIONS ####

# create empty tables to receive the evaluation results:

eval_metrics_eur_full <- eval_metrics_eur_mop <- eval_metrics_inv_full <- eval_metrics_inv_mop <- as.data.frame(matrix(NA, nrow = length(species_list), ncol = 5 * 3))
rownames(eval_metrics_eur_full) <- rownames(eval_metrics_eur_mop) <- rownames(eval_metrics_inv_full) <- rownames(eval_metrics_inv_mop) <- species_list
colnames(eval_metrics_eur_full) <- colnames(eval_metrics_eur_mop) <- colnames(eval_metrics_inv_full) <- colnames(eval_metrics_inv_mop) <- c("AUCratio_glm", "AUCratio_bart","AUCratio_fne", "pROCpval_glm", "pROCpval_bart","pROCpval_fne","omis_glm", "omis_bart","omis_fne","sens_glm", "sens_bart", "sens_fne","spec_glm", "spec_bart","spec_fne")
# eval_metrics_eur_full <- as.data.frame(eval_metrics_eur_full)
# eval_metrics_eur_mop <- as.data.frame(eval_metrics_mop)
head(eval_metrics_eur_full)
head(eval_metrics_inv_mop)


mtp_thresholds <- data.frame(species = species_list, glm_thresh = NA_real_, bart_thresh = NA_real_,fne_thresh = NA_real_)


# fill the tables with the results:

for (species in species_list) {
  
  message("\n", species)
  
  # import native occ+pred and get 2.5% minimum training presence threshold:
  message("computing 2.5% training omission threshold...")
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  glm_pred_nat_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_glm.tif"))
  bart_pred_nat_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/native/", species, "_nat_p_bart.tif"))
  fne_pred_nat_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_nat.tif"))
  glm_pred_occnat <- terra::extract(glm_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  bart_pred_occnat <- terra::extract(bart_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  fne_pred_occnat <- terra::extract(fne_pred_nat_rst, occ_nat[ , c("lon", "lat")])[,2]
  glm_min_2.5_thresh <- quantile(glm_pred_occnat, probs = 0.025, na.rm = TRUE)
  bart_min_2.5_thresh <- quantile(bart_pred_occnat, probs = 0.025, na.rm = TRUE)
  fne_min_2.5_thresh <- quantile(fne_pred_occnat, probs = 0.025, na.rm = TRUE)
  mtp_thresholds[mtp_thresholds$species == species, "glm_thresh"] <- glm_min_2.5_thresh
  mtp_thresholds[mtp_thresholds$species == species, "bart_thresh"] <- bart_min_2.5_thresh
  mtp_thresholds[mtp_thresholds$species == species, "fne_thresh"] <- fne_min_2.5_thresh
  
  # import invasive occurrence data and pred rasters:
  dat_inv <- read.csv(paste0("../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/inv_region_tables/", species, "_inv_region_table.csv"))
  glm_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_glm.tif"))
  bart_pred_inv_rst <- rast(paste0("../pred_rastersBIAS_CLAMPING/EU/invasive/", species, "_inv_p_bart.tif"))
  fne_pred_inv_rst <- rast(paste0("../../0_Nf_modeling/ResultsClimateClamping/", species, "_suitability_map_inv.tif"))
  
  message("evaluating on full European background...")
  bg <- dat_inv
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_full[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_eur_full[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
  message("evaluating on invaded European background...")
  bg <- dat_inv[dat_inv$inv_range == 1, ]
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_full[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  # glm_threshmeasures
  # bart_threshmeasures
  # fne_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_inv_full[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
  
  message("evaluating on MOP European background...")
  bg <- dat_inv[dat_inv$analog == 1, ]
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_eur_mop[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_eur_mop[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
  message("evaluating on MOP invaded background...")
  bg <- dat_inv[dat_inv$inv_range == 1 & dat_inv$analog == 1, ]
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_glm")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(glm_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_bart")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(bart_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  eval_metrics_inv_mop[species, paste0(c("AUCratio", "pROCpval"), "_fne")] <- kuenm_proc(occ.test = bg[bg$presence == 1, c("lon", "lat")], model = raster(fne_pred_inv_rst), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary  # or threshold = 2.5?
  
  glm_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = glm_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = glm_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  bart_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = bart_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = bart_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  fne_threshmeasures <- with(bg, threshMeasures(obs = presence, pred = fne_pred_inv, measures = c("Omission", "Sensitivity", "Specificity"), thresh = fne_min_2.5_thresh, plot = FALSE))$ThreshMeasures[,1]
  
  glm_threshmeasures
  bart_threshmeasures
  fne_threshmeasures
  
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_glm")] <- glm_threshmeasures
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_bart")] <- bart_threshmeasures
  eval_metrics_inv_mop[species, paste0(c("omis", "sens", "spec"), "_fne")] <- fne_threshmeasures
  
}; gc()


# see and plot eval metrics:

mtp_thresholds
range(sort(mtp_thresholds$glm_thresh))
range(sort(mtp_thresholds$bart_thresh))

head(eval_metrics_eur_full)
head(eval_metrics_eur_mop)
head(eval_metrics_inv_full)
head(eval_metrics_inv_mop)

for (species in species_list) {
  barplot(as.matrix(eval_metrics_eur_full[species, ]), las = 2, main = paste(species, "Eur full"))
  barplot(as.matrix(eval_metrics_eur_mop[species, ]), las = 2, main = paste(species, "Eur MOP"))
  barplot(as.matrix(eval_metrics_inv_full[species, ]), las = 2, main = paste(species, "inv full"))
  barplot(as.matrix(eval_metrics_inv_mop[species, ]), las = 2, main = paste(species, "inv MOP"))
}


# export mtp_thresholds and eval metrics:

temp <- as.data.frame(matrix(NA, nrow = length(mtp_thresholds$species), ncol = 4))
temp$V1 <- rep("yes",nrow(temp))
temp$V2 <- rep("climate",nrow(temp))
temp$V3 <- rep("clamping",nrow(temp))
temp$V4 <- rep("2.5",nrow(temp))
colnames(temp) <- c("BIAS","preds","extrapol","threshold")

mtp_thresholds <- cbind(mtp_thresholds,temp)

eval_metrics_eur_mop <- cbind(eval_metrics_eur_mop,temp)
  eval_metrics_eur_mop$mop <- rep("yes_mop",nrow(eval_metrics_eur_mop))
  
eval_metrics_eur_full <- cbind(eval_metrics_eur_full,temp)
  eval_metrics_eur_full$mop <- rep("no",nrow(eval_metrics_eur_full))
  
eval_metrics_inv_full <- cbind(eval_metrics_inv_full,temp)
  eval_metrics_inv_full$mop <- rep("no",nrow(eval_metrics_inv_full))
  
eval_metrics_inv_mop <- cbind(eval_metrics_inv_mop,temp)
  eval_metrics_inv_mop$mop <- rep("yes_mop",nrow(eval_metrics_inv_mop))

write.csv(mtp_thresholds, "../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/mtp_thresholds.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_eur_full), eval_metrics_eur_full), "../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/eval_metrics_eur_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_eur_mop), eval_metrics_eur_mop), "../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/eval_metrics_eur_mop.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_full), eval_metrics_inv_full), "../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/eval_metrics_inv_full.csv", row.names = FALSE)
write.csv(data.frame(species = rownames(eval_metrics_inv_mop), eval_metrics_inv_mop), "../eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/eval_metrics_inv_mop.csv", row.names = FALSE)

remove(cl)



