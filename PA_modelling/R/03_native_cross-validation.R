# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(blockCV)
library(modEvA)


# GET TARGET SPECIES ####

occ_files <- list.files("../../Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# DIVIDE EACH NATIVE RANGE INTO SPATIAL BLOCKS ####

block_size <- 200000

blocks <- vector("list", length(species_list))
names(blocks) <- species_list

for (species in species_list) {
  message("\n", species, "...")

  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  #head(occ_nat)
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  #head(acc_nat)
  species_data <- vect(rbind(occ_nat[ , c("lon", "lat")], acc_nat[ , c("lon", "lat")]))
  
  blocks[[species]] <- spatialBlock(speciesData = as(species_data, "Spatial"), theRange = block_size, k = 5, selection = "systematic")
}


spc <- "agaper"
blocks
blocks[[spc]]
plot(blocks[[spc]]$blocks)
blocks[[spc]]$plots
blocks[[spc]]$folds
blocks[[spc]]$foldID


# GET PREDICTIONS FOR BLOCK CROSS-VALIDATION ####

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the BART predict function to work with data frames instead of just raster layers

dir.create("../pred_CSVs")

for (species in species_list) {
  message("\n", species)
  
  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  #head(occ_nat)
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  #head(acc_nat)
  
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
  # sum(pa_nat$presence == 1)
  # sum(pa_nat$presence == 0)
  
  # add spatial block ID:
  pa_nat$foldID <- blocks[[species]]$foldID
  
  # compute cross-validation models:
  #names(pa_nat)
  #var_names <- names(pa_nat)[grep("^PC", names(pa_nat))]
  var_names <- c("PC1", "PC2")
  form_glm <- as.formula(paste("presence ~", paste(var_names, collapse = "+")))
  folds <- sort(unique(pa_nat$foldID))
  
  for (f in folds) {
    message("fold ", f)
    pa_train <- pa_nat[pa_nat$foldID != f, ]
    
    mod_glm <- glm(formula = form_glm, family = binomial, data = pa_train)
    mod_bart <- bart(x.train = pa_train[ , var_names], y.train = pa_train[ , "presence"], keeptrees = TRUE, verbose = FALSE)

    pa_nat[ , paste0("glm_p_fold", f)] <- predict(mod_glm, pa_nat, type = "response")
    pa_nat[ , paste0("bart_p_fold", f)] <- predict_bart_df(mod_bart, pa_nat)
    #head(pa_nat)
    
    gc()
  }  # end for f
  
  # get prediction for entire native dataset (no folds):
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  mod_bart <- bart(x.train = pa_nat[ , var_names], y.train = pa_nat[ , "presence"], keeptrees = TRUE, verbose = FALSE)
  pa_nat[ , paste0("glm_p")] <- predict(mod_glm, pa_nat, type = "response")
  pa_nat[ , paste0("bart_p")] <- predict_bart_df(mod_bart, pa_nat)
  
  # save result to disk:
  write.csv(pa_nat, paste0("../pred_CSVs/", species, "_pred_crossval.csv"), row.names = FALSE)

}  # end for species


# COMPUTE CROSS-VALIDATION METRICS ####

fold_cols <- grep("_fold", names(dat))
names(dat)[fold_cols]
measures <- c("AUC", "TSS", "MCS")

# create an empty table to receive the cross-validation results:
crossval <- as.data.frame(matrix(nrow = length(folds), ncol = length(measures) * length(names(models))))
colnames(crossval) <- c(outer(names(models), measures, FUN = paste, sep = "_"))
crossval  # for now it's only filled with NAs

par(mfrow = c(5, 2), mar = c(2, 2, 1.1, 1))

for (species in species_list) {
  pred_nat <- read.csv(paste0("../pred_CSVs/", species, "_pred_crossval.csv"))
  

}


for (m in names(models))  for (f in folds) {
  fold_name <- paste0("fold", f)
  fold_col <- names(dat)[grep(paste0(m, "_fold", f), names(dat))]
  fold_dat <- subset(dat, foldID == f)
  crossval[f, paste(m, "AUC", sep = "_")] <- AUC(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = paste(m, "AUC"))
  crossval[f, paste(m, "TSS", sep = "_")] <- threshMeasures(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], thresh = "preval", measures = "TSS", simplif = TRUE, standardize = FALSE, main = paste(m, "TSS"))
  crossval[f, paste(m, "MCS", sep = "_")] <- MillerCalib(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], main = paste(m, "Miller line"))$slope
}
# press the back arrow on the top left of your plotting window to see the fold evaluations for the different models

crossval

# plot the mean cross-validation performance for each model, but note the mean is quite limited information!
crossval_means <- sapply(crossval, mean, na.rm = TRUE)
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1))
barplot(crossval_means, ylim = c(0, max(crossval, na.rm = TRUE)), col = rep(1:length(measures), each = length(models)), las = 2)
abline(h = 1, col = "darkgrey", lty = 2)

# get more information with boxplots of the cross-validation metrics:
boxplot(crossval, col = rep(1:length(measures), each = length(names(models))), las = 2)
abline(h = 1, col = "darkgrey", lty = 2)  # remember Miller calibration slope (MCS) should ideally be close to 1 (not bigger = better)


