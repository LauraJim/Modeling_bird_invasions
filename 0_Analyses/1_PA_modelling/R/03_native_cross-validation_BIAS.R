# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

# LOAD PACKAGES ####

library(terra)
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


# GET SPATIAL BLOCK PREDICTIONS FOR CROSS-VALIDATION ####

dir.create("../eval_metrics/NATIVE/BIAS/5E/CV_blocksBIAS")
dir.create("../eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS")

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the BART predict function to work with data frames instead of raster layers

block_size <- 200000

blocks <- vector("list", length(species_list))
names(blocks) <- species_list


# species <- "agaper"


nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

# species <- "agaper"

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm','blockCV'))  %dopar% {
  #species_list[(grep(species, species_list)+1):length(species_list)]
  message("\n", species)
  
  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  #head(occ_nat)
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  #head(acc_nat)
  
  #import biasgrid
  biasgrid_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_biasgrids/", species, "/", species, "_bias.tif"))
    occ_nat$bias <- terra::extract(biasgrid_nat, occ_nat[c(2:3)])[2]$layer
    acc_nat$bias <- terra::extract(biasgrid_nat, acc_nat[c(2:3)])[2]$layer

  head(occ_nat)
  head(acc_nat)
  
  # import climgrids:
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  #plot(climgrids_nat)
  
  #predictorvars
  mask <- raster(climgrids_nat[[1]]) + raster(climgrids_nat[[2]])  + raster(biasgrid_nat)
  mask <- mask*0
  PC1 <- mask+raster(climgrids_nat[[1]])
  PC2 <- mask+raster(climgrids_nat[[2]])
  bias <- (mask+raster(biasgrid_nat))*0

  
  predictorvars <- stack(PC1,PC2,bias)
  plot(predictorvars)
  names(predictorvars) <- var_names <- c("PC1", "PC2","bias")
  
  # get presence (occupied) and absence (accessible not occupied) pixels:
  occ_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(occ_nat[ , c("lon", "lat")]))
  acc_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(acc_nat[ , c("lon", "lat")]))
  # length(occ_nat_cells)
  # length(acc_nat_cells)
  abs_nat <- acc_nat[!(acc_nat_cells %in% occ_nat_cells), ]
  occ_nat$presence <- 1
  abs_nat$presence <- 0
  pa_nat <- na.omit(rbind(occ_nat, abs_nat))
  # head(pa_nat)
  # sum(pa_nat$presence == 1)
  # sum(pa_nat$presence == 0)
  pa_nat_sv <- vect(pa_nat, geom = c("lon", "lat"), crs = "epsg:4326")
  # nrow(pa_nat_sv)
  
  #block_size <- expanse(convHull(pa_nat_sv)) / 10000000

  message("computing spatial blocks...")
  blocks[[species]] <- spatialBlock(speciesData = as(pa_nat_sv, "Spatial"), species = "presence", theRange = block_size, k = 5, selection = "random", showBlocks = FALSE, seed = grep(species, species_list))
  
  # plot spatial blocks with occ data:
  jpeg(paste0("../eval_metrics/NATIVE/BIAS/5E/CV_blocksBIAS/", species, ".jpg"), width = 500, height = 500)
  plot(trim(climgrids_nat[[1]]), xlim = ext(blocks[[species]]$blocks)[1:2], ylim = ext(blocks[[species]]$blocks)[3:4], col = "wheat2", legend = FALSE, main = species, cex.main = 2)
  plot(blocks[[species]]$blocks, border = "grey40", add = TRUE)
  points(pa_nat[pa_nat$presence == 0, c("lon", "lat")], pch = 20, cex = 0.1, col = "salmon")
  text(blocks[[species]]$blocks, blocks[[species]]$blocks$folds, col = blocks[[species]]$blocks$folds, halo = TRUE)  # col = c("black", "blue", "red", "darkgreen", "purple")
  points(pa_nat[pa_nat$presence == 1, c("lon", "lat")], pch = 20, cex = 0.8, col = "blue")
  legend("topleft", c("presence", "absence"), pch = c(19, 20), col = c("blue", "salmon"), bty = "n")
  dev.off()
  
  # add spatial block ID:
  pa_nat$foldID <- blocks[[species]]$foldID
  
  # compute cross-validation models:
  message("computing fold models...")
  #names(pa_nat)
  #var_names <- names(pa_nat)[grep("^PC", names(pa_nat))]
  var_names <- c("PC1", "PC2","bias")
  form_glm <- as.formula(paste("presence ~", paste(var_names, collapse = "+")))
  folds <- sort(unique(pa_nat$foldID))
  
  for (f in folds) {
    message("fold ", f)
    pa_train <- pa_nat[pa_nat$foldID != f, ]
    
    mod_glm <- glm(formula = form_glm, family = binomial, data = pa_train)
    mod_bart <- bart(x.train = pa_train[ , var_names], y.train = pa_train[ , "presence"], keeptrees = TRUE, verbose = FALSE)
    
    pa_nat0bias <- pa_nat
    pa_nat0bias$bias <- pa_nat0bias$bias*0
    pa_nat[ , paste0("glm_p_fold", f)] <- predict(mod_glm, pa_nat0bias, type = "response")
    tempfile <- as.numeric(raster::extract(predict2.bart(mod_bart, predictorvars),pa_nat[c(2:3)]))
    pa_nat[ , paste0("bart_p_fold", f)] <- tempfile
    #head(pa_nat)
    
    gc()
  }  # end for f
  
  # get prediction for entire native dataset (no folds):
  message("full dataset")
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  mod_bart <- bart(x.train = pa_nat[ , var_names], y.train = pa_nat[ , "presence"], keeptrees = TRUE, verbose = FALSE)
  pa_nat[ , paste0("glm_p")] <- predict(mod_glm, pa_nat0bias, type = "response")
  pa_nat[ , paste0("bart_p")] <- predict_bart_df(mod_bart, pa_nat0bias)
  
  # save result to disk:
  write.csv(pa_nat, paste0("../eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"), row.names = FALSE)
  
}  # end for species

# spc <- "agaper"
# blocks
# blocks[[spc]]
# plot(blocks[[spc]]$blocks)
# blocks[[spc]]$plots
# blocks[[spc]]$folds
# blocks[[spc]]$foldID
# sapply(blocks[[spc]]$folds, function(x) length(x[[1]]))
# length(blocks[[spc]]$foldID)

#blocks
#saveRDS(blocks, "../eval_metrics/NATIVE/BIAS/5E/CV_blocksBIAS/_blocks.rds")

remove(cl)

# COMPUTE CROSS-VALIDATION METRICS ####

dir.create("../eval_metrics/NATIVE/BIAS/5E/crossval_boxplots", recursive = TRUE)
dir.create("../eval_metrics/NATIVE/BIAS/5E/crossval_metrics", recursive = TRUE)

models <- c("glm", "bart")
folds <- 1:5

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  pred_nat <- read.csv(paste0("../eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"))

  # create an empty table to receive this species' cross-validation results:
  crossval <- as.data.frame(matrix(nrow = length(folds), ncol = 3 + length(measures) * length(models)))
  colnames(crossval) <- c("n_pres", "n_abs", "thresh_5trainpres", outer(measures, models, FUN = paste, sep = "_"))

  for (m in models)  for (f in folds) {
    fold_p_col <- paste0(m, "_p_fold", f)
    fold_dat <- pred_nat[pred_nat$foldID == f, ]
    
    thresh <- quantile(pred_nat[pred_nat$presence == 1 & pred_nat$foldID != f, fold_p_col], probs = 0.05, na.rm = TRUE)  # pred for minimum 5% presence in the training data (outside fold)
    glm_pred_nat_rst <- rast(paste0("../pred_rastersBIAS/native/", species, "_nat_p_glm.tif"))
    bart_pred_nat_rst <- rast(paste0("../pred_rastersBIAS/native/", species, "_nat_p_bart.tif"))
    
    proc <- kuenm_proc(occ.test = fold_dat[fold_dat$presence == 1, c("lon", "lat")], model = raster(get(paste0(m, "_pred_nat_rst"))), threshold = 5, rand.percent = 50, iterations = 1000)$pROC_summary
    threshmeas <- threshMeasures(obs = fold_dat$presence, pred = fold_dat[ , fold_p_col], measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
    
    crossval[f, "n_pres"] <- sum(pred_nat[pred_nat$foldID == f, "presence"] == 1, na.rm = TRUE)
    crossval[f, "n_abs"] <- sum(pred_nat[pred_nat$foldID == f, "presence"] == 0, na.rm = TRUE)
    crossval[f, paste0("thresh", m)] <- thresh
    crossval[f, paste0(c("AUCratio_", "pROCpval_"), m)] <- proc
    crossval[f, paste0(c("omis_", "sens_", "spec_"), m)] <- threshmeas
  }
  temp <- data.frame(fold = as.integer(rownames(crossval)), crossval)
  write.csv(temp, paste0("../eval_metrics/NATIVE/BIAS/5E/crossval_metrics/crossval_metrics_nat_", species, ".csv"), row.names = FALSE)
  
  # boxplots of the cross-validation metrics:
  jpeg(paste0("../eval_metrics/NATIVE/BIAS/5E/crossval_boxplots/", species, ".jpg"))
  par(mar = c(7.2, 3, 2, 1))
  boxplot(crossval[ , 3:ncol(crossval)], las = 2, main = species)
  abline(h = 1, col = "darkgrey", lty = 2)
  dev.off()
  
  print(species)
  print(crossval)
  gc()
}

remove(cl)

################################################################################################################################################################################

# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# GET SPATIAL BLOCK PREDICTIONS FOR CROSS-VALIDATION ####

dir.create("../eval_metrics/NATIVE/BIAS/2.5E/CV_blocksBIAS")
dir.create("../eval_metrics/NATIVE/BIAS/2.5E/pred_CSVsBIAS")

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the BART predict function to work with data frames instead of raster layers

block_size <- 200000

blocks <- vector("list", length(species_list))
names(blocks) <- species_list


# species <- "agaper"


nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

# species <- "agaper"

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm','blockCV'))  %dopar% {
  #species_list[(grep(species, species_list)+1):length(species_list)]
  message("\n", species)
  
  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/", species, "_native.csv"))
  #head(occ_nat)
  acc_nat <- read.csv(paste0("../../0_Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  #head(acc_nat)
  
  #import biasgrid
  biasgrid_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_biasgrids/", species, "/", species, "_bias.tif"))
  occ_nat$bias <- terra::extract(biasgrid_nat, occ_nat[c(2:3)])[2]$layer
  acc_nat$bias <- terra::extract(biasgrid_nat, acc_nat[c(2:3)])[2]$layer

  head(occ_nat)
  head(acc_nat)
  
  # import climgrids:
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  #plot(climgrids_nat)
  
  #predictorvars
  mask <- raster(climgrids_nat[[1]]) + raster(climgrids_nat[[2]]) + raster(biasgrid_nat)
  mask <- mask*0
  PC1 <- mask+raster(climgrids_nat[[1]])
  PC2 <- mask+raster(climgrids_nat[[2]])
  bias <- (mask+raster(biasgrid_nat))*0
  
  predictorvars <- stack(PC1,PC2,bias)
  plot(predictorvars)
  names(predictorvars) <- var_names <- c("PC1", "PC2","bias")
  
  # get presence (occupied) and absence (accessible not occupied) pixels:
  occ_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(occ_nat[ , c("lon", "lat")]))
  acc_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(acc_nat[ , c("lon", "lat")]))
  # length(occ_nat_cells)
  # length(acc_nat_cells)
  abs_nat <- acc_nat[!(acc_nat_cells %in% occ_nat_cells), ]
  occ_nat$presence <- 1
  abs_nat$presence <- 0
  pa_nat <- na.omit(rbind(occ_nat, abs_nat))
  # head(pa_nat)
  # sum(pa_nat$presence == 1)
  # sum(pa_nat$presence == 0)
  pa_nat_sv <- vect(pa_nat, geom = c("lon", "lat"), crs = "epsg:4326")
  # nrow(pa_nat_sv)
  
  #block_size <- expanse(convHull(pa_nat_sv)) / 10000000
  
  message("computing spatial blocks...")
  blocks[[species]] <- spatialBlock(speciesData = as(pa_nat_sv, "Spatial"), species = "presence", theRange = block_size, k = 5, selection = "random", showBlocks = FALSE, seed = grep(species, species_list))
  
  # plot spatial blocks with occ data:
  jpeg(paste0("../eval_metrics/NATIVE/BIAS/2.5E/CV_blocksBIAS/", species, ".jpg"), width = 500, height = 500)
  plot(trim(climgrids_nat[[1]]), xlim = ext(blocks[[species]]$blocks)[1:2], ylim = ext(blocks[[species]]$blocks)[3:4], col = "wheat2", legend = FALSE, main = species, cex.main = 2)
  plot(blocks[[species]]$blocks, border = "grey40", add = TRUE)
  points(pa_nat[pa_nat$presence == 0, c("lon", "lat")], pch = 20, cex = 0.1, col = "salmon")
  text(blocks[[species]]$blocks, blocks[[species]]$blocks$folds, col = blocks[[species]]$blocks$folds, halo = TRUE)  # col = c("black", "blue", "red", "darkgreen", "purple")
  points(pa_nat[pa_nat$presence == 1, c("lon", "lat")], pch = 20, cex = 0.8, col = "blue")
  legend("topleft", c("presence", "absence"), pch = c(19, 20), col = c("blue", "salmon"), bty = "n")
  dev.off()
  
  # add spatial block ID:
  pa_nat$foldID <- blocks[[species]]$foldID
  
  # compute cross-validation models:
  message("computing fold models...")
  #names(pa_nat)
  #var_names <- names(pa_nat)[grep("^PC", names(pa_nat))]
  var_names <- c("PC1", "PC2","bias")
  form_glm <- as.formula(paste("presence ~", paste(var_names, collapse = "+")))
  folds <- sort(unique(pa_nat$foldID))
  
  for (f in folds) {
    message("fold ", f)
    pa_train <- pa_nat[pa_nat$foldID != f, ]
    
    mod_glm <- glm(formula = form_glm, family = binomial, data = pa_train)
    mod_bart <- bart(x.train = pa_train[ , var_names], y.train = pa_train[ , "presence"], keeptrees = TRUE, verbose = FALSE)
    
    pa_nat0bias <- pa_nat
    pa_nat0bias$bias <- pa_nat0bias$bias*0
    pa_nat[ , paste0("glm_p_fold", f)] <- predict(mod_glm, pa_nat0bias, type = "response")
    tempfile <- as.numeric(raster::extract(predict2.bart(mod_bart, predictorvars),pa_nat[c(2:3)]))
    pa_nat[ , paste0("bart_p_fold", f)] <- tempfile
    #head(pa_nat)
    
    gc()
  }  # end for f
  
  # get prediction for entire native dataset (no folds):
  message("full dataset")
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  mod_bart <- bart(x.train = pa_nat[ , var_names], y.train = pa_nat[ , "presence"], keeptrees = TRUE, verbose = FALSE)
  pa_nat[ , paste0("glm_p")] <- predict(mod_glm, pa_nat0bias, type = "response")
  pa_nat[ , paste0("bart_p")] <- predict_bart_df(mod_bart, pa_nat0bias)
  
  # save result to disk:
  write.csv(pa_nat, paste0("../eval_metrics/NATIVE/BIAS/2.5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"), row.names = FALSE)
  
}  # end for species

# spc <- "agaper"
# blocks
# blocks[[spc]]
# plot(blocks[[spc]]$blocks)
# blocks[[spc]]$plots
# blocks[[spc]]$folds
# blocks[[spc]]$foldID
# sapply(blocks[[spc]]$folds, function(x) length(x[[1]]))
# length(blocks[[spc]]$foldID)

#blocks
#saveRDS(blocks, "../eval_metrics/NATIVE/BIAS/2.5E/CV_blocksBIAS/_blocks.rds")

remove(cl)

# COMPUTE CROSS-VALIDATION METRICS ####

dir.create("../eval_metrics/NATIVE/BIAS/2.5E/crossval_boxplots", recursive = TRUE)
dir.create("../eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics", recursive = TRUE)

models <- c("glm", "bart")
folds <- 1:5

nc = 20
cl = makeCluster(nc)
registerDoParallel(cl)

measures <- c("AUCratio", "pROCpval", "omis", "sens", "spec")

foreach (species = species_list,.packages=c('terra','embarcadero','fuzzySim','modEvA','kuenm'))  %dopar% {
  pred_nat <- read.csv(paste0("../eval_metrics/NATIVE/BIAS/2.5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"))
  
  # create an empty table to receive this species' cross-validation results:
  crossval <- as.data.frame(matrix(nrow = length(folds), ncol = 3 + length(measures) * length(models)))
  colnames(crossval) <- c("n_pres", "n_abs", "thresh_2.5trainpres", outer(measures, models, FUN = paste, sep = "_"))
  
  for (m in models)  for (f in folds) {
    fold_p_col <- paste0(m, "_p_fold", f)
    fold_dat <- pred_nat[pred_nat$foldID == f, ]
    
    thresh <- quantile(pred_nat[pred_nat$presence == 1 & pred_nat$foldID != f, fold_p_col], probs = 0.025, na.rm = TRUE)  # pred for minimum 2.5% presence in the training data (outside fold)
    glm_pred_nat_rst <- rast(paste0("../pred_rastersBIAS/native/", species, "_nat_p_glm.tif"))
    bart_pred_nat_rst <- rast(paste0("../pred_rastersBIAS/native/", species, "_nat_p_bart.tif"))
    
    proc <- kuenm_proc(occ.test = fold_dat[fold_dat$presence == 1, c("lon", "lat")], model = raster(get(paste0(m, "_pred_nat_rst"))), threshold = 2.5, rand.percent = 50, iterations = 1000)$pROC_summary
    threshmeas <- threshMeasures(obs = fold_dat$presence, pred = fold_dat[ , fold_p_col], measures = c("Omission", "Sensitivity", "Specificity"), thresh = thresh, plot = FALSE)$ThreshMeasures[,1]
    
    crossval[f, "n_pres"] <- sum(pred_nat[pred_nat$foldID == f, "presence"] == 1, na.rm = TRUE)
    crossval[f, "n_abs"] <- sum(pred_nat[pred_nat$foldID == f, "presence"] == 0, na.rm = TRUE)
    crossval[f, paste0("thresh", m)] <- thresh
    crossval[f, paste0(c("AUCratio_", "pROCpval_"), m)] <- proc
    crossval[f, paste0(c("omis_", "sens_", "spec_"), m)] <- threshmeas
  }
  temp <- data.frame(fold = as.integer(rownames(crossval)), crossval)
  write.csv(temp, paste0("../eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics/crossval_metrics_nat_", species, ".csv"), row.names = FALSE)
  
  # boxplots of the cross-validation metrics:
  jpeg(paste0("../eval_metrics/NATIVE/BIAS/2.5E/crossval_boxplots/", species, ".jpg"))
  par(mar = c(7.2, 3, 2, 1))
  boxplot(crossval[ , 3:ncol(crossval)], las = 2, main = species)
  abline(h = 1, col = "darkgrey", lty = 2)
  dev.off()
  
  print(species)
  print(crossval)
  gc()
}

remove(cl)
