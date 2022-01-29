# LOAD PACKAGES ####

library(terra)
library(maps)


# CREATE OUTPUT FOLDER ####

dir.create("../pred_JPGs", recursive = TRUE)


# GET PREDICTION RASTER MAP FILE NAMES ####

pred_files_nat <- list.files("../pred_rasters/native", full.names = TRUE)
pred_files_nat <- pred_files_nat[grep("_f_", pred_files_nat)]
pred_files_nat

pred_files_inv <- list.files("../pred_rasters/invasive", full.names = TRUE)
pred_files_inv <- pred_files_inv[grep("_f_", pred_files_inv)]
pred_files_inv


# IMPORT EACH RASTER AND SAVE AS JPG ####

for (f in c(pred_files_nat, pred_files_inv)) {
  rst <- rast(f)
  filename <- basename(tools::file_path_sans_ext(f))
  species <- sapply(strsplit(filename, "_"), `[`, 1)
  jpeg(paste0("../pred_JPGs/", filename, ".jpg"), width = 800, height = 430)
  maps::map("world", col = "grey", mar = c(1, 1, 2, 6))
  plot(rst, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(basename(tools::file_path_sans_ext(f)))
  dev.off()
}
