
pred_files <- list.files("../pred_rasters/native", full.names = TRUE)
pred_files <- pred_files[grep("_f.tif", pred_files)]
pred_files

dir.create("../pred_JPGs/")

for (f in pred_files) {
  rst <- rast(f)
  filename <- basename(tools::file_path_sans_ext(f))
  species <- sapply(strsplit(filename, "_"), `[`, 1)
  jpeg(paste0("../pred_JPGs/", filename, ".jpg"), width = 800, height = 600)
  maps::map("world", col = "grey", mar = c(1, 1, 2, 6))
  plot(rst, col = hcl.colors(100), range = c(0, 1), add = TRUE)
  title(basename(tools::file_path_sans_ext(f)))
  dev.off()
}
