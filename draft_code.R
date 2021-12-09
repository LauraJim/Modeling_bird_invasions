source(".\\Nf_modeling\\GEspace.R")
source(".\\Nf_modeling\\nicheG.R")

# read table of a species occurrence with environmental data points and table 
#  with random background points that contain environmental data
species <- read.csv("./Nf_modeling/occurrences-v2/plomel_native.csv",header=T)[,-1]
ranpoints <- read.csv("./Nf_modeling/accessible-areas-v2/plomel_native_range.csv",header=T)[,-1]

# Plot datasets in both Geographic and Environmental spaces
GE.space(bckgrnd=ranpoints, GE.occ=species)


# Project model into the native range
species.id <- "estmel"
pc.native <- stack(paste0("./bioreg_climgrids/",paste0(species.id,"/"),
                          species.id,"_ClimGridsPCA.native.tif"))

estmel.wn <- niche.G(Estck = pc.native, mu = c(3.00512808, 7.036996838), 
                 Sigma = matrix(c(3.367508327, -2.28366811, -2.28366811, 9.342821904),ncol=2))

x11()
plot(estmel.wn,xlim=c(-30,60),ylim=c(-35,30))

writeRaster(estmel.wn,"./Nf_modeling/Results/estmel_suitability_map.tif", overwrite = T)
