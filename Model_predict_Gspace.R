source(".\\Nf_modeling\\GEspace.R")
source(".\\Nf_modeling\\nicheG.R")

# Read table with estimated parameters of niche model
mles <- read.csv("./Nf_modeling/Results-v3/mle_allspecies_v3.csv",header=T)


# Project models to geographic space in the native and invaded regions

for (k in 1:nrow(mles)) {
  # Select species
  species.id <- mles[k,1]
  
  # Project model into the native range
  pc.native <- stack(paste0("./bioreg_occ_dist_climgrids/",paste0(species.id,"/"),
                            species.id,"_ClimGridsPCA.native.tif"))
  
  suit.wn.nat <- niche.G(Estck = pc.native, mu = c(mles[k,3],mles[k,4]), 
                       Sigma = matrix(c(mles[k,5], mles[k,6], mles[k,6], 
                                        mles[k,7]),ncol=2))
  writeRaster(suit.wn.nat,paste0("./Nf_modeling/Results-v3/",species.id,
                             "_suitability_map_nat.tif"), overwrite = T)
  # Project model into the invaded region
  pc.invasive <- stack(paste0("./bioreg_occ_dist_climgrids/",paste0(species.id,"/"),
                              species.id,"_ClimGridsPCA.invasive.tif"))
  
  suit.wn.inv <- niche.G(Estck = pc.invasive, mu = c(mles[k,3],mles[k,4]), 
                     Sigma = matrix(c(mles[k,5], mles[k,6], mles[k,6], 
                                      mles[k,7]),ncol=2))
  
  writeRaster(suit.wn.inv,paste0("./Nf_modeling/Results-v3/",species.id,
                             "_suitability_map_inv.tif"), overwrite = T)
}

# Useful code -----------------------------------
# To visualize the datasets in both Geographic and Environmental spaces

# read table of a species occurrence with environmental data points and table 
#  with random background points that contain environmental data
species <- read.csv("./Nf_modeling/occurrences-v2/plomel_native.csv",header=T)[,-1]
ranpoints <- read.csv("./Nf_modeling/accessible-areas-v2/plomel_native_range.csv",header=T)[,-1]

# Plot datasets in both Geographic and Environmental spaces
GE.space(bckgrnd=ranpoints, GE.occ=species)


