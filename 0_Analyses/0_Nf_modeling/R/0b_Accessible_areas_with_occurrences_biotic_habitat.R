# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##libraries##
library(rgdal)
library(raster)
library(dismo)
library(ggplot2)
library(ecospat)
library(rgeos)

# Read species IDs
bird.species <- list.dirs(path = "../../../1_data/ranges", full.names = FALSE, recursive = FALSE)

##occurrence data
native.occs <- readOGR("../../../1_data/occurrences","native.range")
invasive.occs <- readOGR("../../../1_data/occurrences","invasive.range")

## relevant regions
europe.default <- readOGR("../../../1_data/biogeo","europe")
bioregs <- readOGR("../../../1_data/biogeo","newRealms")

# summary information
sp.nocc.ini <- vector("numeric",length = length(bird.species))
sp.nocc.end <- sp.nocc.ini
sp.nocc.rare <- sp.nocc.ini
sp.nocc.inv <- sp.nocc.ini
sp.nocc.inv.rare <- sp.nocc.ini


# BEFORE RUNNING: change file paths!!!
# See lines: 73,75,89,91,100,108,117

for (x in 1:length(bird.species)){ 
  species.id <- bird.species[x]
  
# 1) select species occurrences from GBIF dataset
  native.bird <- native.occs[native.occs$species==species.id,]
  # initial sample size
  sp.nocc.ini[x] <- nrow(native.bird)
  invasive.bird <- invasive.occs[invasive.occs$species==species.id,]  
  sp.nocc.inv[x] <- nrow(invasive.bird)  
  
# 2) restrict GBIF data to occurrences within (close to) native range boundaries
  native.range.shp <- readOGR(paste0("../../../1_Data/ranges/",species.id),species.id)
  
  #http://datazone.birdlife.org/species/spcdistPOS
  #subset BirdLife shape files
  #keep native and reintroduced ranges
  native.range.shp <- native.range.shp[which(native.range.shp$ORIGIN == 1 |native.range.shp$ORIGIN == 2),]  
  #keep resident and breeding season ranges
  native.range.shp <- native.range.shp[which(native.range.shp$SEASONAL == 1 |native.range.shp$SEASONAL == 2),]
  #buffer with 0.5?C
  native.range.shp <- gBuffer(native.range.shp,width=0.5)                                          
  #plot(native.range.shp)
  
  #keep only occurrences that are within the buffered native range
  native.bird.BL <- raster::intersect(native.bird,native.range.shp)                                             
  # sample size after restriction
  sp.nocc.end[x] <- nrow(native.bird.BL)
  
# 3) rarefy occurrrence data 
  native.bird.rare <- sp::remove.duplicates(native.bird.BL, zero = 50) # WGS 84 units are meters
  # data frame with longitude and latitude coordinates
  native.bird.rare.ll <- native.bird.rare@data[,1:2] # to get data
  # sample size of thinned dataset
  sp.nocc.rare[x] <- nrow(native.bird.rare.ll)
  invasive.bird.rare <- sp::remove.duplicates(invasive.bird, zero = 50)
  invasive.bird.rare.ll <- invasive.bird.rare@data[,1:2]
  # number of occurrences in the invasive area after thinning
  sp.nocc.inv.rare[x] <- nrow(invasive.bird.rare.ll)
    
  #native.bird.rare <- humboldt.occ.rarefy(in.pts=native.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km")
  #invasive.bird.rare <- humboldt.occ.rarefy(in.pts=invasive.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km") 
  
# 4) Read PCA layers
  pc.native <- stack(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_allgrids/",paste0(species.id,"/"),
                            species.id,"_allGridsPCA.native.tif"))
  pc.invasive <- stack(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_allgrids/",paste0(species.id,"/"),
                              species.id,"_allGridsPCA.invasive.tif"))

# 5) Add environmental data to rarefied occurrences    
  native.bird.rare.vars <- data.frame(extract(pc.native,native.bird.rare))
  invasive.bird.rare.vars <- data.frame(extract(pc.invasive,invasive.bird.rare))
  # paste data frames (longitude, latitude, PC1, PC2)
  native.bird.rare.vars <- na.omit(data.frame(native.bird.rare.ll,native.bird.rare.vars))
  invasive.bird.rare.vars <- na.omit(data.frame(invasive.bird.rare.ll,invasive.bird.rare.vars))
  # adjust column names
  colnames(native.bird.rare.vars) <- c("lon","lat","PC1","PC2")
  colnames(invasive.bird.rare.vars) <- c("lon","lat","PC1","PC2")

# 6) Save tables
  write.csv(native.bird.rare.vars,paste0("../../../0_Analyses/0_Nf_modeling/occurrencesClimBioticHabitat/",
                                         species.id,"_native.csv"))
  write.csv(invasive.bird.rare.vars,paste0("../../../0_Analyses/0_Nf_modeling/occurrencesClimBioticHabitat/",
                                         species.id,"_invasive.csv"))
 
# 7) Select random points inside accessible area and add PC values
  N <- 10000        
  rp.aa <- data.frame(randomPoints(pc.native,N,p=native.bird.rare,excludep=T))
  rp.aa.vars <- data.frame(extract(pc.native,rp.aa))
  rp.aa.vars <- na.omit(data.frame(rp.aa,rp.aa.vars))
  colnames(rp.aa.vars) <- c("lon","lat","PC1","PC2")
  write.csv(rp.aa.vars,paste0("../../../0_Analyses/0_Nf_modeling/accessible_areas_ClimBioticHabitat/",
                                         species.id,"_native_range.csv"))
  
# 8) Create table with random points inside the invaded region
  rp.ir <- data.frame(randomPoints(pc.invasive,N/2,p=invasive.bird.rare,excludep=T))
  rp.ir.vars <- data.frame(extract(pc.invasive,rp.ir))
  rp.ir.vars <- na.omit(data.frame(rp.ir,rp.ir.vars))
  colnames(rp.ir.vars) <- c("lon","lat","PC1","PC2")
  write.csv(rp.ir.vars,paste0("../../../0_Analyses/0_Nf_modeling/accessible_areas_ClimBioticHabitat/",
                              species.id,"_invaded_region.csv"))
} #end of for-loop

# Save summary information
sum.table <- cbind(bird.species,sp.nocc.ini,sp.nocc.end,sp.nocc.rare,
                   sp.nocc.inv,sp.nocc.inv.rare)
colnames(sum.table) <- c("speciesID","N.initial","N.cleaned","N.rarefied",
                         "N.invasive","N.inv.rarefied")
write.csv(sum.table,"../../../0_Analyses/0_Nf_modeling/allspecies_samplesizes_ClimBioticHabitat.csv",row.names = F)


### END