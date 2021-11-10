##libraries##
library(rgdal)
library(raster)
library(dismo)
#library(sdm)
library(ggplot2)
library(ecospat)

# Read species IDs
bird.species <- list.dirs(path = "./ranges", full.names = FALSE, recursive = FALSE)

##occurrence data
native.occs <- readOGR("./occurrences","native.range")
invasive.occs <- readOGR("./occurrences","invasive.range")

## relevant regions
europe.default <- readOGR("./biogeo","europe")
bioregs <- readOGR("./biogeo","newRealms")

# summary information
sp.nocc.ini <- vector("numeric",length = length(bird.species))
sp.nocc.end <- sp.nocc.ini
sp.nocc.rare <- sp.nocc.ini
sp.nocc.inv <- sp.nocc.ini
sp.nocc.inv.rare <- sp.nocc.ini

for (x in 1:length(bird.species)){ 
  x=2
  species.id <- bird.species[x]
  
# 1) select species occurrences from GBIF dataset
  native.bird <- native.occs[native.occs$species==species.id,]
  # initial sample size
  sp.nocc.ini[x] <- nrow(native.bird)
  invasive.bird <- invasive.occs[invasive.occs$species==species.id,]  
  sp.nocc.inv[x] <- nrow(invasive.bird)  
  
# 2) restrict GBIF data to occurrences in (close to native range boundaries)
  native.range.shp <- readOGR(paste0("./ranges/",species.id),species.id)
  
  #http://datazone.birdlife.org/species/spcdistPOS
  #subset BirdLife shape files
  #keep native and reintroduced ranges
  native.range.shp <- native.range.shp[which(native.range.shp$ORIGIN == 1 |native.range.shp$ORIGIN == 2),]  
  #keep resident and breeding season ranges
  native.range.shp <- native.range.shp[which(native.range.shp$SEASONAL == 1 |native.range.shp$SEASONAL == 2),]
  #buffer with 0.5?C
  native.range.shp <- buffer(native.range.shp,width=0.5,dissolve=TRUE)                                          
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
  pc.native <- raster(paste0("./bioreg_climgrids/",paste0(bird.species[x],"/"),
                             bird.species[x],"_ClimGridsPCA.native.tif"))
  pc.invasive <- raster(paste0("./bioreg_climgrids/",paste0(bird.species[x],"/"),
                               bird.species[x],"_ClimGridsPCA.invasive.tif"))

# 5) Add environmental data to rarefied occurrences    
  native.bird.rare.vars <- data.frame(extract(pc.native,native.bird.rare))
  invasive.bird.rare.vars <- data.frame(extract(pc.invasive,invasive.bird.rare))
  # paste data frames (longitude, latitude, PC1, PC2)
  native.bird.rare.vars <- na.omit(data.frame(native.bird.rare.ll,native.bird.rare.vars))
  invasive.bird.rare.vars <- na.omit(data.frame(invasive.bird.rare.ll,invasive.bird.rare.vars))
  head(native.bird.rare.vars)
  dim(native.bird.rare.vars)
  
    #select custom native range background: bioregions intersecting with occurrences
  bioreg.birds <- raster::intersect(bioregs,native.bird.BL)
  #plot(bioreg.birds)
  #plot(native.range.shp,add=TRUE,col="blue")
  #plot(native.bird.BL,add=TRUE,col="red")
  
  #select background data WITH ENVIRONMENTAL DATA 
  n.bg <- 1000    #number of background points to be selected randomly but excluding presence locations    
  bg.data <- data.frame(randomPoints(ClimGridsPCA.native,n.bg,p=native.bird,excludep=TRUE))
  bg.data.vars <- data.frame(extract(ClimGridsPCA.native,bg.data))
  bg.data.vars$species <- rep(0,nrow(bg.data.vars))
  bg.data.vars <- na.omit(data.frame(bg.data.vars,bg.data))
  head(bg.data.vars)
  nrow(bg.data.vars)
  
  species.data <- rbind(native.bird.rare.vars,bg.data.vars)    
  head(species.data)
  
}

# Save summary information
sum.table <- cbind(bird.species,sp.nocc.ini,sp.nocc.end,sp.nocc.rare,
                   sp.nocc.invasive,sp.nocc.inv.rare)
colnames(sum.table) <- c("speciesID","N.initial","N.cleaned","N.rarefied",
                         "N.invasive","N.inv.rarefied")
### END