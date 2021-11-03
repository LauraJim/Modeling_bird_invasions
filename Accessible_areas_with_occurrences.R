##libraries##
library(rgdal)
library(raster)
library(RStoolbox)
library(humboldt)
library(dismo)
library(sdm)
library(ggplot2)
library(ecospat)
library(sp)

# Read species IDs
bird.species <- list.dirs(path = "./ranges", full.names = FALSE, recursive = FALSE)

##occurrence data
native.occs <- readOGR("./occurrences","native.range")
#invasive.occs <- readOGR("./occurrences","invasive.range")

## relevant regions
#europe.default <- readOGR("./biogeo","europe")
bioregs <- readOGR("./biogeo","newRealms")

# summary information
sp.nocc.ini <- vector("numeric",length = length(bird.species))
sp.nocc.end <- sp.nocc.ini

for (x in 1:length(bird.species)){ 
  species.id <- bird.species[1]
  
  # select species occurrences
  native.bird <- native.occs[native.occs$species==species.id,]
  sp.nocc.ini[1] <- nrow(native.bird)
  #invasive.bird <- invasive.occs[invasive.occs$species==species.id,]  
  #nrow(invasive.bird)  
  
  #restrict GBIF data to occurrences in (close to native range boundaries)
  native.range.shp <- readOGR(paste0("./ranges/",species.id),species.id)
  
  #http://datazone.birdlife.org/species/spcdistPOS
  #subset BirdLife shape files
  native.range.shp <- native.range.shp[which(native.range.shp$ORIGIN == 1 |native.range.shp$ORIGIN == 2),]      #keep native and reintroduced ranges
  native.range.shp <- native.range.shp[which(native.range.shp$SEASONAL == 1 |native.range.shp$SEASONAL == 2),]  #keep resident and breeding season ranges
  native.range.shp <- buffer(native.range.shp,width=0.5,dissolve=TRUE)                                          #buffer with 0.5?C
  #plot(native.range.shp)
  
  native.bird.BL <- raster::intersect(native.bird,native.range.shp)                                             #keep only occurrences that are within the buffered native range
  sp.nocc.end[1] <- nrow(native.bird.BL)
  
  #select custom native range background: bioregions intersecting with occurrences
  bioreg.birds <- raster::intersect(bioregs,native.bird.BL)
  #plot(bioreg.birds)
  #plot(native.range.shp,add=TRUE,col="blue")
  #plot(native.bird.BL,add=TRUE,col="red")
  
  #rarefy occurrrence data 
  native.bird.rare <- humboldt.occ.rarefy(in.pts=native.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km")
  invasive.bird.rare <- humboldt.occ.rarefy(in.pts=invasive.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km") 
  
  #rarefied native range occurrence data WITH ENVIRONMENTAL DATA    
  native.bird.rare.vars <-  data.frame(extract(ClimGridsPCA.native,native.bird.rare[c(1:2)]))
  native.bird.rare.vars$species <- rep(1,nrow(native.bird.rare.vars))
  native.bird.rare.vars <- na.omit(data.frame(native.bird.rare.vars,native.bird.rare[c(1:2)]))
  head(native.bird.rare.vars)
  nrow(native.bird.rare.vars)
  
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