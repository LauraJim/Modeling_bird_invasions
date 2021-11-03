##libraries##
library(rgdal)
library(raster)
library(RStoolbox)
library(humboldt)
library(dismo)
library(sdm)
library(groupdata2)
library(modEvA)
library(parallel)
library(foreach)
library(plyr)
library(ggplot2)
library(ecospat)
library(ntbox)
library(doParallel)
library(sp)

bird.species <- list.dirs(path = "C:/Users/ststrubb/Dropbox/SDM/ranges", full.names = FALSE, recursive = FALSE)

for (x in bird.species){ 
  species.id <- x

##occurrence data
native.occs <- readOGR("C:/Users/ststrubb/Dropbox/SDM/occurrences","native.range")
invasive.occs <- readOGR("C:/Users/ststrubb/Dropbox/SDM/occurrences","invasive.range")

##for accessible areas
europe.default <- readOGR("C:/Users/ststrubb/Dropbox/SDM/biogeo","europe")
bioregs <- readOGR("C:/Users/ststrubb/Dropbox/SDM/biogeo","newRealms")

#climate grids (5mins)
datafiles <- list.files(path="C:/Users/ststrubb/Dropbox/SDM/worldclimv2", 
                        pattern =".tif$", full.names=TRUE)
  ClimGrids <- stack()
  for(i in 1:NROW(datafiles)){
    tempraster <- raster(datafiles[i])
    ClimGrids <- stack(ClimGrids,tempraster)
  }

#select species.id
native.bird <- native.occs[native.occs$species==species.id,]
  nrow(native.bird)
invasive.bird <- invasive.occs[invasive.occs$species==species.id,]  
  nrow(invasive.bird)  

#restrict GBIF data to occurrences in (close to native range boundaries)
native.range.shp <- readOGR(paste0("C:/Users/ststrubb/Dropbox/SDM/RANGES/",species.id),species.id)

  #http://datazone.birdlife.org/species/spcdistPOS
  #subset BirdLife shape files
  native.range.shp <- native.range.shp[which(native.range.shp$ORIGIN == 1 |native.range.shp$ORIGIN == 2),]      #keep native and reintroduced ranges
  native.range.shp <- native.range.shp[which(native.range.shp$SEASONAL == 1 |native.range.shp$SEASONAL == 2),]  #keep resident and breeding season ranges
  native.range.shp <- buffer(native.range.shp,width=0.5,dissolve=TRUE)                                          #buffer with 0.5°C
  plot(native.range.shp)
  
  native.bird.BL <- raster::intersect(native.bird,native.range.shp)                                             #keep only occurrences that are within the buffered native range
      nrow(native.bird.BL)

#select custom native range background: bioregions intersecting with occurrences
bioreg.birds <- raster::intersect(bioregs,native.bird.BL)
      plot(bioreg.birds)
      plot(native.range.shp,add=TRUE,col="blue")
      plot(native.bird.BL,add=TRUE,col="red")
      
#create predictor grids     
cutout <- bind(bioreg.birds,europe.default)  
ClimGrids.masked <- mask(ClimGrids,cutout)  

#do PCA on all EU + native range climate variables and select first 2 axes  
ClimGridsPCA.all <- rasterPCA(ClimGrids.masked,nComp=2,spca=TRUE)
  plot(ClimGridsPCA.all$map)
  
  #PCA loadings
  loadings <- data.frame(round(ClimGridsPCA.all$model$loadings[,1:2],2)) # top 2 loadings
  loadings$variables <- row.names(loadings)
  
  #PCA VarImportance
  x <- summary(ClimGridsPCA.all$model)
    vars <- x$sdev^2
    vars <- vars/sum(vars)
    VarImportance <- data.frame(rbind("Standard deviation" = x$sdev, "Proportion of Variance" = vars,"Cumulative Proportion" = cumsum(vars))) 
    VarImportance$indices <- row.names(VarImportance)
    VarImportance <- VarImportance[c(ncol(VarImportance),1:(ncol(VarImportance)-1))]
             
  #create native and invasive range grids
  ClimGridsPCA <- stack(ClimGridsPCA.all$map$PC1,ClimGridsPCA.all$map$PC2)
    plot(ClimGridsPCA)
  
  ClimGridsPCA.native <- mask(ClimGridsPCA,bioreg.birds)
  ClimGridsPCA.invasive <- mask(ClimGridsPCA,europe.default)
  
setwd("C:/Users/ststrubb/Dropbox/SDM/bioreg_climgrids")
  dir.create(species.id)
  write.table(loadings, file.path(paste0("./",species.id), paste0(species.id,"_loadings",".txt")), sep = "\t",row.names=FALSE)
  write.table(VarImportance, file.path(paste0("./",species.id), paste0(species.id,"_VarImportance",".txt")), sep = "\t",row.names=FALSE)
  
  writeRaster(ClimGridsPCA.native, file.path(paste0("./",species.id), paste0(species.id,"_ClimGridsPCA.native",".tif")),overwrite=TRUE)
  writeRaster(ClimGridsPCA.invasive, file.path(paste0("./",species.id), paste0(species.id,"_ClimGridsPCA.invasive",".tif")),overwrite=TRUE)
rm(list = ls())
}