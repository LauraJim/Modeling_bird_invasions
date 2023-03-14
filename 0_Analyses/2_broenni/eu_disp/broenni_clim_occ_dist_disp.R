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
library(adehabitatHR)
library(rgeos)
library(geobuffer)

bird.species <- list.dirs(path = "E:/SDM/ranges", full.names = FALSE, recursive = FALSE)
#bird.species <- bird.species[c(1:3)]

niche.results <- data.frame(matrix(vector(), 0, 10,
                                        dimnames=list(c(), c("sim1", "sim2", "equ","niche_overlap_F","niche_overlap_T","stability.index","expansion.index","unfilling.index","species","threshold"))),
                                 stringsAsFactors=F)

thresholds <- c(0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1)

for (x in bird.species){ 
    species.id <- x

############################################################################  
###CREATE SHAPEFILES FOR DELINEATING NATIVE AND INVASIVE BACKGROUND AREAS###
############################################################################  
    regions <- readOGR("E:/SDM/bioreg_occ_dist_climgrids/backgrounds",species.id)
    #europe <- readOGR("E:/SDM/biogeo","europe")
    
    ##occurrence data
    native.occs <- readOGR("E:/SDM/occurrences","native.range")
    invasive.occs <- readOGR("E:/SDM/occurrences","invasive.range")
    
    #climate grids (5mins)
    datafiles <- list.files(path="E:/SDM/worldclimv2", 
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
native.range.shp <- readOGR(paste0("E:/SDM/ranges/",species.id),species.id)

#http://datazone.birdlife.org/species/spcdistPOS
#subset BirdLife shape files
native.range.shp <- native.range.shp[which(native.range.shp$ORIGIN == 1 |native.range.shp$ORIGIN == 2),]      #keep native and reintroduced ranges
  native.range.shp <- native.range.shp[which(native.range.shp$SEASONAL == 1 |native.range.shp$SEASONAL == 2),]  #keep resident and breeding season ranges
  native.range.shp <- gBuffer(native.range.shp,width=0.5)                                          #buffer with 0.5°C
  plot(native.range.shp)
  
native.bird.BL <- raster::intersect(native.bird,native.range.shp)                                             #keep only occurrences that are within the buffered native range
      nrow(native.bird.BL)
      native.bird.BL.r <- sp::remove.duplicates(native.bird.BL, zero = 50)
      native.bird.BL.r <- humboldt.occ.rarefy(in.pts=native.bird.BL.r,colxy=1:2, rarefy.dist = 50, rarefy.units = "km")
      class(native.bird.BL.r)
      head(native.bird.BL.r)
      plot(native.bird.BL.r$x,native.bird.BL.r$y)
      
#create occurence distance around native range occurrence  
     gcDist <- sp::spDists(as.matrix(native.bird.BL.r[c(1:2)]),longlat = FALSE)
     gcDist[gcDist == 0] <- NA 
     buf_distance <- mean(gcDist,na.rm=TRUE)
     buf_distance
     
     native.bird.BL.r <- SpatialPoints(native.bird.BL.r[c(1:2)])
     native.bird.BL.occ_dist <- gBuffer(native.bird.BL.r,width=buf_distance)
     
     plot(regions)
     plot(native.range.shp,add=TRUE,col='blue')
     plot(native.bird.BL.occ_dist,add=TRUE,col='red')
     plot(native.bird.BL,add=TRUE)
     
     bioreg_occ_dist <- raster::intersect(regions,native.bird.BL.occ_dist)
     plot(regions)
     plot(bioreg_occ_dist,add=TRUE,col='red')
     plot(native.bird.BL,add=TRUE,col='yellow')

######################################################        
###CREATE NATIVE AND INVASIVE RANGE PREDICTOR GRIDS###
######################################################
invasive.range.shp <- readOGR(paste("E:/SDM/ranges/",species.id,sep=""),paste("for_",species.id,"_europe_inv",sep=""))

#create predictor grids     
cutout <- bind(bioreg_occ_dist,invasive.range.shp)  
ClimGrids.masked <- mask(ClimGrids,cutout)  

ClimGrids.native <- mask(ClimGrids.masked,native.bird.BL.occ_dist)
ClimGrids.invasive <- mask(ClimGrids.masked,invasive.range.shp)

#create species occurrence climate data
native.bird.occ <- na.exclude(raster::extract(ClimGrids.native,native.bird[c(1:2)]))
  head(native.bird.occ)

invasive.bird.occ <- na.exclude(raster::extract(ClimGrids.invasive,invasive.bird[c(1:2)]))
  head(invasive.bird.occ)
  
#create background climate data
native.bird.bkg <- rasterToPoints(ClimGrids.native)[,c(3:21)]
  head(native.bird.bkg)
  
invasive.bird.bkg <- rasterToPoints(ClimGrids.invasive)[,c(3:21)]
  head(invasive.bird.bkg)  

##############################  
####niche quantifications#####
##############################  
  
  #calibration of PCA-env 
  pca.env <-dudi.pca(rbind(native.bird.bkg,invasive.bird.bkg), center = T, scale = T, scannf = F, nf = 2)
  
  # predict the scores on the PCA axes
  scores.bkg<- pca.env$li
  scores.bkg.native<- suprow(pca.env,native.bird.bkg)$lisup
  scores.bkg.invasive<- suprow(pca.env,invasive.bird.bkg)$lisup
  scores.occ.native<- suprow(pca.env,native.bird.occ)$lisup
  scores.occ.invasive<- suprow(pca.env,invasive.bird.occ)$lisup
  
  thresholds <- seq(0,1,0.05)
  for (x in thresholds){
    threshold <- x
  
  # calculation of occurrence density
  z1<- ecospat.grid.clim.dyn(scores.bkg,scores.bkg.native,scores.occ.native,R=100)
  z2<- ecospat.grid.clim.dyn(scores.bkg,scores.bkg.invasive,scores.occ.invasive,R=100)
  equ<-ecospat.niche.equivalency.test(z1,z2,rep=100,ncores=30)
  
  # test of niche equivalency and similarity
  sim1<-ecospat.niche.similarity.test(z1,z2,rep=100,alternative = "greater",rand.type = 1) #niches randomly shifted in both area
  sim2<-ecospat.niche.similarity.test(z1,z2,rep=100,alternative = "greater",rand.type = 2) #niche randomly shifted only in invaded area
  
  # overlap corrected by availabilty of background conditions
  niche_overlap_F <- ecospat.niche.overlap(z1,z2,cor=T) 
  # uncorrected overlap
  niche_overlap_T <- ecospat.niche.overlap(z1,z2,cor=F) 
  
################################
#### niche visualizations  #####
################################
  
  # occurrence density plots
  #setwd("E:/SDMbroenni")
  #dir.create(species.id)
  
  #mypath <- file.path("E:/SDMbroenni",species.id,paste("native range occurrende density", species.id, ".png", sep = ""))
  #png(file=mypath, width = 225, height = 225, units='mm', res = 300)
  #ecospat.plot.niche(z1,title="PCA-env - EU niche",name.axis1="PC1",name.axis2="PC2")
  #dev.off()
  
  #mypath <- file.path("E:/SDMbroenni",species.id,paste("invasive range occurrence density", species.id, ".png", sep = ""))
  #png(file=mypath, width = 225, height = 225, units='mm', res = 300)
  #ecospat.plot.niche(z2,title="PCA-env - NA niche",name.axis1="PC1",name.axis2="PC2")
  #dev.off()
  
  # contribution of original variables
  #mypath <- file.path("E:/SDMbroenni",species.id,paste("PCA contributions", species.id, ".png", sep = ""))
  #png(file=mypath, width = 225, height = 225, units='mm', res = 300)
  #ecospat.plot.contrib(pca.env$co,pca.env$eig)
  #dev.off()
  
  # niche tests plots
  #ecospat.plot.overlap.test(equ,"D","Equivalency")
  #ecospat.plot.overlap.test(sim1,"D","Similarity 1<->2") #niches randomly shifted in both areas
  #ecospat.plot.overlap.test(sim2,"D","Similarity 1->2") #niche randomly shifted only in invaded area
  
  #mypath <- file.path("E:/SDMbroenni",species.id,paste("niche overlap native range density", species.id, ".png", sep = ""))
  #png(file=mypath, width = 225, height = 225, units='mm', res = 300)
  #ecospat.plot.niche.dyn(z1,z2,quant=0.5,name.axis1="PC1",name.axis2="PC2",interest=1)
  #dev.off()
  
  #mypath <- file.path("E:/SDMbroenni",species.id,paste("niche overlap invasive range density", species.id, ".png", sep = ""))
  #png(file=mypath, width = 225, height = 225, units='mm', res = 300)
  #ecospat.plot.niche.dyn(z1,z2,quant=0.5,name.axis1="PC1",name.axis2="PC2",interest=2)
  #dev.off()
  
  #niche dynamics
  ind1 <- ecospat.niche.dyn.index(z1, z2,intersection=threshold)
  
  expansion.index <- ind1$dynamic.index.w[1]
  stability.index <- ind1$dynamic.index.w[2]
  unfilling.index <-ind1$dynamic.index.w[3]
  
  
  results <- data.frame(t(rbind(sim1$p.D,sim2$p.D,equ$p.D, niche_overlap_F$D, niche_overlap_T$D, stability.index, expansion.index, unfilling.index,species.id,threshold)))
  results
  niche.results <- rbind(results, niche.results)
  }
}

write.table(niche.results,"niche_results_th_env.txt",sep="\t",row.names=FALSE)

