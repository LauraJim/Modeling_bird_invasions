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

#################
###IMPORT DATA###
#################

##occurrence data
native.occs <- readOGR("E:/SDM/occurrences","native.range")
invasive.occs <- readOGR("E:/SDM/occurrences","invasive.range")

##for accessible areas
europe.default <- readOGR("E:/SDM/biogeo","europe")
europe.invashist <- readOGR("E:/SDM/biogeo","invasive_range_acces_area")
bioregs <- readOGR("E:/SDM/biogeo","newRealms")

#climate grids (5mins)
datafiles <- list.files(path="E:/SDM/worldclimv2", 
                        pattern =".tif$", full.names=TRUE)
  ClimGrids <- stack()
  for(i in 1:NROW(datafiles)){
    tempraster <- raster(datafiles[i])
    ClimGrids <- stack(ClimGrids,tempraster)
  }

#########################  
###EXAMPLE WITH esttro### (Estrilda troglodytes)
######################### 
  
######################  
###DATA PREPARATION###  
######################  
  
#select esttro
native.bird <- native.occs[native.occs$species=="esttro",]
  nrow(native.bird)
invasive.bird <- invasive.occs[invasive.occs$species=="esttro",]  
  nrow(invasive.bird)  

#restrict GBIF data to occurrences in (close to native range boundaries)
native.range.shp <- readOGR("E:/SDM/RANGES/esttro","esttro")
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
ClimGridsPCA <- rasterPCA(ClimGrids.masked,nComp=2,spca=TRUE)
  plot(ClimGridsPCA$map)
  summary(ClimGridsPCA$model)
  ClimGridsPCA <- stack(ClimGridsPCA$map$PC1,ClimGridsPCA$map$PC2)
  plot(ClimGridsPCA)
  
  ClimGridsPCA.native <- mask(ClimGridsPCA,bioreg.birds)
    plot(ClimGridsPCA.native)
  ClimGridsPCA.invasive <- mask(ClimGridsPCA,europe.default)
    plot(ClimGridsPCA.invasive)

#rarefy occurrrence data 
    native.bird.rare <- humboldt.occ.rarefy(in.pts=native.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km")
    invasive.bird.rare <- humboldt.occ.rarefy(in.pts=invasive.bird,colxy=1:2, rarefy.dist = 50, rarefy.units = "km") 
  
#################################  
###SOME DISTRIBUTION MODELLING###
#################################

###        
#some more data formatting and preparations
###
    
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
  
#create folds for 5-fold cross validation
n.folds <- 5
species.data <- data.frame(fold(data =species.data, k = n.folds, method = "n_rand"))
  head(species.data)
  
###
#SDM: BRT  
###
technique <- "brt"
  
cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
#run 'sdm' framework with BRT separately for each fold/partition
brt.temp.results<- foreach (fold.iter= 1:n.folds,.packages=c("sdm","modEvA","parallel")) %dopar% {             
    tempdata <- species.data[!species.data$.fold ==fold.iter,]
    nrow(tempdata[tempdata$species ==1,])
    nrow(tempdata[tempdata$species ==0,])
    
    d <- sdmData(species~PC1+PC2,train=tempdata[tempdata$species ==1,][c(1:5)],
                 bg=tempdata[tempdata$species ==0,][c(1:5)])
    m <- sdm(species~PC1+PC2,data=d,methods=c('brt'),
             parallelSettings=list(ncore=floor(cores/n.folds), method='parallel'),
             modelSettings=list(brt=list(n.trees=10000,distribution ="bernoulli",interaction.depth = 5)))
    
    m.predict.native <- predict(m,ClimGridsPCA.native,parallelSettings=list(ncore=floor(cores/2*n.folds), method='parallel'))
    m.predict.invasive <- predict(m,ClimGridsPCA.invasive,parallelSettings=list(ncore=floor(cores/2*n.folds), method='parallel'))
    
    response.curves <- data.frame(getResponseCurve(m,1)@response$PC1,getResponseCurve(m,1)@response$PC2)
    
    variable.importance <- cbind(data.frame(getVarImp(m)@varImportance),rep(fold.iter,2))
    
    results <- list(m.predict.native,response.curves,variable.importance,m.predict.invasive)
  }
  stopCluster(cl)
  rm(cl) 
  
#extract spatial prediction of HS across native range  
  brt.native.preds <- crop(stack(brt.temp.results[[1]][[1]],brt.temp.results[[2]][[1]],brt.temp.results[[3]][[1]],
                            brt.temp.results[[4]][[1]],brt.temp.results[[5]][[1]]),bioreg.birds)
    plot(brt.native.preds)
    #just for the fun of visualization, a simple unweighted mean and the occurrences used
    plot(mean(brt.native.preds))
    plot(SpatialPoints(coords = cbind(native.bird.rare$x,native.bird.rare$y)),add=TRUE)
  
  brt.invasive.preds <- crop(stack(brt.temp.results[[1]][[4]],brt.temp.results[[2]][[4]],brt.temp.results[[3]][[4]],
                            brt.temp.results[[4]][[4]],brt.temp.results[[5]][[4]]),europe.default)
    plot(brt.invasive.preds)
    #just for the fun of visualization, a simple unweighted mean and the occurrences used
    plot(mean(brt.invasive.preds))
    plot(SpatialPoints(coords = cbind(invasive.bird.rare$x,invasive.bird.rare$y)),add=TRUE)
  
###  
#check response curves
###
track <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
order <- c(seq(1,100),seq(1,100),seq(1,100),seq(1,100),seq(1,100))
  
response.curves <- rbind(cbind(brt.temp.results[[1]][[2]]),
                         cbind(brt.temp.results[[2]][[2]]),
                         cbind(brt.temp.results[[3]][[2]]),
                         cbind(brt.temp.results[[4]][[2]]),
                         cbind(brt.temp.results[[5]][[2]]))
  response.curves$track <- track
  response.curves$order <- order
  head(response.curves)
  
response.curves.summary <- ddply(response.curves, c("order"), summarise,
                                 PC1    = mean(PC1),
                                 PC1hs = mean(brt_ID.1),
                                 PC1hs.sd = sd(brt_ID.1),
                                 PC2 = mean(PC2),
                                 PC2hs = mean(brt_ID.1.1),
                                 PC2hs.sd = sd(brt_ID.1.1))

  response.curves.summary$lowPC1 <- response.curves.summary$PC1hs-qnorm(.975)*(response.curves.summary$PC1hs.sd/sqrt(5))
  response.curves.summary$highPC1 <- response.curves.summary$PC1hs+qnorm(.975)*(response.curves.summary$PC1hs.sd/sqrt(5))
  
  response.curves.summary$lowPC2 <- response.curves.summary$PC2hs-qnorm(.975)*(response.curves.summary$PC2hs.sd/sqrt(5))
  response.curves.summary$highPC2 <- response.curves.summary$PC2hs+qnorm(.975)*(response.curves.summary$PC2hs.sd/sqrt(5))
  head(response.curves.summary)
  
pPC1 <- plot.response <- ggplot(response.curves.summary, aes(x = PC1, y = PC1hs, group = 1)) + 
    geom_line(col='red') + theme_minimal() +
    geom_ribbon(aes(ymin = lowPC1, ymax = highPC1), alpha = 0.1) + ylim(0,1) + xlim(minValue(ClimGridsPCA$PC1),maxValue(ClimGridsPCA$PC1))+
    annotate("rect", xmin = minValue(ClimGridsPCA.native$PC1), xmax = maxValue(ClimGridsPCA.native$PC1), ymin = 0.8, ymax = 1, alpha = .3,fill = "red") + 
    annotate("rect", xmin = minValue(ClimGridsPCA.invasive$PC1), xmax = maxValue(ClimGridsPCA.invasive$PC1), ymin = 0.85, ymax = 0.95, alpha = .6,fill = "black")
  pPC1
  
pPC2 <- plot.response <- ggplot(response.curves.summary, aes(x = PC2, y = PC2hs, group = 1)) + 
    geom_line(col='red') + theme_minimal() +
    geom_ribbon(aes(ymin = lowPC2, ymax = highPC2), alpha = 0.1) + ylim(0,1)+ xlim(minValue(ClimGridsPCA$PC2),maxValue(ClimGridsPCA$PC2))+
    annotate("rect", xmin = minValue(ClimGridsPCA.native$PC2), xmax = maxValue(ClimGridsPCA.native$PC2), ymin = 0.8, ymax = 1, alpha = .3,fill = "red") + 
    annotate("rect", xmin = minValue(ClimGridsPCA.invasive$PC2), xmax = maxValue(ClimGridsPCA.invasive$PC2), ymin = 0.85, ymax = 0.95, alpha = .6,fill = "black")
pPC2

###
#check variable importance
###
    var.importance <- rbind(brt.temp.results[[1]][[3]],brt.temp.results[[2]][[3]],brt.temp.results[[3]][[3]],
                          brt.temp.results[[4]][[3]],brt.temp.results[[5]][[3]])
  
  var.importance.summary <- ddply(var.importance, c("variables"), summarise,
                                  scorTest    = mean(corTest),
                                  scorTest.sd = sd(corTest),
                                  sAUCtest = mean(AUCtest),
                                  sAUCtest.sd = sd(AUCtest))
  var.importance.summary                                
  
### 
#model evaluation for the native range  
###
cores <- detectCores()
  cl <- makeCluster(n.folds)
  registerDoParallel(cl)
  
for.native.eval<- foreach (fold.iter= 1:n.folds,.packages=c("modEvA","ntbox","raster","ecospat")) %dopar% {
    
    #take the data in the current fold and extract HS values from the corresponding native range prediction
    tempdata.indep <- species.data[species.data$.fold ==fold.iter,]
    nrow(tempdata.indep[tempdata.indep$species ==1,]) #check number of occurrences in first fold
    nrow(tempdata.indep[tempdata.indep$species ==0,]) #check numnber of bg data in first fold
    PresBackg.hs <- extract(brt.native.preds[[fold.iter]],tempdata.indep[c(4:5)])
    PresBackg.hs <- data.frame(PresBackg.hs,tempdata.indep$species)
    colnames(PresBackg.hs) <- c("hs","species")
    
    #AUC
    AUC <- AUC(obs=PresBackg.hs$species,pred=PresBackg.hs$hs)$AUC
    
    #TSS (max)
    TSS <- ecospat.max.tss(PresBackg.hs$hs,PresBackg.hs$species)$max.TSS
    
    #SENSITIVITY (at max TSS)???
    temp1 <- AUC(obs=PresBackg.hs$species,pred=PresBackg.hs$hs)
    temp2 <- ecospat.max.tss(PresBackg.hs$hs,PresBackg.hs$species)$max.threshold
    SENS <- temp1$thresholds[temp1$thresholds$thresholds==temp2,]$sensitivity
    
    #SPECIFICITY (at max TSS)
    SPEC <- temp1$thresholds[temp1$thresholds$thresholds==temp2,]$specificity
    
    #pROC (5% omission)
    partial_roc5 <- pROC(continuous_mod=brt.native.preds[[fold.iter]],
                         test_data = tempdata.indep[tempdata.indep$species ==1,][c(4:5)],
                         n_iter=1000,E_percent=5,
                         boost_percent=50,
                         parallel=FALSE)
    partial_roc5 <- data.frame(t(partial_roc5$pROC_summary))
    colnames(partial_roc5) <- c("AUCat5","AUCratioat5","Pat5")
    
    eval <- data.frame(AUC,TSS,SENS,SPEC,partial_roc5,technique,fold.iter)
    eval
    
  }    
  
  stopCluster(cl)
  rm(cl)      
  
  native.eval <- rbind(data.frame(t(unlist(for.native.eval[1]))),data.frame(t(unlist(for.native.eval[2]))),
                       data.frame(t(unlist(for.native.eval[3]))),data.frame(t(unlist(for.native.eval[4]))),
                       data.frame(t(unlist(for.native.eval[5]))))
  native.eval
  
###
#model evaluation for the native range [europe.default as background]
###
  brt.invasive.preds
  
#rarefied invasive range occurrence data WITH ENVIRONMENTAL DATA    
  invasive.bird.rare.vars <-  data.frame(extract(ClimGridsPCA.invasive,invasive.bird.rare[c(1:2)]))
  invasive.bird.rare.vars$species <- rep(1,nrow(invasive.bird.rare.vars))
  invasive.bird.rare.vars <- na.omit(data.frame(invasive.bird.rare.vars,invasive.bird.rare[c(1:2)]))
  head(invasive.bird.rare.vars)
  nrow(invasive.bird.rare.vars)
  
#select background data WITH ENVIRONMENTAL DATA 
  n.bg <- 1000    #number of background points to be selected randomly but excluding presence locations    
  bg.data <- data.frame(randomPoints(ClimGridsPCA.invasive,n.bg,p=invasive.bird,excludep=TRUE))
  bg.data.vars <- data.frame(extract(ClimGridsPCA.invasive,bg.data))
  bg.data.vars$species <- rep(0,nrow(bg.data.vars))
  bg.data.vars <- na.omit(data.frame(bg.data.vars,bg.data))
  head(bg.data.vars)
  nrow(bg.data.vars)
  
  species.invasive.data <- rbind(invasive.bird.rare.vars,bg.data.vars)    
  head(species.invasive.data)

#quick evaluation of predictive accuracy based on the unweighted mean of the 5 folds
      
    #extract hs values corresponding to invasive occurrences and pseudo-absences
    PresBackg.hs <- extract(mean(brt.invasive.preds),species.invasive.data[c(4:5)])
    PresBackg.hs <- data.frame(PresBackg.hs,species.invasive.data$species)
    colnames(PresBackg.hs) <- c("hs","species")
    
    #find minimum training presence threshold across the native range
    mtpt <- extract(mean(brt.native.preds),native.bird.rare[c(1:2)])  #minimum training presence I think is a bit extreme but let's give it a go
      boxplot(mtpt)
      mtpt.min <- round(min(mtpt,na.rm=TRUE),2)
      mtpt.min
      
    #run AUC to get list of statistics
    temp1 <- AUC(obs=PresBackg.hs$species,pred=PresBackg.hs$hs)
    
    #get statistics
    SENS <- temp1$thresholds[temp1$thresholds$thresholds==mtpt.min,]$sensitivity
    SPEC <- temp1$thresholds[temp1$thresholds$thresholds==mtpt.min,]$specificity  
    TSS <- SENS + SPEC +-1
    invasive.eval <- cbind(SENS, SPEC, TSS)
    invasive.eval
    
    #make map
    mreclas <- c(0, mtpt.min,0,  mtpt.min,1, 1)
    rclmat <- matrix(mreclas, ncol=3, byrow=TRUE)
    rclmat
    
    InvasPresAbs <- reclassify(mean(brt.invasive.preds),rclmat)
      plot(InvasPresAbs)
      plot(SpatialPoints(coords = cbind(invasive.bird.rare$x,invasive.bird.rare$y)),add=TRUE)
    