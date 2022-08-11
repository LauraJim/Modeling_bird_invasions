# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(ggplot2)
library(raster)
library(cowplot)
library(ggpubr)
library(rgdal)
library(rgbif)
library(humboldt)
library(ggimage)
library(png)
library(patchwork)

# CREATE OUTPUT FOLDER ####

dir.create("../model_objects/", recursive = TRUE)


# GET TARGET SPECIES ####

occ_files <- list.files("../../Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list <- "estast"
species_list


# COMPUTE PA (PRESENCE/ABSENCE) MODELS ####

  species <- species_list
  message("\n", species, "\npreparing data...")
  
  # import native occurrences and accessible (background) points:
  occ_nat <- read.csv(paste0("../../Nf_modeling/occurrences-v3/", species, "_native.csv"))
  head(occ_nat)
  acc_nat <- read.csv(paste0("../../Nf_modeling/accessible-areas-v3/", species, "_native_range.csv"))
  head(acc_nat)
  
  #get BirdLife data
  #myspecies <- c("Estrilda astrild")
  #gbif_data <- occ_data(scientificName = myspecies, hasCoordinate = TRUE, limit = 100000)
  #for.estast2022 <- na.omit(as.data.frame(gbif_data$data)[c(2:4)])
  #for.estast2022 <- for.estast2022[c(1,3,2)]
  #colnames(for.estast2022) <- c("sp","x","y")
  #for.estast2022$sp <- rep("estast",nrow(for.estast2022))
  #head(for.estast2022)
  #estast2022<-humboldt.occ.rarefy(in.pts=for.estast2022,colxy=2:3, rarefy.dist = 100, rarefy.units = "km")
  #plot(estast2022$x,estast2022$y)
  #write.table(estast2022,"estast2022.txt",sep="\t",row.names=FALSE)
  #estast2022 <- read.table("estast2022.txt",sep="\t",h=T)
  #head(estast2022)
  #plot(estast2022$x,estast2022$y)
  
  #estast2022 <- occ_nat
  #estast2022 <- estast2022[c(2:3)]
  #colnames(estast2022) <- c("x","y")
  estast.points <- SpatialPointsDataFrame(coords = c(occ_nat[,c("lon", "lat")]),
                                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                                          data = occ_nat)
  #get native range BirdLife   
  birdlife.raw<-readOGR("D:/Ststrubb/OneDrive - UGent/Projects/MC_paper/Modeling_bird_invasions/ranges/estast","estast")
  birdlife.raw <- birdlife.raw[birdlife.raw$ORIGIN==1,]
  birdlife <- fortify(birdlife.raw)
  
  estast.points <- spTransform(estast.points, CRS(proj4string(birdlife.raw)))
  estast.points <- raster::intersect(estast.points,birdlife.raw)
  plot(estast.points)
  head(estast.points)
  estast.points.spatV <- vect(estast.points)
  
  # import climgrids:
  climgrids_nat <- rast(paste0("../../bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  #plot(climgrids_nat)
  
  #extract data
  temp1 <- terra::extract(climgrids_nat,estast.points.spatV,xy=TRUE)
  temp1$X <- seq(1:nrow(temp1))
  occ_nat_new <- temp1[c(6,4:5,2:3)]
  colnames(occ_nat_new)<- names(occ_nat)
  head(occ_nat_new)
  
  occ_nat <- na.omit(occ_nat_new)
  
  # get presence (occupied) and absence (accessible not occupied) pixels:
  occ_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(occ_nat[ , c("lon", "lat")]))
  acc_nat_cells <- terra::cellFromXY(object = climgrids_nat, xy = as.matrix(acc_nat[ , c("lon", "lat")]))
  # length(occ_nat_cells)
  # length(acc_nat_cells)
  abs_nat <- acc_nat[!(acc_nat_cells %in% occ_nat_cells), ]
  occ_nat$presence <- 1
  abs_nat$presence <- 0
  pa_nat <- rbind(occ_nat, abs_nat)
  # head(pa_nat)
  # sum(pa_nat$presence == 1)
  # sum(pa_nat$presence == 0)
  
  
  # compute presence-absence models:
  
  message("computing models...")
  #names(pa_nat)
  var_names <- names(pa_nat)[grep("PC", names(pa_nat))]
  #var_names
  
  form_glm <- as.formula(paste("presence ~", paste(var_names, collapse = "+")))
  mod_glm <- glm(formula = form_glm, family = binomial, data = pa_nat)
  #summary(mod_glm)
  saveRDS(mod_glm, paste0("../model_objects/", species, "_mod_glm.rds"))
  
  set.seed(grep(species, species_list))
  mod_bart <- bart(x.train = pa_nat[ , var_names], y.train = pa_nat[ , "presence"], keeptrees = TRUE, verbose = FALSE)
  #summary(mod_bart)
  invisible(mod_bart$fit$state)
  saveRDS(mod_bart, paste0("../model_objects/", species, "_mod_bart.rds"))
  
  gc()
  # end for spc

###################################################################################
###RESPONSE CURVES#################################################################  
###################################################################################  
climgrids_nat <- rast(paste0("../../bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
climgrids_inv <- rast(paste0("../../bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.invasive.tif"))

PC1 <- data.frame(merge(climgrids_nat[[1]],climgrids_inv[[1]]))
PC2 <- data.frame(merge(climgrids_nat[[2]],climgrids_inv[[2]]))

#
###FOR PC1### 
#
forPC1 <- data.frame(seq(min(data.frame(PC1)),max(data.frame(PC1)),by=0.1))
  colnames(forPC1)<- "PC1"
  forPC1$PC2 <- rep(mean(PC2$estast_ClimGridsPCA.native_2))
  head(forPC1)
  #GLM
  glm_p <- data.frame(predict(mod_glm,forPC1, type = "response",se=TRUE))
    forPC1$prob <- glm_p$fit
    forPC1$prob_se <- glm_p$se.fit
    forPC1$method <- rep("GLM",nrow(forPC1))
    head(forPC1)
  #BART
    x1 <- forPC1[c(1:2)]
    x1$x <- seq(1:nrow(x1))
    x1$y <- seq(1:nrow(x1))
    raster.x <- rasterFromXYZ(x1[c(3:4,1)])
    raster.y <- rasterFromXYZ(x1[c(3:4,2)])  
  bart_p <- predict(mod_bart, raster::stack(raster.x,raster.y),quantiles = c(0.05,0.95))
    temp <-  forPC1[c(1:2)]
    temp$bart_prob <- na.omit(as.data.frame(bart_p))$layer.1
    upper <- na.omit(as.data.frame(bart_p))$layer.3
    lower <- na.omit(as.data.frame(bart_p))$layer.2
    temp$bart_prob_se <- (upper-lower)/3.92
    colnames(temp) <- c("PC1","PC2","prob","prob_se")
    temp$method <- rep("BART",nrow(temp))
    head(temp)
    
    final.PC1 <- rbind(forPC1,temp)
    head(final.PC1) 
  
#
###FOR PC2### 
#
  forPC2 <- data.frame(seq(min(data.frame(PC2)),max(data.frame(PC2)),by=0.1))
  colnames(forPC2)<- "PC2"
  forPC2$PC1 <- rep(mean(PC1$estast_ClimGridsPCA.native_1))
  head(forPC2)
  summary(forPC2)
  #GLM
  glm_p <- data.frame(predict(mod_glm,forPC2, type = "response",se=TRUE))
  forPC2$prob <- glm_p$fit
  forPC2$prob_se <- glm_p$se.fit
  forPC2$method <- rep("GLM",nrow(forPC2))
  head(forPC2)
  #BART
  x1 <- forPC2[c(1:2)]
  x1$x <- seq(1:nrow(x1))
  x1$y <- seq(1:nrow(x1))
  raster.x <- rasterFromXYZ(x1[c(3:4,1)])
  raster.y <- rasterFromXYZ(x1[c(3:4,2)])  
  bart_p <- predict(mod_bart, raster::stack(raster.x,raster.y),quantiles = c(0.2,0.8))
  temp <-  forPC2[c(1:2)]
  temp$bart_prob <- na.omit(as.data.frame(bart_p[[1]]))$layer.1
  upper <- na.omit(as.data.frame(bart_p))$layer.3
  lower <- na.omit(as.data.frame(bart_p))$layer.2
  temp$bart_prob_se <- (upper-lower)/3.92
  colnames(temp) <- c("PC2","PC1","prob","prob_se")
  temp$method <- rep("BART",nrow(temp))
  head(temp)

  final.PC2 <- rbind(forPC2,temp)
  head(final.PC2) 
  
#make graphs
maxval <- max(max(final.PC1$prob),max(final.PC2$prob))

plotPC1 <- ggplot(final.PC1,aes(x = PC1, y = prob, linetype = method)) +  geom_line(colour="black",size=0.1)
  plotPC1 <- plotPC1 + ylim(c(-0.15,maxval))
  plotPC1 <- plotPC1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plotPC1 <- plotPC1 + theme(legend.position = "none")  + ylab("suitability") 
  plotPC1 <- plotPC1 + theme(axis.title = element_text(size = 8)) + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  plotPC1 <- plotPC1 + annotate("rect", xmin = min(data.frame(climgrids_nat)[1],na.rm=TRUE), xmax = max(data.frame(climgrids_nat)[1],na.rm=TRUE), ymin = -0.06, ymax = -0.11, alpha = .3)
  plotPC1 <- plotPC1 + annotate("rect", xmin = min(data.frame(climgrids_inv)[1],na.rm=TRUE), xmax = max(data.frame(climgrids_inv)[1],na.rm=TRUE), ymin = -0.11, ymax = -0.15, alpha = .9)
  plotPC1
  
plotPC2 <- ggplot(final.PC2,aes(x = PC2, y = prob, linetype = method)) +  geom_line(colour="black",size=0.1)
  plotPC2 <- plotPC2 + ylim(c(-0.15,maxval))
  plotPC2 <- plotPC2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plotPC2 <- plotPC2 + theme(legend.position = "none") + ylab("") 
  plotPC2 <- plotPC2 + theme(axis.title = element_text(size = 6)) + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  plotPC2 <- plotPC2 + annotate("rect", xmin = min(data.frame(climgrids_nat)[2],na.rm=TRUE), xmax = max(data.frame(climgrids_nat)[2],na.rm=TRUE), ymin = -0.06, ymax = -0.11, alpha = .3)
  plotPC2 <- plotPC2 + annotate("rect", xmin = min(data.frame(climgrids_inv)[2],na.rm=TRUE), xmax = max(data.frame(climgrids_inv)[2],na.rm=TRUE), ymin = -0.11, ymax = -0.15, alpha = .9)
  plotPC2 
  
###################################################################################
###MAPS############################################################################  
###################################################################################
setwd("D:/Ststrubb/OneDrive - UGent/Projects/MC_paper/Modeling_bird_invasions/PA_modelling")
#
###READ IN DATA###
#
#mask
  mask <- raster(paste0("../../Modeling_bird_invasions/overall_eval/mask/","mask.tif"))
#predicted probabilities grid
  filenames_glm <- list.files(paste0("../../Modeling_bird_invasions/PA_modelling/pred_rasters/invasive/"),pattern="*p_glm.tif",all.files=TRUE, full.names=FALSE)
  rasters_glm <- mask + stack(paste0("../../Modeling_bird_invasions/PA_modelling/pred_rasters/invasive/",filenames_glm))
  filenames_bart <- list.files(paste0("../../Modeling_bird_invasions/PA_modelling/pred_rasters/invasive/"),pattern="*p_bart.tif",all.files=TRUE, full.names=FALSE)
  rasters_bart <- mask + stack(paste0("../../Modeling_bird_invasions/PA_modelling/pred_rasters/invasive/",filenames_bart))
 
#read in thresholds  
  tresholds5E <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/5E/DS_eval_metrics/", "mtp_thresholds.csv"))
#read in occurrence data
  occurrences <- read.csv(paste0("../../Modeling_bird_invasions/Nf_modeling/occurrences-v3/estast_invasive.csv"))
  occurrences_sparse <-  humboldt.occ.rarefy(in.pts= occurrences,colxy=2:3, rarefy.dist = 75, rarefy.units = "km")             
 
###################################################################################################################################
###RICHNESS MAPS###################################################################################################################
###################################################################################################################################
i <- 9

colors1  <- c("gray","black")
#
###GLM###
#
thresholds5E_glm <- tresholds5E[c(1:2)]
  raster_glm_5E_pa <- stack()
    m5E_glm <- c(0, thresholds5E_glm$glm_thresh[i], 0,  thresholds5E_glm$glm_thresh[i], 1, 1)
    rclmat_5E_glm <- matrix(m5E_glm, ncol=3, byrow=TRUE)  
    reclas_raster <- reclassify(rasters_glm[[i]], rclmat_5E_glm)
    raster_glm_5E_pa <- stack(raster_glm_5E_pa,reclas_raster)
    
  glm_richness_5E <-  raster_glm_5E_pa
    glm_richness_5E_points = rasterToPoints(glm_richness_5E)
    glm_richness_5E_df = data.frame(glm_richness_5E_points)
    colnames(glm_richness_5E_df) <- c("x","y","probs")
    glm_richness_5E_df$probs <- as.factor(glm_richness_5E_df$probs)
    head(glm_richness_5E_df)
  
  p_glm5E <- ggplot() + geom_raster(data=glm_richness_5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = "black")
  p_glm5E <- p_glm5E + theme(legend.position = "none")
  p_glm5E <- p_glm5E + geom_point(data=occurrences_sparse[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8)
  p_glm5E

#
###BART###
#
thresholds5E_bart <- tresholds5E[c(1,3)]
  raster_bart_5E_pa <- stack()
    m5E_bart <- c(0, thresholds5E_bart$bart_thresh[i], 0,  thresholds5E_bart$bart_thresh[i], 1, 1)
    rclmat_5E_bart <- matrix(m5E_bart, ncol=3, byrow=TRUE)  
    reclas_raster <- reclassify(rasters_bart[[i]], rclmat_5E_bart)
    raster_bart_5E_pa <- stack(raster_bart_5E_pa,reclas_raster)
  
  bart_richness_5E <-  raster_bart_5E_pa
    bart_richness_5E_points = rasterToPoints(bart_richness_5E)
    bart_richness_5E_df = data.frame(bart_richness_5E_points)
    colnames(bart_richness_5E_df) <- c("x","y","probs")
    bart_richness_5E_df$probs <- as.factor(bart_richness_5E_df$probs)
    head(bart_richness_5E_df)

  p_bart5E <- ggplot() + geom_raster(data=bart_richness_5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = "black")
    p_bart5E <- p_bart5E + theme(legend.position = "none")
    p_bart5E <- p_bart5E + geom_point(data=occurrences_sparse[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8)
    p_bart5E

#
###FNE###
#
setwd('..')
source(".\\Nf_modeling\\GEspace.R")
source(".\\Nf_modeling\\nicheG.R")
    
  # Read table with estimated parameters of niche model
  mles <- read.csv("D:/Ststrubb/OneDrive - UGent/Projects/MC_paper/Modeling_bird_invasions/Nf_modeling/Results-v3/mle_allspecies_v3.csv",header=T)
    
    
  # Project models to geographic space in the native and invaded regions
    
  k<-9
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

  
  filenames_fne <- list.files(paste0("../Modeling_bird_invasions/Nf_modeling/Results-v3/"),pattern="*suitability_map_inv.tif",all.files=TRUE, full.names=FALSE)
  rasters_fne<- mask + stack(paste0("../Modeling_bird_invasions/Nf_modeling/Results-v3/",filenames_fne))
    
thresholds5E_fne <- tresholds5E[c(1:4)]

  raster_fne_5E_pa <- stack()
    m5E_fne <- c(0, thresholds5E_fne$fne_thresh[i], 0,  thresholds5E_fne$fne_thresh[i], 1, 1)
    rclmat_5E_fne <- matrix(m5E_fne, ncol=3, byrow=TRUE)  
    reclas_raster <- reclassify(rasters_fne[[i]], rclmat_5E_fne)
    raster_fne_5E_pa <- stack(raster_fne_5E_pa,reclas_raster)
 
  fne_richness_5E <-  raster_fne_5E_pa
    fne_richness_5E_points = rasterToPoints(fne_richness_5E)
    fne_richness_5E_df = data.frame(fne_richness_5E_points)
    colnames(fne_richness_5E_df) <- c("x","y","probs")
    fne_richness_5E_df$probs <- as.factor(fne_richness_5E_df$probs)
    head(fne_richness_5E_df)
  
  p_fne5E <- ggplot() + geom_raster(data=fne_richness_5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
    p_fne5E <- p_fne5E + theme(legend.position = "none")
    p_fne5E <- p_fne5E + geom_point(data=occurrences_sparse[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8)
    p_fne5E
    
#    
###NicheMapperSpecies
#    
nm_species <- raster(paste0("../Modeling_bird_invasions/NicheMapper/","estast_yearly_energetics_pa.tif"))    
  nm_species <- nm_species+mask
  nm_species
  
  NM_points = rasterToPoints(nm_species)
    NM_points_df = data.frame(NM_points)
    colnames(NM_points_df) <- c("x","y","probs")
    NM_points_df$probs <- as.factor(NM_points_df$probs)
    head(NM_points_df)
    
  NM_fig <- ggplot() + geom_raster(data=NM_points_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
    NM_fig <- NM_fig + theme(legend.position = "none")
    NM_fig <- NM_fig + geom_point(data=occurrences_sparse[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8)
    NM_fig
    
###################################################################################################################################
###FNE ELLIPS######################################################################################################################
###################################################################################################################################    
# Read functions
source(".\\Nf_modeling\\fit_wn_maha_model.R")
    
# Packages
library(scales)
    
# Read species IDs
sp.id <- read.csv("./Nf_modeling/allspecies_samplesizes_v3.csv",header=T)[,c(1,4,6)]
    
# Summary table with estimated parameters
  mle.summary <- matrix(0,nrow = nrow(sp.id),ncol = 8)
  colnames(mle.summary) <- c("speciesID", "N", "mu1", "mu2", "sigma11",
                               "sigma12", "sigma22", "def.pos")
  # Confidence levels of ellipses in the figures
  lvs <- c(0.25,0.5,0.75,0.95)
  # Set colorpalette: 
  colpal <- c(alpha("grey70",0.7), alpha("gold3",0.7), "purple3", "grey10",
                "brown")

  j<-9
  # 1) Sample sizes
  n.nat <- nrow(occ_nat)
  n.inv <- sp.id[j,3]
  
  # 2) Read occurrence data
  # native range, used to fit the models
  f.occ.nat <- paste0("./Nf_modeling/occurrences-v3/",sp.id[j,1],"_native.csv")
  sp.occ.nat <- occ_nat[c(2:5)]
  # invasive range, used to evaluate the models
  f.occ.inv <- paste0("./Nf_modeling/occurrences-v3/",sp.id[j,1],"_invasive.csv")
  sp.occ.inv <- read.csv(f.occ.inv,header=T)[,-1]
  
  # 3) Read tables with random samples in the native and invaded areas
  # native range, used to fit the models
  f.rs.nat <- paste0("./Nf_modeling/accessible-areas-v3/",sp.id[j,1],
                     "_native_range.csv")
  sp.rs.nat <- read.csv(f.rs.nat,header=T)[,-1]
  # invasive range, used to evaluate the models
  f.rs.inv <- paste0("./Nf_modeling/accessible-areas-v3/",sp.id[j,1],
                     "_invaded_region.csv")
  sp.rs.inv <- read.csv(f.rs.inv,header=T)[,-1]
  
  # 4) Apply functions to estimate parameters
  mle <- fitNiche(E.occ = sp.occ.nat[,3:4], E.samM = sp.rs.nat[,3:4])
  
  # 5)

    # 5.1) Save estimated parameters in the summary table
    mle.summary[j,] <- c(sp.id[j,1],n.nat,mle$wn.mu,
                         as.vector(mle$wn.sigma)[-2],mle$dp)
    
    # 5.2) Define ellipses using these estimates
    ellis <- list()
    for(i in 1:length(lvs)){
      ellis[[i]] <- ellipse::ellipse(x=mle$wn.sigma, centre=as.numeric(mle$wn.mu), level=lvs[i])
    }
    
    # Visualization of results
    
    # 5.3) PLOT
    # plot will be saved as .png
    png(paste0("./Nf_modeling/",sp.id[j,1],"_modelfit.png"),
        width = 1800, height = 1800, res = 300, pointsize = 8)
    # x11()
    # background points from accessible area, using invaded regions to set plot limits
    t1<- plot(rbind(sp.rs.nat[,3:4],sp.rs.inv[,3:4]),col=colpal[1],pch=1,
         xlab="PC1", ylab="PC2",
         cex.lab=1.3,main="Environmental space")
    # add points from invaded region
    points(sp.rs.inv[,3:4],col=colpal[2],pch=1) 
    # ellipse wn
    for(x in 1:length(lvs)){lines(ellis[[x]],col=colpal[3],lwd=2)} 
    # add centre of estimated niche
    points(matrix(mle$wn.mu,ncol=2),col=colpal[3],pch=19,cex=1.5)
    # add presence points used to fit model
    points(sp.occ.nat[,3:4],col=colpal[4],pch=19,cex=1.3) 
    # add presence points use to evaluate model
    points(sp.occ.inv[,3:4],col=colpal[5],pch=17,cex=1.3)
    # figure's legend
    legend("bottomleft",legend = c(paste("Species ID:",sp.id[j,1]),"Native range",
                                   "Invaded region","Occurrences",
                                   "Recorded invasions", "Fitted model"),
           pch=c(NA,19,19,19,17,NA),col = c("white", colpal[c(1:2,4:5,3)]),
           lwd=c(rep(NA,5),2),bty = "n")
    # finish saving png
    dev.off()

#get data for plotting ellipses        
native.background <- rbind(sp.rs.nat[,3:4],sp.rs.inv[,3:4])
   native.background$area <- rep("native",nrow(native.background))
   head(native.background)
  
  invasive.background <- sp.rs.inv[,3:4]  
  invasive.background$area <- rep("invasive",nrow(invasive.background))
  head(invasive.background)
  
  background <- rbind(native.background, invasive.background)
  
  niche.center <- data.frame(matrix(mle$wn.mu,ncol=2))
    colnames(niche.center) <- c("PC1","PC2")
    niche.center$center <- "center"
  
#plot ellips
length(lvs)
  #https://htmlcolorcodes.com/

  test <- ggplot() + geom_point(data=background, aes(x=PC1, y=PC2,color=area,shape=area),size=0.5)
  test <- test + scale_color_manual(values = c("native" = "#fff333", "invasive" = "#abb2b9"))
  test <- test + scale_shape_manual(values = c("native" = 16, "invasive" = 16))
  test <- test + scale_size_manual(values = c("native" = 50, "invasive" = 16))
  test <- test  + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  test <- test + geom_point(data=sp.occ.nat[,3:4], aes(x=PC1, y=PC2),col="blue",size=0.3)
  test <- test + geom_point(data=sp.occ.inv[,3:4], aes(x=PC1, y=PC2),col='red',size=0.4)
  test <- test + theme(legend.position = "none") + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  test <- test + geom_path(data = data.frame(ellis[[1]]), aes(x = PC1 , y = PC2),size=0.02)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[2]]), aes(x = PC1 , y = PC2),size=0.02)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[3]]), aes(x = PC1 , y = PC2),size=0.02)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[4]]), aes(x = PC1 , y = PC2),size=0.02)  #add niche ellips
  test <- test + geom_point(data=niche.center, aes(x=PC1, y=PC2),shape=8,size=1)   #add niche center
  test <- test + theme(axis.title = element_text(size = 6))
  test
  
###################################################################################################################################
###NATIVE RANGE####################################################################################################################
###################################################################################################################################    
  #get background map
  nativemap <- rast("D:/Ststrubb/OneDrive - UGent/Projects/MC_paper/Modeling_bird_invasions/ranges/estast/estast_ClimGridsPCA.native.tif")
  native.map_points <-  terra::as.data.frame(nativemap*0,xy=TRUE)
  colnames(native.map_points) <- c("x","y","values")
  native.map_points$values <- as.factor(native.map_points$values)
  head(native.map_points)
  
  sp.occ.nat_sparse <- humboldt.occ.rarefy(in.pts= sp.occ.nat,colxy=2:3, rarefy.dist =30, rarefy.units = "km")  
  
native.map <- ggplot() + geom_raster(data=native.map_points,aes(x=x,y=y,fill=values))+theme_void()+scale_fill_manual(values = "#fff333")
  native.map <- native.map + theme(legend.position = "none")
  native.map <- native.map + geom_polygon(data=birdlife, aes(long, lat,group=group),col="black",fill="NA",size=0.3)
  native.map <- native.map + geom_point(data=sp.occ.nat_sparse[c(1:2)],aes(x=lon,y=lat),colour="blue",size=0.3,shape=19)
  native.map <- native.map
  native.map
  
  
####################
###GATHER FIGURES###  
####################  
plotPC1
plotPC2
p_glm5E
p_bart5E
p_fne5E
NM_fig
test
native.map
  
plotPCX <- plotPC1 + theme(legend.position = c(0.8,0.9),
                           legend.title= element_blank(),
                           legend.text=element_text(size=6),
                           legend.margin = margin(0, 0, 0, 0),
                           legend.spacing.x = unit(0, "mm"),
                           legend.spacing.y = unit(0.1, "mm"),
                           legend.key.height=unit(0.2, "cm"))
leg <- get_legend(plotPCX)


p1a <- ggarrange(p_glm5E, p_bart5E, labels = c("GLM","BART"),ncol = 2, nrow = 1,font.label = list(size = 8, color = "black"))
  p1a

p1b <- ggarrange(plotPC1, plotPC2,labels = c("",""),ncol = 2, nrow = 1,font.label = list(size = 8, color = "black"))
  p1b <- p1b 
  p1b
  
  p1 <- ggarrange(p1a,p1b,nrow=2,heights=c(1,0.75))
  p1 <- p1 + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
  p1 <- p1 + theme(plot.background = element_rect(colour = "black", fill=NA, size=0.7))
  p1
  
p2 <- ggarrange(p_fne5E, test,labels = c("FNE",""),ncol = 2, nrow = 1,font.label = list(size = 8, color = "black")) 
  p2 <- p2 + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
  p2 <- p2 + theme(plot.background = element_rect(colour = "black", fill=NA, size=0.7))
  p2 
  
p3 <- ggarrange(NM_fig,labels=c("NicheMapper"),hjust=-0.17,font.label = list(size = 8, color = "black"))
  p3 <- p3 + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
  p3 <- p3 + theme(plot.background = element_rect(colour = "black", fill=NA, size=0.7))
  p3 
  
p4 <- ggarrange(native.map,labels=c("native range"),hjust=-1.5,font.label = list(size = 8, color = "black"))
  p4 <- p4 + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
  p4 <- p4 + theme(plot.background = element_rect(colour = "black", fill=NA, size=0.))
  p4   
  
p34 <-  ggarrange(p3, p4,ncol = 2, nrow = 1) 
  
  
p.final <- ggarrange(p1,p2,p34,nrow=3,heights = c(1.5, 1,1))+ geom_subview(x=.40, y=.72, subview=leg)
  p.final <- p.final +theme(plot.background = element_rect(colour = "black", fill=NA, size=1))
  p.final <- p.final + annotate("text", label = "sensitivity: 1.00", x = 0.125, y = 0.97,size=2)
  p.final <- p.final + annotate("text", label = "sensitivity: 1.00", x = 0.625, y = 0.97,size=2)
  p.final <- p.final + annotate("text", label = "sensitivity: 1.00", x = 0.12, y = 0.54,size=2)
  p.final <- p.final + annotate("text", label = "sensitivity: 0.55", x = 0.122, y = 0.252,size=2)
  p.final
      
        
  
ggsave(filename="Estast_MultiPanel.tiff",
         plot=p.final,
         device='tiff',
         width=89,
         height=89*2,
         units="mm")       
 