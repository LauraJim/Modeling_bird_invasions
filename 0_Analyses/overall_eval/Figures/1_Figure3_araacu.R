# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(embarcadero)
library(fuzzySim)
library(modEvA)
library(maps)
library(foreach)
library(doParallel)
library(ggimage)


# CREATE OUTPUT FOLDERS ####

# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list


# GET MODEL PREDICTION RASTERS FOR EACH SPECIES ####

species <- "araacu"
  
  #pa_nat$glm_p <- predict(mod_glm, pa_nat, type = "response")
  #pred_bart <- predict_bart_df(mod_bart, pa_nat, quantiles = c(0.025, 0.975))
  #head(pred_bart)
  #pred_bart$uncert <- pred_bart[ , 3] - pred_bart[ , 2]
  #names(pred_bart) <- c("bart_p", "bart_q0025", "bart_q0975", "bart_uncert")
  #pa_nat <- data.frame(pa_nat, pred_bart)
  #head(pa_nat)
  
  #pa_nat_sv <- vect(pa_nat, geom = c("lon", "lat"))
  #par(mfrow = c(3, 1))
  #plot(pa_nat_sv, "presence", cex = 0.1, col = c(NA, "black"))
  #plot(pa_nat_sv, "glm_p", cex = 0.1)
  #plot(pa_nat_sv, "bart_p", cex = 0.1)
  
  gc()
  
  climgrids_nat <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.native.tif"))
  climgrids_nat$bias <- (climgrids_nat[[1]]*0)
  
  climgrids_inv <- rast(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/", species, "/", species, "_ClimGridsPCA.invasive.tif"))
  climgrids_inv$bias <- (climgrids_inv[[1]]*0)
  
  mod_glm <- readRDS(paste0("../../1_PA_modelling/model_objectsBIAS/", species, "_mod_glm.rds"))
  mod_bart <- readRDS(paste0("../../1_PA_modelling/model_objectsBIAS/", species, "_mod_bart.rds"))
  
  message("\n", species, "\npredicting probability on raster maps...")
  var_names <- names(mod_glm$coefficients)[-1]
  names(climgrids_nat) <- names(climgrids_inv) <- var_names
  
  message(" - GLM...")
  glm_p_nat <- predict(climgrids_nat, mod_glm, type = "response")
  glm_p_inv <- predict(climgrids_inv, mod_glm, type = "response")
  message(" - BART...")
  bart_p_nat <- predict(mod_bart, raster::stack(climgrids_nat))
  bart_p_nat <- rast(bart_p_nat)
  gc()
  bart_p_inv <- predict(mod_bart, raster::stack(climgrids_inv))
  bart_p_inv <- rast(bart_p_inv)
  
  message("converting to favourability...")
  glm_f_nat <- fuzzySim::Fav(pred = glm_p_nat, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_nat <- fuzzySim::Fav(pred = bart_p_nat, sample.preval = modEvA::prevalence(model = mod_bart))
  glm_f_inv <- fuzzySim::Fav(pred = glm_p_inv, sample.preval = modEvA::prevalence(model = mod_glm))
  bart_f_inv <- fuzzySim::Fav(pred = bart_p_inv, sample.preval = modEvA::prevalence(model = mod_bart))
  
###################################################################################
###RESPONSE CURVES#################################################################  
###################################################################################  
PC1 <- data.frame(merge(climgrids_nat[[1]],climgrids_inv[[1]]))
  PC2 <- data.frame(merge(climgrids_nat[[2]],climgrids_inv[[2]]))
  
  #
  ###FOR PC1### 
  #
  forPC1 <- data.frame(seq(min(data.frame(PC1)),max(data.frame(PC1)),by=0.1))
  colnames(forPC1)<- "PC1"
  forPC1$PC2 <- rep(mean(PC2$PC2))
  forPC1$bias <- rep(0,nrow(forPC1))
  head(forPC1)
  #GLM
  glm_p <- data.frame(predict(mod_glm,forPC1, type = "response",se=TRUE))
  forPC1$prob <- glm_p$fit
  forPC1$prob_se <- glm_p$se.fit
  forPC1$method <- rep("GLM",nrow(forPC1))
  forPC1 <- forPC1[c(1:2,4:6)]
  head(forPC1)
  #BART
  x1 <- forPC1[c(1:2)]
  x1$x <- seq(1:nrow(x1))
  x1$y <- seq(1:nrow(x1))
  raster.x <- rasterFromXYZ(x1[c(3:4,1)])
  raster.y <- rasterFromXYZ(x1[c(3:4,2)]) 
  raster.z <- raster.y*0
  preds <- raster::stack(raster.x,raster.y,raster.z)
  names(preds) <- c("PC1","PC2","bias")
  bart_p <- predict(mod_bart,preds ,quantiles = c(0.05,0.95))
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
  forPC2$PC1 <- rep(mean(PC1$PC1))
  head(forPC2)
  forPC2$bias <- rep(0,nrow(forPC2))
  summary(forPC2)
  #GLM
  glm_p <- data.frame(predict(mod_glm,forPC2, type = "response",se=TRUE))
  forPC2$prob <- glm_p$fit
  forPC2$prob_se <- glm_p$se.fit
  forPC2$method <- rep("GLM",nrow(forPC2))
  forPC2 <- forPC2[c(1:2,4:6)]
  head(forPC2)
  #BART
  x1 <- forPC2[c(1:2)]
  x1$x <- seq(1:nrow(x1))
  x1$y <- seq(1:nrow(x1))
  raster.x <- rasterFromXYZ(x1[c(3:4,1)])
  raster.y <- rasterFromXYZ(x1[c(3:4,2)]) 
  raster.z <- raster.y*0
  preds <- raster::stack(raster.x,raster.y,raster.z)
  names(preds) <- c("PC1","PC2","bias")
  bart_p <- predict(mod_bart, preds,quantiles = c(0.2,0.8))
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
  plotPC1 <- plotPC1 + ylim(c(-0.015,maxval))
  plotPC1 <- plotPC1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plotPC1 <- plotPC1 + theme(legend.position = "none")  + ylab("suitability") 
  plotPC1 <- plotPC1 + theme(axis.title = element_text(size = 8)) + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  plotPC1 <- plotPC1 + annotate("rect", xmin = min(data.frame(climgrids_nat)[1],na.rm=TRUE), xmax = max(data.frame(climgrids_nat)[1],na.rm=TRUE), ymin = -0.006, ymax = -0.011, alpha = .3)
  plotPC1 <- plotPC1 + annotate("rect", xmin = min(data.frame(climgrids_inv)[1],na.rm=TRUE), xmax = max(data.frame(climgrids_inv)[1],na.rm=TRUE), ymin = -0.011, ymax = -0.015, alpha = .9)
  plotPC1
  
  plotPC2 <- ggplot(final.PC2,aes(x = PC2, y = prob, linetype = method)) +  geom_line(colour="black",size=0.1)
  plotPC2 <- plotPC2 + ylim(c(-0.015,maxval))
  plotPC2 <- plotPC2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  plotPC2 <- plotPC2 + theme(legend.position = "none") + ylab("") 
  plotPC2 <- plotPC2 + theme(axis.title = element_text(size = 6)) + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  plotPC2 <- plotPC2 + annotate("rect", xmin = min(data.frame(climgrids_nat)[2],na.rm=TRUE), xmax = max(data.frame(climgrids_nat)[2],na.rm=TRUE), ymin = -0.006, ymax = -0.011, alpha = .3)
  plotPC2 <- plotPC2 + annotate("rect", xmin = min(data.frame(climgrids_inv)[2],na.rm=TRUE), xmax = max(data.frame(climgrids_inv)[2],na.rm=TRUE), ymin = -0.011, ymax = -0.015, alpha = .9)
  plotPC2 
  
###################################################################################
###MAPS############################################################################  
###################################################################################
#
###READ IN DATA###
#

#mask
  mask <- raster(paste0("../../../0_Analyses/overall_eval/Figures/","mask.tif"))

  #read in thresholds
  threshold_2.5E <- read.csv(paste0("../../1_PA_modelling/eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics/crossval_metrics_nat_araacu.csv"))
  
  #read in occurrence data
  occurrences <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/araacu_invasive.csv"))
  
  ###################################################################################################################################
  ###RICHNESS MAPS###################################################################################################################
  ###################################################################################################################################
  i <- 7
  
  colors1  <- c("gray","black")
  #
  ###GLM###
  #
  raster_glm_2.5E_pa <- stack()
  m2.5E_glm <- c(0, mean(threshold_2.5E$threshglm), 0,  mean(threshold_2.5E$threshglm), 1, 1)
  rclmat_2.5E_glm <- matrix(m2.5E_glm, ncol=3, byrow=TRUE)  
  reclas_raster <- reclassify(raster(glm_p_inv), rclmat_2.5E_glm)
  raster_glm_2.5E_pa <- stack(raster_glm_2.5E_pa,reclas_raster)
  
  glm_richness_2.5E <-  raster_glm_2.5E_pa
  glm_richness_2.5E_points = rasterToPoints(glm_richness_2.5E)
  glm_richness_2.5E_df = data.frame(glm_richness_2.5E_points)
  colnames(glm_richness_2.5E_df) <- c("x","y","probs")
  glm_richness_2.5E_df$probs <- as.factor(glm_richness_2.5E_df$probs)
  head(glm_richness_2.5E_df)
  
  p_glm5E <- ggplot() + geom_raster(data=glm_richness_2.5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
  p_glm5E <- p_glm5E + theme(legend.position = "none")
  p_glm5E <- p_glm5E + geom_point(data=occurrences[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8,stroke=0)
  p_glm5E
  
  #
  ###BART###
  #
  raster_bart_2.5E_pa <- stack()
  m2.5E_bart <- c(0, mean(threshold_2.5E$threshbart), 0,  mean(threshold_2.5E$threshbart), 1, 1)
  rclmat_2.5E_bart <- matrix(m2.5E_bart, ncol=3, byrow=TRUE)  
  reclas_raster <- reclassify(raster(bart_p_inv), rclmat_2.5E_bart)
  raster_bart_2.5E_pa <- stack(raster_bart_2.5E_pa,reclas_raster)
  
  bart_richness_2.5E <-  raster_bart_2.5E_pa
  bart_richness_2.5E_points = rasterToPoints(bart_richness_2.5E)
  bart_richness_2.5E_df = data.frame(bart_richness_2.5E_points)
  colnames(bart_richness_2.5E_df) <- c("x","y","probs")
  bart_richness_2.5E_df$probs <- as.factor(bart_richness_2.5E_df$probs)
  head(bart_richness_2.5E_df)
  
  p_bart5E <- ggplot() + geom_raster(data=bart_richness_2.5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
  p_bart5E <- p_bart5E + theme(legend.position = "none")
  p_bart5E <- p_bart5E + geom_point(data=occurrences[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8,stroke=0)
  p_bart5E
  
  #
  ###FNE###
  #
  source("../../0_Nf_modeling/GEspace.R")
  source("../../0_Nf_modeling/nicheG.R")
  
  # Read table with estimated parameters of niche model
  mles <- read.csv("../../0_Nf_modeling/ResultsClimate/mle_allspecies_v3.csv",header=T)
  
  
  # Project models to geographic space in the native and invaded regions
  
  k<-7
  # Select species
  species.id <- mles[k,1]
  # Project model into the native range
  pc.native <- stack(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/",paste0(species.id,"/"),
                            species.id,"_ClimGridsPCA.native.tif"))
  suit.wn.nat <- niche.G(Estck = pc.native, mu = c(mles[k,3],mles[k,4]), 
                         Sigma = matrix(c(mles[k,5], mles[k,6], mles[k,6], 
                                          mles[k,7]),ncol=2))
  
  # Project model into the invaded region
  pc.invasive <- stack(paste0("../../../1_Data/sdm_predictors/bioreg_occ_dist_climgrids/",paste0(species.id,"/"),
                              species.id,"_ClimGridsPCA.invasive.tif"))
  
  suit.wn.inv <- niche.G(Estck = pc.invasive, mu = c(mles[k,3],mles[k,4]), 
                         Sigma = matrix(c(mles[k,5], mles[k,6], mles[k,6], 
                                          mles[k,7]),ncol=2))
 
  filenames_fne <- list.files(paste0("../../0_Nf_modeling/ResultsClimate/"),pattern="*suitability_map_inv.tif",all.files=TRUE, full.names=FALSE)
  rasters_fne<- mask + stack(paste0("../../0_Nf_modeling/ResultsClimate/",filenames_fne))
  
  #get threshold
  native_range <- raster(paste0("../../../0_Analyses/0_Nf_modeling/ResultsClimate/",species,"_suitability_map_nat.tif"))
  native_range_occs <- read.csv(paste0("../../../0_Analyses/0_Nf_modeling/occurrences-v3/",species,"_native.csv"))
  native_range_suitabilities <- raster::extract(native_range,native_range_occs[c(2:3)])
  
  pred_nat <- read.csv(paste0("../../../0_Analyses/1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/pred_CSVsBIAS/", species, "_pred_crossval.csv"))[c(2:3,7)]
  pred_nat$fne_suit <- raster::extract(native_range,pred_nat[c(1:2)]) 
  thresh <- quantile(native_range_suitabilities, probs = 0.025, na.rm = TRUE)  # pred for minimum 5% presence in the training data
  
  
  raster_fne_2.5E_pa <- stack()
  m2.5E_fne <- c(0, thresh, 0,  thresh, 1, 1)
  rclmat_2.5E_fne <- matrix(m2.5E_fne, ncol=3, byrow=TRUE)  
  reclas_raster <- reclassify(suit.wn.inv, rclmat_2.5E_fne)
  raster_fne_2.5E_pa <- stack(raster_fne_2.5E_pa,reclas_raster)
  
  fne_richness_2.5E <-  raster_fne_2.5E_pa
  fne_richness_2.5E_points = rasterToPoints(fne_richness_2.5E)
  fne_richness_2.5E_df = data.frame(fne_richness_2.5E_points)
  colnames(fne_richness_2.5E_df) <- c("x","y","probs")
  fne_richness_2.5E_df$probs <- as.factor(fne_richness_2.5E_df$probs)
  head(fne_richness_2.5E_df)
  
  p_fne5E <- ggplot() + geom_raster(data=fne_richness_2.5E_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
  p_fne5E <- p_fne5E + theme(legend.position = "none")
  p_fne5E <- p_fne5E + geom_point(data=occurrences[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8,stroke=0)
  p_fne5E
  
  #    
  ###NicheMapperSpecies
  #    
  nm_species <- raster(paste0("../../4_NicheMapper/","ARAACU_yearly_energetics_pa.tif"))    
  nm_species <- nm_species+mask
  nm_species
  
  NM_points = rasterToPoints(nm_species)
  NM_points_df = data.frame(NM_points)
  colnames(NM_points_df) <- c("x","y","probs")
  NM_points_df$probs <- as.factor(NM_points_df$probs)
  head(NM_points_df)
  
  NM_fig <- ggplot() + geom_raster(data=NM_points_df,aes(x=x,y=y,fill=probs))+theme_void()+scale_fill_manual(values = colors1)
  NM_fig <- NM_fig + theme(legend.position = "none")
  NM_fig <- NM_fig + geom_point(data=occurrences[c(2:3)],aes(x=lon,y=lat),colour="red",size=0.8,stroke=0)
  NM_fig
  
  ###################################################################################################################################
  ###FNE ELLIPS######################################################################################################################
  ###################################################################################################################################    
  # Read functions
  source("../../0_Nf_modeling/fit_wn_maha_model.R")
  
  # Packages
  library(scales)
  
  # Read species IDs
  sp.id <- read.csv("../../0_Nf_modeling/allspecies_samplesizes_v3.csv",header=T)[,c(1,4,6)]
  
  # Summary table with estimated parameters
  mle.summary <- matrix(0,nrow = nrow(sp.id),ncol = 8)
  colnames(mle.summary) <- c("speciesID", "N", "mu1", "mu2", "sigma11",
                             "sigma12", "sigma22", "def.pos")
  # Confidence levels of ellipses in the figures
  lvs <- c(0.25,0.5,0.75,0.95)
  # Set colorpalette: 
  colpal <- c(alpha("grey70",0.7), alpha("gold3",0.7), "purple3", "grey10",
              "brown")
              
  j<-7  
  # 1) Sample sizes
  n.nat <- sp.id[j,2]
  n.inv <- sp.id[j,3]
  
  # 2) Read occurrence data
  # native range, used to fit the models
  sp.occ.nat <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/",sp.id[j,1],"_native.csv"))
  sp.occ.nat <- sp.occ.nat[c(2:5)]
  
  # invasive range, used to evaluate the models
  sp.occ.inv <- read.csv(paste0("../../0_Nf_modeling/occurrences-v3/",sp.id[j,1],"_invasive.csv"))
  sp.occ.inv <- sp.occ.inv[c(2:5)]
  
  # 3) Read tables with random samples in the native and invaded areas
  # native range, used to fit the models
  f.rs.nat <- paste0("../../0_Nf_modeling/accessible-areas-v3/",sp.id[j,1],"_native_range.csv")
  sp.rs.nat <- read.csv(f.rs.nat,header=T)[,-1]
  # invasive range, used to evaluate the models
  f.rs.inv <- paste0("../../0_Nf_modeling/accessible-areas-v3/",sp.id[j,1],
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
  png(paste0("../../0_Nf_modeling/",sp.id[j,1],"_modelfit.png"),
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
  test <- test + geom_point(data=sp.occ.nat[,3:4], aes(x=PC1, y=PC2),col="blue",size=0.5,stroke=0,shape=19)
  test <- test + geom_point(data=sp.occ.inv[,3:4], aes(x=PC1, y=PC2),col='red',size=0.7,stroke=0,shape=19)
  test <- test + theme(legend.position = "none") + theme(axis.text.x = element_text(size=6)) + theme(axis.text.y = element_text(size=6))
  test <- test + geom_path(data = data.frame(ellis[[1]]), aes(x = PC1 , y = PC2),size=0.02,linewidth=0.25)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[2]]), aes(x = PC1 , y = PC2),size=0.02,linewidth=0.25)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[3]]), aes(x = PC1 , y = PC2),size=0.02,linewidth=0.25)  #add niche ellips
  test <- test + geom_path(data = data.frame(ellis[[4]]), aes(x = PC1 , y = PC2),size=0.02,linewidth=0.25)  #add niche ellips
  test <- test + geom_point(data=niche.center, aes(x=PC1, y=PC2),shape=8,size=1)   #add niche center
  test <- test + theme(axis.title = element_text(size = 6))
  test
  
  ###################################################################################################################################
  ###NATIVE RANGE####################################################################################################################
  ###################################################################################################################################    
  #get background map
  nativemap <- raster(paste0("../../../1_Data/ranges/araacu/nativemap.tif"))
  
  native.map_points <-  terra::as.data.frame(nativemap*0,xy=TRUE)
  colnames(native.map_points) <- c("x","y","values")
  native.map_points$values <- as.factor(native.map_points$values)
  head(native.map_points)
  
  #get native range BirdLife   
  birdlife.raw<-readOGR(paste0("../../../1_Data/ranges/araacu"),"araacu")
  birdlife.raw <- birdlife.raw[birdlife.raw$ORIGIN==1,]
  birdlife <- fortify(birdlife.raw)
  
  native.map <- ggplot() + geom_raster(data=native.map_points,aes(x=x,y=y,fill=values))+theme_void()+scale_fill_manual(values = "#fff333",na.value="white")
  native.map <- native.map + theme(legend.position = "none")
  native.map <- native.map + geom_polygon(data=birdlife, aes(long, lat,group=group),col="black",fill="NA",size=0.3)
  native.map <- native.map + geom_point(data=sp.occ.nat[c(1:2)],aes(x=lon,y=lat),colour="blue",size=0.5,stroke=0,shape=19)
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
  p4 <- p4 + theme(plot.background = element_rect(colour = "black", fill=NA, size=0.7))
  p4   
  
  p34 <-  ggarrange(p3, p4,ncol = 2, nrow = 1) 
  
  
  p.final <- ggarrange(p1,p2,p34,nrow=3,heights = c(1.5, 1,1))+ geom_subview(x=.40, y=.72, subview=leg)
  p.final <- p.final +theme(plot.background = element_rect(colour = "black", fill=NA, size=1))
  p.final <- p.final + annotate("text", label = "sensitivity: 1.00", x = 0.125, y = 0.97,size=2)
  p.final <- p.final + annotate("text", label = "specificity: 0.30", x = 0.125, y = 0.955,size=2)
  
  p.final <- p.final + annotate("text", label = "sensitivity: 0.91", x = 0.625, y = 0.97,size=2)
  p.final <- p.final + annotate("text", label = "specificity: 0.84", x = 0.625, y = 0.955,size=2)
  
  
  p.final <- p.final + annotate("text", label = "sensitivity: 0.91", x = 0.12, y = 0.54,size=2)
  p.final <- p.final + annotate("text", label = "specificity: 0.85", x = 0.12, y = 0.525,size=2)
  
  
  p.final <- p.final + annotate("text", label = "sensitivity: 0.87", x = 0.122, y = 0.252,size=2)
  p.final <- p.final + annotate("text", label = "specificity: 0.92", x = 0.122, y = 0.237,size=2)
  
  
  p.final
  
  
  
  ggsave(filename="Araacu_MultiPanel.tiff",
         plot=p.final,
         device='tiff',
         width=89,
         height=89*2,
         units="mm")       
  
  