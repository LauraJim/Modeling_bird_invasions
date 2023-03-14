# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# LOAD PACKAGES ####
library(raster)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(terra)
###################################################################################################################################
###READ IN DATA####################################################################################################################
###################################################################################################################################
#mask
mask <- raster("mask.tif")


#predicted probabilities grid
setwd("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/pred_rastersBIAS_CLAMPING/EU/invasive")
  glm_BIAS_CLAMPING_files <- list.files(pattern="*p_glm.tif",all.files=TRUE, full.names=FALSE)
  glm_BIAS_CLAMPING <- stack(glm_BIAS_CLAMPING_files)
  names(glm_BIAS_CLAMPING) <- glm_BIAS_CLAMPING_files
  
setwd("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/pred_rastersBIAS/EU/invasive")
  glm_BIAS_EXTRAPOL_files <- list.files(pattern="*p_glm.tif",all.files=TRUE, full.names=FALSE)
  glm_BIAS_EXTRAPOL <- stack(glm_BIAS_EXTRAPOL_files)
  names(glm_BIAS_EXTRAPOL) <- glm_BIAS_EXTRAPOL_files  

setwd("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/pred_rastersBIAS_CLAMPING/EU/invasive")
  bart_BIAS_CLAMPING_files <- list.files(pattern="*p_bart.tif",all.files=TRUE, full.names=FALSE)
  bart_BIAS_CLAMPING <- stack(bart_BIAS_CLAMPING_files)
  names(bart_BIAS_CLAMPING) <- bart_BIAS_CLAMPING_files
  
setwd("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/pred_rastersBIAS/EU/invasive")
  bart_BIAS_EXTRAPOL_files <- list.files(pattern="*p_bart.tif",all.files=TRUE, full.names=FALSE)
  bart_BIAS_EXTRAPOL <- stack(bart_BIAS_EXTRAPOL_files)
  names(bart_BIAS_EXTRAPOL) <- bart_BIAS_EXTRAPOL_files   
  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd() 
  
#read in thresholds  
tresholds5E <- read.csv("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/5E/mtp_thresholds.csv")
tresholds2.5E <- read.csv("E:/Modeling_bird_invasions/0_Analyses/1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/2.5E/mtp_thresholds.csv")


###################################################################################################################################
###RICHNESS MAPS###################################################################################################################
###################################################################################################################################
#
###GLM#############################################################################################################################
#
thresholds5E_glm <- tresholds5E[c(1:2)]
thresholds2.5E_glm <- tresholds2.5E[c(1:2)]

#5E CLAMPING
dir.create("BIAS/GLM_BIAS_CLAMPING_5E")
raster_glm_5E_pa <- stack()
   for (i in 1:20){
      m5E_glm <- c(0, thresholds5E_glm$glm_thresh[i], 0,  thresholds5E_glm$glm_thresh[i], 1, 1)
         rclmat_5E_glm <- matrix(m5E_glm, ncol=3, byrow=TRUE)  
         reclas_raster <- reclassify(glm_BIAS_CLAMPING[[i]], rclmat_5E_glm)
         jpeg(paste0("BIAS/GLM_BIAS_CLAMPING_5E/", names(glm_BIAS_CLAMPING[[i]]), ".jpg"), width = 430, height = 430)
         plot(trim(rast(reclas_raster)))
         dev.off()
         raster_glm_5E_pa <- stack(raster_glm_5E_pa,reclas_raster)
   }

  glm_richness_5E <-  sum(raster_glm_5E_pa) 
   glm_richness_5E_points = rasterToPoints(glm_richness_5E)
   glm_richness_5E_df = data.frame(glm_richness_5E_points)
   head(glm_richness_5E_df) #breaks will be set to column "layer"
   glm_richness_5E_df$layer[1] <- 20
   glm_richness_5E_df$layer[2] <- 1
   
   p_glm5E <- ggplot() + geom_raster(data=glm_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))
   p_glm5E <- p_glm5E+theme(legend.title=element_blank())
   p_glm5E
   
   ggsave(filename="GLM_BIAS_CLAMPING_5E.tiff",
          plot=p_glm5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
#2.5E CLAMPING
   dir.create("BIAS/GLM_BIAS_CLAMPING_2.5E")
   raster_glm_2.5E_pa <- stack()
   for (i in 1:20){
     m2.5E_glm <- c(0, thresholds2.5E_glm$glm_thresh[i], 0,  thresholds2.5E_glm$glm_thresh[i], 1, 1)
     rclmat_2.5E_glm <- matrix(m2.5E_glm, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(glm_BIAS_CLAMPING[[i]], rclmat_2.5E_glm)
     jpeg(paste0("BIAS/GLM_BIAS_CLAMPING_2.5E/", names(glm_BIAS_CLAMPING[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_glm_2.5E_pa <- stack(raster_glm_2.5E_pa,reclas_raster)
   }
   
   glm_richness_2.5E <-  sum(raster_glm_2.5E_pa) 
   glm_richness_2.5E_points = rasterToPoints(glm_richness_2.5E)
   glm_richness_2.5E_df = data.frame(glm_richness_2.5E_points)
   head(glm_richness_2.5E_df) #breaks will be set to column "layer"
   glm_richness_2.5E_df$layer[1] <- 20
   glm_richness_2.5E_df$layer[2] <- 1
   
   p_glm2.5E <- ggplot() + geom_raster(data=glm_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_glm2.5E <- p_glm2.5E+theme(legend.title=element_blank())
   p_glm2.5E
   
   ggsave(filename="GLM_BIAS_CLAMPING_2.5E.tiff",
          plot=p_glm5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")   
   
#5E EXTRAPOL
dir.create("BIAS/GLM_BIAS_EXTRAPOL_5E")
 raster_glm_5E_pa <- stack()
   for (i in 1:20){
     m5E_glm <- c(0, thresholds5E_glm$glm_thresh[i], 0,  thresholds5E_glm$glm_thresh[i], 1, 1)
     rclmat_5E_glm <- matrix(m5E_glm, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(glm_BIAS_EXTRAPOL[[i]], rclmat_5E_glm)
     jpeg(paste0("BIAS/GLM_BIAS_EXTRAPOL_5E/", names(glm_BIAS_EXTRAPOL[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_glm_5E_pa <- stack(raster_glm_5E_pa,reclas_raster)
   }
   
   glm_richness_5E <-  sum(raster_glm_5E_pa) 
   glm_richness_5E_points = rasterToPoints(glm_richness_5E)
   glm_richness_5E_df = data.frame(glm_richness_5E_points)
   head(glm_richness_5E_df) #breaks will be set to column "layer"
   glm_richness_5E_df$layer[1] <- 20
   glm_richness_5E_df$layer[2] <- 1
   
   p_glm5E <- ggplot() + geom_raster(data=glm_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_glm5E <- p_glm5E+theme(legend.title=element_blank())
   p_glm5E
   
   ggsave(filename="GLM_BIAS_EXTRAPOL_5E.tiff",
          plot=p_glm5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
#2.5E EXTRAPOL
   dir.create("BIAS/GLM_BIAS_EXTRAPOL_2.5E")
   raster_glm_2.5E_pa <- stack()
   for (i in 1:20){
     m2.5E_glm <- c(0, thresholds2.5E_glm$glm_thresh[i], 0,  thresholds2.5E_glm$glm_thresh[i], 1, 1)
     rclmat_2.5E_glm <- matrix(m2.5E_glm, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(glm_BIAS_EXTRAPOL[[i]], rclmat_2.5E_glm)
     jpeg(paste0("BIAS/GLM_BIAS_EXTRAPOL_2.5E/", names(glm_BIAS_EXTRAPOL[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_glm_2.5E_pa <- stack(raster_glm_2.5E_pa,reclas_raster)
   }
   
   glm_richness_2.5E <-  sum(raster_glm_2.5E_pa) 
   glm_richness_2.5E_points = rasterToPoints(glm_richness_2.5E)
   glm_richness_2.5E_df = data.frame(glm_richness_2.5E_points)
   head(glm_richness_2.5E_df) #breaks will be set to column "layer"
   glm_richness_2.5E_df$layer[1] <- 20
   glm_richness_2.5E_df$layer[2] <- 1
   
   p_glm2.5E <- ggplot() + geom_raster(data=glm_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_glm2.5E <- p_glm2.5E+theme(legend.title=element_blank())
   p_glm2.5E
   
   ggsave(filename="GLM_BIAS_EXTRAPOL_2.5E.tiff",
          plot=p_glm5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")      
   
#
###BART#############################################################################################################################
#
thresholds5E_bart <- tresholds5E[c(1,3)]
thresholds2.5E_bart <- tresholds2.5E[c(1,3)]
   
   #5E CLAMPING
   dir.create("BIAS/BART_BIAS_CLAMPING_5E")
   raster_bart_5E_pa <- stack()
   for (i in 1:20){
     m5E_bart <- c(0, thresholds5E_bart$bart_thresh[i], 0,  thresholds5E_bart$bart_thresh[i], 1, 1)
     rclmat_5E_bart <- matrix(m5E_bart, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(bart_BIAS_CLAMPING[[i]], rclmat_5E_bart)
     jpeg(paste0("BIAS./BART_BIAS_CLAMPING_5E/", names(bart_BIAS_CLAMPING[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_bart_5E_pa <- stack(raster_bart_5E_pa,reclas_raster)
   }
   
   bart_richness_5E <-  sum(raster_bart_5E_pa) 
   bart_richness_5E_points = rasterToPoints(bart_richness_5E)
   bart_richness_5E_df = data.frame(bart_richness_5E_points)
   head(bart_richness_5E_df) #breaks will be set to column "layer"
   bart_richness_5E_df$layer[1] <- 20
   bart_richness_5E_df$layer[2] <- 1
   
   p_bart5E <- ggplot() + geom_raster(data=bart_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_bart5E <- p_bart5E+theme(legend.title=element_blank())
   p_bart5E
   
   ggsave(filename="BART_BIAS_CLAMPING_5E.tiff",
          plot=p_bart5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
   #2.5E CLAMPING
   dir.create("BIAS/BART_BIAS_CLAMPING_2.5E")
   raster_bart_2.5E_pa <- stack()
   for (i in 1:20){
     m2.5E_bart <- c(0, thresholds2.5E_bart$bart_thresh[i], 0,  thresholds2.5E_bart$bart_thresh[i], 1, 1)
     rclmat_2.5E_bart <- matrix(m2.5E_bart, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(bart_BIAS_CLAMPING[[i]], rclmat_2.5E_bart)
     jpeg(paste0("BIAS./BART_BIAS_CLAMPING_2.5E/", names(bart_BIAS_CLAMPING[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_bart_2.5E_pa <- stack(raster_bart_2.5E_pa,reclas_raster)
   }
   
   bart_richness_2.5E <-  sum(raster_bart_2.5E_pa) 
   bart_richness_2.5E_points = rasterToPoints(bart_richness_2.5E)
   bart_richness_2.5E_df = data.frame(bart_richness_2.5E_points)
   head(bart_richness_2.5E_df) #breaks will be set to column "layer"
   bart_richness_2.5E_df$layer[1] <- 20
   bart_richness_2.5E_df$layer[2] <- 1
   
   p_bart2.5E <- ggplot() + geom_raster(data=bart_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_bart2.5E <- p_bart2.5E+theme(legend.title=element_blank())
   p_bart2.5E
   
   ggsave(filename="BART_BIAS_CLAMPING_2.5E.tiff",
          plot=p_bart5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")   
   
   #5E EXTRAPOL
   dir.create("BIAS/BART_BIAS_EXTRAPOL_5E")
   raster_bart_5E_pa <- stack()
   for (i in 1:20){
     m5E_bart <- c(0, thresholds5E_bart$bart_thresh[i], 0,  thresholds5E_bart$bart_thresh[i], 1, 1)
     rclmat_5E_bart <- matrix(m5E_bart, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(bart_BIAS_EXTRAPOL[[i]], rclmat_5E_bart)
     jpeg(paste0("BIAS./BART_BIAS_EXTRAPOL_5E/", names(bart_BIAS_EXTRAPOL[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_bart_5E_pa <- stack(raster_bart_5E_pa,reclas_raster)
   }
   
   bart_richness_5E <-  sum(raster_bart_5E_pa) 
   bart_richness_5E_points = rasterToPoints(bart_richness_5E)
   bart_richness_5E_df = data.frame(bart_richness_5E_points)
   head(bart_richness_5E_df) #breaks will be set to column "layer"
   bart_richness_5E_df$layer[1] <- 20
   bart_richness_5E_df$layer[2] <- 1
   
   p_bart5E <- ggplot() + geom_raster(data=bart_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_bart5E <- p_bart5E+theme(legend.title=element_blank())
   p_bart5E
   
   ggsave(filename="BART_BIAS_EXTRAPOL_5E.tiff",
          plot=p_bart5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
   #2.5E EXTRAPOL
   dir.create("BIAS/BART_BIAS_EXTRAPOL_2.5E")
   raster_bart_2.5E_pa <- stack()
   for (i in 1:20){
     m2.5E_bart <- c(0, thresholds2.5E_bart$bart_thresh[i], 0,  thresholds2.5E_bart$bart_thresh[i], 1, 1)
     rclmat_2.5E_bart <- matrix(m2.5E_bart, ncol=3, byrow=TRUE)  
     reclas_raster <- reclassify(bart_BIAS_EXTRAPOL[[i]], rclmat_2.5E_bart)
     jpeg(paste0("BIAS./BART_BIAS_EXTRAPOL_2.5E/", names(bart_BIAS_EXTRAPOL[[i]]), ".jpg"), width = 430, height = 430)
     plot(trim(rast(reclas_raster)))
     dev.off()
     raster_bart_2.5E_pa <- stack(raster_bart_2.5E_pa,reclas_raster)
   }
   
   bart_richness_2.5E <-  sum(raster_bart_2.5E_pa) 
   bart_richness_2.5E_points = rasterToPoints(bart_richness_2.5E)
   bart_richness_2.5E_df = data.frame(bart_richness_2.5E_points)
   head(bart_richness_2.5E_df) #breaks will be set to column "layer"
   bart_richness_2.5E_df$layer[1] <- 20
   bart_richness_2.5E_df$layer[2] <- 1
   
   p_bart2.5E <- ggplot() + geom_raster(data=bart_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
     scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
     guides(fill=guide_legend(title="species richness"))
   p_bart2.5E <- p_bart2.5E+theme(legend.title=element_blank())
   p_bart2.5E
   
   ggsave(filename="BART_BIAS_EXTRAPOL_2.5E.tiff",
          plot=p_bart5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")      
   
   
   
###############################################################################################################################################   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
#2.5E
raster_glm_2.5E_pa <- stack()
   for (i in 1:20){
   m2.5E_glm <- c(0, thresholds2.5E_glm$glm_thresh[i], 0,  thresholds2.5E_glm$glm_thresh[i], 1, 1)
      rclmat_2.5E_glm <- matrix(m2.5E_glm, ncol=3, byrow=TRUE)  
      reclas_raster <- reclassify(rasters_glm[[i]], rclmat_2.5E_glm)
      raster_glm_2.5E_pa <- stack(raster_glm_2.5E_pa,reclas_raster)
   }
   glm_richness_2.5E <-  sum(raster_glm_2.5E_pa)      
   
   glm_richness_2.5E_points = rasterToPoints(glm_richness_2.5E)
   glm_richness_2.5E_df = data.frame(glm_richness_2.5E_points)
   head(glm_richness_2.5E_df) #breaks will be set to column "layer"
   glm_richness_2.5E_df$layer[1] <- 20
   glm_richness_2.5E_df$layer[2] <- 1
   
   p_glm2.5E <- ggplot() + geom_raster(data=glm_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))+theme(legend.title=element_blank())
   p_glm2.5E
   
   ggsave(filename="Fig_Revision_Xb2_sp_rich(glm2.5E).tiff",
          plot=p_glm2.5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")       

#
###BART#############################################################################################################################
#
thresholds5E_bart <- tresholds5E[c(1,3)]
thresholds2.5E_bart <- tresholds2.5E[c(1,3)]

#5E
raster_bart_5E_pa <- stack()
   for (i in 1:20){
      m5E_bart <- c(0, thresholds5E_bart$bart_thresh[i], 0,  thresholds5E_bart$bart_thresh[i], 1, 1)
      rclmat_5E_bart <- matrix(m5E_bart, ncol=3, byrow=TRUE)  
      reclas_raster <- reclassify(rasters_bart[[i]], rclmat_5E_bart)
      raster_bart_5E_pa <- stack(raster_bart_5E_pa,reclas_raster)
   }
   bart_richness_5E <-  sum(raster_bart_5E_pa)      
   
   bart_richness_5E_points = rasterToPoints(bart_richness_5E)
   bart_richness_5E_df = data.frame(bart_richness_5E_points)
   head(bart_richness_5E_df) #breaks will be set to column "layer"
   bart_richness_5E_df$layer[1] <- 20
   bart_richness_5E_df$layer[2] <- 1
   
   p_bart5E <- ggplot() + geom_raster(data=bart_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))+theme(legend.title=element_blank())
   p_bart5E
   
   ggsave(filename="Fig_Revision_Xc1_sp_rich(bart5E).tiff",
          plot=p_bart5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
#2.5E
raster_bart_2.5E_pa <- stack()
   for (i in 1:20){
      m2.5E_bart <- c(0, thresholds2.5E_bart$bart_thresh[i], 0,  thresholds2.5E_bart$bart_thresh[i], 1, 1)
      rclmat_2.5E_bart <- matrix(m2.5E_bart, ncol=3, byrow=TRUE)  
      reclas_raster <- reclassify(rasters_bart[[i]], rclmat_2.5E_bart)
      raster_bart_2.5E_pa <- stack(raster_bart_2.5E_pa,reclas_raster)
   }
   bart_richness_2.5E <-  sum(raster_bart_2.5E_pa)      
   
   bart_richness_2.5E_points = rasterToPoints(bart_richness_2.5E)
   bart_richness_2.5E_df = data.frame(bart_richness_2.5E_points)
   head(bart_richness_2.5E_df) #breaks will be set to column "layer"
   bart_richness_2.5E_df$layer[1] <- 20
   bart_richness_2.5E_df$layer[2] <- 1
   
   p_bart2.5E <- ggplot() + geom_raster(data=bart_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))+theme(legend.title=element_blank())
   p_bart2.5E
   
   ggsave(filename="Fig_Revision_Xc2_sp_rich(bart2.5E).tiff",
          plot=p_bart2.5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")       
   
#
###FNE#############################################################################################################################
#
thresholds5E_fne <- tresholds5E[c(1:4)]
thresholds2.5E_fne <- tresholds2.5E[c(1:4)]
   
#5E
raster_fne_5E_pa <- stack()
   for (i in 1:20){
      m5E_fne <- c(0, thresholds5E_fne$fne_thresh[i], 0,  thresholds5E_fne$fne_thresh[i], 1, 1)
      rclmat_5E_fne <- matrix(m5E_fne, ncol=3, byrow=TRUE)  
      reclas_raster <- reclassify(rasters_fne[[i]], rclmat_5E_fne)
      raster_fne_5E_pa <- stack(raster_fne_5E_pa,reclas_raster)
   }
   fne_richness_5E <-  sum(raster_fne_5E_pa)      
   
   fne_richness_5E_points = rasterToPoints(fne_richness_5E)
   fne_richness_5E_df = data.frame(fne_richness_5E_points)
   head(fne_richness_5E_df) #breaks will be set to column "layer"
   fne_richness_5E_df$layer[1] <- 20
   fne_richness_5E_df$layer[2] <- 1
   
   p_fne5E <- ggplot() + geom_raster(data=fne_richness_5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))+theme(legend.title=element_blank())
   p_fne5E
   
   ggsave(filename="Fig_Revision_Xd1_sp_rich(fne5E).tiff",
          plot=p_fne5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")
   
#2.5E
raster_fne_2.5E_pa <- stack()
   for (i in 1:20){
      m2.5E_fne <- c(0, thresholds2.5E_fne$fne_thresh[i], 0,  thresholds2.5E_fne$fne_thresh[i], 1, 1)
      rclmat_2.5E_fne <- matrix(m2.5E_fne, ncol=3, byrow=TRUE)  
      reclas_raster <- reclassify(rasters_fne[[i]], rclmat_2.5E_fne)
      raster_fne_2.5E_pa <- stack(raster_fne_2.5E_pa,reclas_raster)
   }
   fne_richness_2.5E <-  sum(raster_fne_2.5E_pa)      
   
   fne_richness_2.5E_points = rasterToPoints(fne_richness_2.5E)
   fne_richness_2.5E_df = data.frame(fne_richness_2.5E_points)
   head(fne_richness_2.5E_df) #breaks will be set to column "layer"
   fne_richness_2.5E_df$layer[1] <- 20
   fne_richness_2.5E_df$layer[2] <- 1
   
   p_fne2.5E <- ggplot() + geom_raster(data=fne_richness_2.5E_df,aes(x=x,y=y,fill=layer))+theme_void()+
      scale_fill_viridis_c(option = "D",breaks=seq(1:20))+
      guides(fill=guide_legend(title="species richness"))
   p_fne2.5E <- p_fne2.5E+theme(legend.title=element_blank())
   p_fne2.5E
   
   ggsave(filename="Fig_Revision_Xd2_sp_rich(fne2.5E).tiff",
          plot=p_fne2.5E,
          device='tiff',
          width=210,
          height=297,
          units="mm")       

############################################################      
###COMBINE IN A SINGLE PLOT AND STANDARDIZE RANGE 0 to 20###
############################################################
all_sdm_richness <- ggarrange(p_glm5E,p_glm2.5E,p_bart5E,p_bart5E,p_fne5E,p_fne5E,
                  labels = c("GLM(5E)","GLM(2.5E)","BART(5E)","BART(2.5E)","FNE(5E)","FNE(2.5E"), 
                  ncol = 2, nrow = 3,common.legend = TRUE, legend="right")
all_sdm_richness
ggsave(filename="Fig_Revision_Xb_sp_rich(SDM).tiff",
       plot=all_sdm_richness,
       device='tiff',
       width=210,
       height=297,
       units="mm")       



