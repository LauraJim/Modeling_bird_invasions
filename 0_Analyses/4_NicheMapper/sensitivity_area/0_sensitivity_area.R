# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# LOAD PACKAGES ####

library(terra)
library(modEvA)
library(kuenm)
library(raster)
library(enmSdm)
#library(ecospat)
library(foreach)
library(doParallel)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(ggeffects)
library(ggpubr)

# GET TARGET SPECIES ####

occ_files <- list.files("../../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

n.iters <- data.frame(c(1,1,3,3,3,6,2,3,3,1,3,5,1,2,3,1,3,3,2,3))
  colnames(n.iters) <- "iter"
  n.iters$species <- species_list

# READ IN SENSITIVITY FILES
data_area <- data.frame(matrix(ncol = 24, nrow = 0))
for (species in species_list){
  iters <- n.iters[n.iters$species ==species,]$iter
  tempspecies <- data.frame(matrix(ncol = 24, nrow = 0))
  for (i in 1:iters){
  param <- read.table(paste0("./", species,"/",i,"/","parameters.txt"),sep="\t",h=T)
    param$iter <- rep(i,nrow(param))
  vals <- read.table(paste0("./", species,"/",i,"/","all.sensi.txt"),sep="\t",h=T)
    vals$iter <- rep(i,nrow(param))
  temp <- cbind(param,vals)
  temp$species <- rep(species,nrow(temp))
  tempspecies <- rbind(tempspecies,temp)
  }
  data_area <- rbind(data_area,tempspecies)
  }
  
head(data_area)
tail(data_area)

# ANALYSIS: MIXED BETARAGRESSION
  data_area$speci.MN <- replace(data_area$speci.MN, data_area$speci.MN ==0,0.001) 
  data_area$speci.MN <- replace(data_area$speci.MN, data_area$speci.MN ==1,0.999) 
  data_area$speci.MN <- (data_area$speci.MN -1)*-1 

  m <- lmer(BMR ~body.mass + (1|species),data=data_area)
  data_area$BMRres <- residuals(m)

data_area$Tb <- scale(data_area$Tb)
data_area$body.mass <- scale(data_area$body.mass)
data_area$fd <- scale(data_area$fd)
data_area$fl <- scale(data_area$fl)
data_area$BMRres <- scale(data_area$BMRres)

m1 <- glmmTMB(speci.MN ~Tb + body.mass  + fd + fl + BMRres +(1|species/iter), data_area, family="beta_family")
  summary(m1)
  #ggplot(data_area,aes(x=Tb,y=speci.MN))+geom_point()+geom_smooth(method=loess,color='red')
  
  
## BMR
  mydf <- ggpredict(m1, terms = "BMRres",allow.new.levels=TRUE)
  fm1 <- ggplot(mydf,aes(x,predicted)) 
    fm1 <- fm1 + geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
    scale_y_continuous(limits = c(0,1), oob = scales::squish) +
    ylab("% Europe predicted as suitable\n") + theme(axis.title.y = element_text(size = 15))+
    xlab("\n BMR") + theme(axis.title.x = element_text(size = 15)) + theme_classic()
  fm1

## BODY MASS
  mydf <- ggpredict(m1, terms = "body.mass",allow.new.levels=TRUE)
  fm2 <- ggplot(mydf,aes(x,predicted)) 
  fm2 <- fm2 + geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
    scale_y_continuous(limits = c(0,1), oob = scales::squish) +
    ylab("% Europe predicted as suitable\n") + theme(axis.title.y = element_text(size = 15))+
    xlab("\n body mass") + theme(axis.title.x = element_text(size = 15)) + theme_classic()
  fm2  
  
## FEATHER DEPTH
  mydf <- ggpredict(m1, terms = "fd",allow.new.levels=TRUE)
  fm3 <- ggplot(mydf,aes(x,predicted)) 
  fm3 <- fm3 + geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
    scale_y_continuous(limits = c(0,1), oob = scales::squish) +
    ylab("% Europe predicted as suitable\n") + theme(axis.title.y = element_text(size = 15))+
    xlab("\n feather depth") + theme(axis.title.x = element_text(size = 15)) + theme_classic()
  fm3   
  
## FEATHER LENGTH
  mydf <- ggpredict(m1, terms = "fl",allow.new.levels=TRUE)
  fm4 <- ggplot(mydf,aes(x,predicted)) 
  fm4 <- fm4 + geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
    scale_y_continuous(limits = c(0,1), oob = scales::squish) +
    ylab("% Europe predicted as suitable\n") + theme(axis.title.y = element_text(size = 15))+
    xlab("\n feather length") + theme(axis.title.x = element_text(size = 15)) + theme_classic()
  fm4
  
## BODY TEMPERATURE
  mydf <- ggpredict(m1, terms = "Tb",allow.new.levels=TRUE)
  fm5 <- ggplot(mydf,aes(x,predicted)) 
  fm5 <- fm5 + geom_line(linewidth = 1.5) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
    scale_y_continuous(limits = c(0,1), oob = scales::squish) +
    ylab("% Europe predicted as suitable\n") + theme(axis.title.y = element_text(size = 15))+
    xlab("\n body temperature") + theme(axis.title.x = element_text(size = 15)) + theme_classic()
  fm5  
  
############################################################      
###COMBINE IN A SINGLE PLOT AND STANDARDIZE RANGE 0 to 20###
############################################################
all_figs <- ggarrange(fm1,fm2,fm3,fm4,fm5,
                                labels = c("A","B","C","D","E","F"), 
                                ncol = 2, nrow = 3,common.legend = TRUE, legend="right")
all_figs

  ggsave(filename="all_figs_sensitivity_area.png",
         plot=all_figs,
         device='png',
         width=210,
         height=297,
         units="mm")   
  
  
