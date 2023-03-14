# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(lme4)
library(lmerTest)
library(effects)
library(reshape2)
library(dplyr)
library(glmmTMB)
library(emmeans)

# GET TARGET SPECIES ####

occ_files <- list.files("../0_Nf_modeling/occurrences-v3")
occ_files

species_list <- unique(sapply(strsplit(basename(tools::file_path_sans_ext(occ_files)), "_"), `[`, 1))
species_list

###########################################################################################################################################################
####AUC####################################################################################################################################################
###########################################################################################################################################################
#2.5E
#GET AUC data
#CLIMATE
AUC_native_clim_sdm <- data.frame(matrix(nrow = 0, ncol = 5)) 
for (species in species_list) {
 AUCs <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(5,10)]
  AUCs <- data.frame(t(apply(AUCs,2,mean,na.rm=TRUE)))
  AUCs$preds <- "clim"
  AUCs$species <- species
  AUCs$range <- "native"
  AUC_native_clim_sdm <- rbind(AUC_native_clim_sdm,AUCs)
}
  AUC_native_clim_sdm$threshold <- rep(2.5,nrow(AUC_native_clim_sdm))
  
AUC_native_clim_fne <- data.frame(matrix(nrow = 0, ncol = 9)) 
for (species in species_list) {
  AUCs <- read.csv(paste0("../0_Nf_modeling/eval_metrics/NATIVE/BIAS/2.5E/", "fne_native_",species,".csv"))
  AUC_native_clim_fne <- rbind(AUC_native_clim_fne,AUCs)
  }
  AUC_native_clim_fne$algorithm <- rep("fne",nrow(AUC_native_clim_fne))


#CLIMATE_BIOTIC_HABITAT
AUC_native_clim_hb_sdm <- data.frame(matrix(nrow = 0, ncol = 5)) 
for (species in species_list) {
  AUCs <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(5,10)]
  AUCs <- data.frame(t(apply(AUCs,2,mean,na.rm=TRUE)))
  AUCs$preds <- "clim_hb"
  AUCs$species <- species
  AUCs$range <- "native"
  AUC_native_clim_hb_sdm <- rbind(AUC_native_clim_hb_sdm,AUCs)
}
  AUC_native_clim_hb_sdm$threshold <- rep(2.5,nrow(AUC_native_clim_hb_sdm))

  
AUC_native_clim_hb_fne <- data.frame(matrix(nrow = 0, ncol = 9)) 
  for (species in species_list) {
    AUCs <- read.csv(paste0("../0_Nf_modeling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/2.5E/", "fne_native_",species,".csv"))
    AUC_native_clim_hb_fne <- rbind(AUC_native_clim_hb_fne,AUCs)
  }  
  AUC_native_clim_hb_fne$algorithm <- rep("fne",nrow(AUC_native_clim_hb_fne))
  
AUCnative <- rbind(AUC_native_clim_sdm,AUC_native_clim_hb_sdm)
  head(AUCnative)
  AUCnative <- melt(AUCnative,id.vars=c("preds","species","range","threshold"))
  colnames(AUCnative) <- c("preds","species","range","threshold","algorithm","AUC")
  head(AUCnative)
  
  temp_fne <- rbind(AUC_native_clim_fne ,AUC_native_clim_hb_fne )
  temp_fne <- temp_fne[c(6:10,1)]
  colnames(temp_fne) <- c("preds","species","range","threshold","algorithm","AUC")
  
  AUCnative_2.5E <- rbind(AUCnative,temp_fne)
  
#5E
#GET AUC data
#CLIMATE
AUC_native_clim_sdm <- data.frame(matrix(nrow = 0, ncol = 5)) 
for (species in species_list) {
    AUCs <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(5,10)]
    AUCs <- data.frame(t(apply(AUCs,2,mean,na.rm=TRUE)))
    AUCs$preds <- "clim"
    AUCs$species <- species
    AUCs$range <- "native"
    AUC_native_clim_sdm <- rbind(AUC_native_clim_sdm,AUCs)
  }
  AUC_native_clim_sdm$threshold <- rep(5,nrow(AUC_native_clim_sdm))
  
AUC_native_clim_fne <- data.frame(matrix(nrow = 0, ncol = 9)) 
for (species in species_list) {
    AUCs <- read.csv(paste0("../0_Nf_modeling/eval_metrics/NATIVE/BIAS/5E/", "fne_native_",species,".csv"))
    AUC_native_clim_fne <- rbind(AUC_native_clim_fne,AUCs)
  }
AUC_native_clim_fne$algorithm <- rep("fne",nrow(AUC_native_clim_hb_fne))
  
  
#CLIMATE_BIOTIC_HABITAT
AUC_native_clim_hb_sdm <- data.frame(matrix(nrow = 0, ncol = 5)) 
for (species in species_list) {
  AUCs <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(5,10)]
  AUCs <- data.frame(t(apply(AUCs,2,mean,na.rm=TRUE)))
  AUCs$preds <- "clim_hb"
  AUCs$species <- species
  AUCs$range <- "native"
  AUC_native_clim_hb_sdm <- rbind(AUC_native_clim_hb_sdm,AUCs)
}
  AUC_native_clim_hb_sdm$threshold <- rep(5,nrow(AUC_native_clim_hb_sdm))
  
  
AUC_native_clim_hb_fne <- data.frame(matrix(nrow = 0, ncol = 9)) 
for (species in species_list) {
    AUCs <- read.csv(paste0("../0_Nf_modeling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/", "fne_native_",species,".csv"))
    AUC_native_clim_hb_fne <- rbind(AUC_native_clim_hb_fne,AUCs)
}  
AUC_native_clim_hb_fne$algorithm <- rep("fne",nrow(AUC_native_clim_hb_fne))
  
AUCnative <- rbind(AUC_native_clim_sdm,AUC_native_clim_hb_sdm)
head(AUCnative)
AUCnative <- melt(AUCnative,id.vars=c("preds","species","range","threshold"))
colnames(AUCnative) <- c("preds","species","range","threshold","algorithm","AUC")
head(AUCnative)
  
temp_fne <- rbind(AUC_native_clim_fne ,AUC_native_clim_hb_fne )
temp_fne <- temp_fne[c(6:10,1)]
colnames(temp_fne) <- c("preds","species","range","threshold","algorithm","AUC")
  
AUCnative_5E <- rbind(AUCnative,temp_fne)

AUCnative <- rbind(AUCnative_5E,AUCnative_2.5E)

#
##COMPARE NATIVE RANGE AUCs##
#
AUCnative$threshold <- as.factor(AUCnative$threshold)
m1 <- lmer(AUC ~(preds+threshold+algorithm)^3 + (1|species),data=AUCnative)
m1 <- lmer(AUC ~(preds+threshold+algorithm)^2 + (1|species),data=AUCnative)
m1 <- lmer(AUC ~preds+threshold+algorithm + preds:algorithm + (1|species),data=AUCnative)
  anova(m1)
  plot(allEffects(m1))
  allEffects(m1)
  emmeans (m1,  ~ threshold)
    pairs(emmeans (m1,  ~ threshold))
  emmeans (m1,  ~ preds|algorithm)
    pairs(emmeans (m1,  ~ preds|algorithm))
  emmeans (m1,  ~ algorithm|preds)
    pairs(emmeans (m1,  ~ algorithm|preds))
    
    #get raw mean and SD
    five_perc <- AUCnative[AUCnative$threshold ==5,]
      mean(five_perc$AUC)
      sd(five_perc$AUC)
    two.five_perc <- AUCnative[AUCnative$threshold ==2.5,]
      mean(two.five_perc$AUC)
      sd(two.five_perc$AUC)
      
    notFNE <- AUCnative[!AUCnative$algorithm=="fne",] 
      #CLIM
      notFNEclim <- notFNE[notFNE$preds=="clim",]
        notFNEclimBART <- notFNEclim[notFNEclim$algorithm=="AUCratio_bart",]
          mean(notFNEclimBART$AUC)
          sd(notFNEclimBART$AUC)
        notFNEclimGLM <- notFNEclim[notFNEclim$algorithm=="AUCratio_glm",]
          mean(notFNEclimGLM$AUC)
          sd(notFNEclimGLM$AUC)
      #CBH
      notFNEchb <- notFNE[notFNE$preds=="clim_hb",]
        notFNEchbBART <- notFNEchb[notFNEchb$algorithm=="AUCratio_bart",]
          mean(notFNEchbBART$AUC)
          sd(notFNEchbBART$AUC)
        notFNEchbGLM <- notFNEchb[notFNEchb$algorithm=="AUCratio_glm",]
          mean(notFNEchbGLM$AUC)
          sd(notFNEchbGLM$AUC)
          
      FNE <- AUCnative[AUCnative$algorithm=="fne",]
        FNEclim <- FNE[FNE$preds=="clim",]
          mean(FNEclim$AUC)
          sd(FNEclim$AUC)
        FNEchb <- FNE[FNE$preds=="clim_hb",]
          mean(FNEchb$AUC)
          sd(FNEchb$AUC)

#############################################################################################################################################################
#GET SPEC data###############################################################################################################################################
#############################################################################################################################################################    
#2.5E
SPEC_native_clim <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SPEC <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(9,14)]
      SPEC <- data.frame(t(apply(SPEC,2,mean,na.rm=TRUE)))
      SPEC$preds <- "clim"
      SPEC$species <- species
      SPEC$range <- "native"
      SPEC_native_clim <- rbind(SPEC_native_clim,SPEC)
    }
    
SPEC_native_clim_hb <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SPEC <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(9,14)]
      SPEC <- data.frame(t(apply(SPEC,2,mean,na.rm=TRUE)))
      SPEC$preds <- "clim_hb"
      SPEC$species <- species
      SPEC$range <- "native"
      SPEC_native_clim_hb <- rbind(SPEC_native_clim_hb,SPEC)
    }
    
  head(SPEC_native_clim)
  head(SPEC_native_clim_hb)
    
SPECnative2.5E <- rbind(SPEC_native_clim,SPEC_native_clim_hb)
    head(SPECnative2.5E)
    SPECnative2.5E <- melt(SPECnative2.5E,id.vars=c("preds","species","range"))
    colnames(SPECnative2.5E) <- c("preds","species","range","algorithm","SPEC")
    SPECnative2.5E$threshold <- rep(2.5,nrow(SPECnative2.5E))
    head(SPECnative2.5E)
    
#5E
SPEC_native_clim <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SPEC <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(9,14)]
      SPEC <- data.frame(t(apply(SPEC,2,mean,na.rm=TRUE)))
      SPEC$preds <- "clim"
      SPEC$species <- species
      SPEC$range <- "native"
      SPEC_native_clim <- rbind(SPEC_native_clim,SPEC)
    }
    
SPEC_native_clim_hb <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SPEC <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(9,14)]
      SPEC <- data.frame(t(apply(SPEC,2,mean,na.rm=TRUE)))
      SPEC$preds <- "clim_hb"
      SPEC$species <- species
      SPEC$range <- "native"
      SPEC_native_clim_hb <- rbind(SPEC_native_clim_hb,SPEC)
    }
    
head(SPEC_native_clim)
    head(SPEC_native_clim_hb)
    
    SPECnative5E <- rbind(SPEC_native_clim,SPEC_native_clim_hb)
    head(SPECnative5E)
    SPECnative5E <- melt(SPECnative5E,id.vars=c("preds","species","range"))
    colnames(SPECnative5E) <- c("preds","species","range","algorithm","SPEC")
    SPECnative5E$threshold <- rep(5,nrow(SPECnative5E))
    head(SPECnative5E)
    
SPECnative <- rbind(SPECnative2.5E,SPECnative5E)
    SPECnative$SPEC <- replace(SPECnative$SPEC, SPECnative$SPEC ==0,0.001) 
    SPECnative$SPEC <- replace(SPECnative$SPEC, SPECnative$SPEC ==1,0.999) 

SPECnative$threshold <- as.factor(SPECnative$threshold)        
m2 <- glmmTMB(SPEC ~ (preds+threshold+algorithm)^3 + (1|species), SPECnative, family=beta_family(link="logit"))
m2 <- glmmTMB(SPEC ~ (preds+threshold+algorithm)^2 + (1|species), SPECnative, family=beta_family(link="logit"))
m2 <- glmmTMB(SPEC ~ preds+threshold+algorithm+preds:algorithm + (1|species), SPECnative, family=beta_family(link="logit"))
  car::Anova(m2) 
  plot(allEffects(m2))
    allEffects(m2)
    emmeans (m2,  ~ preds|algorithm,type="response")
    pairs(emmeans (m2,  ~ preds|algorithm,type="response"))

    notFNE <- SPECnative[!SPECnative$algorithm=="fne",] 
    #CLIM
    notFNEclim <- notFNE[notFNE$preds=="clim",]
    notFNEclimBART <- notFNEclim[notFNEclim$algorithm=="spec_bart",]
      mean(notFNEclimBART$SPEC)
      sd(notFNEclimBART$SPEC)
    notFNEclimGLM <- notFNEclim[notFNEclim$algorithm=="spec_glm",]
      mean(notFNEclimGLM$SPEC)
      sd(notFNEclimGLM$SPEC)
    #CBH
    notFNEchb <- notFNE[notFNE$preds=="clim_hb",]
    notFNEchbBART <- notFNEchb[notFNEchb$algorithm=="spec_bart",]
      mean(notFNEchbBART$SPEC)
      sd(notFNEchbBART$SPEC)
    notFNEchbGLM <- notFNEchb[notFNEchb$algorithm=="spec_glm",]
      mean(notFNEchbGLM$SPEC)
      sd(notFNEchbGLM$SPEC)
    
    
    
#############################################################################################################################################################
#GET SENS data###############################################################################################################################################
#############################################################################################################################################################    
#2.5E
SENS_native_clim <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SENS <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(8,13)]
      SENS <- data.frame(t(apply(SENS,2,mean,na.rm=TRUE)))
      SENS$preds <- "clim"
      SENS$species <- species
      SENS$range <- "native"
      SENS_native_clim <- rbind(SENS_native_clim,SENS)
    }
    
SENS_native_clim_hb <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SENS <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/2.5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(8,13)]
      SENS <- data.frame(t(apply(SENS,2,mean,na.rm=TRUE)))
      SENS$preds <- "clim_hb"
      SENS$species <- species
      SENS$range <- "native"
      SENS_native_clim_hb <- rbind(SENS_native_clim_hb,SENS)
    }
    
  head(SENS_native_clim)
  head(SENS_native_clim_hb)
    
  SENSnative2.5E <- rbind(SENS_native_clim,SENS_native_clim_hb)
  head(SENSnative2.5E)
  SENSnative2.5E <- melt(SENSnative2.5E,id.vars=c("preds","species","range"))
    colnames(SENSnative2.5E) <- c("preds","species","range","algorithm","SENS")
    SENSnative2.5E$threshold <- rep(2.5,nrow(SENSnative2.5E))
    head(SENSnative2.5E)
    
#5E
SENS_native_clim <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SENS <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(8,13)]
      SENS <- data.frame(t(apply(SENS,2,mean,na.rm=TRUE)))
      SENS$preds <- "clim"
      SENS$species <- species
      SENS$range <- "native"
      SENS_native_clim <- rbind(SENS_native_clim,SENS)
    }
    
SENS_native_clim_hb <- data.frame(matrix(nrow = 0, ncol = 5)) 
    for (species in species_list) {
      SENS <- read.csv(paste0("../1_PA_modelling/eval_metrics/NATIVE/BIAS_BIOTIC_HABITAT/5E/crossval_metrics/", "crossval_metrics_nat_",species,".csv"))[c(8,13)]
      SENS <- data.frame(t(apply(SENS,2,mean,na.rm=TRUE)))
      SENS$preds <- "clim_hb"
      SENS$species <- species
      SENS$range <- "native"
      SENS_native_clim_hb <- rbind(SENS_native_clim_hb,SENS)
    }
    
head(SENS_native_clim)
head(SENS_native_clim_hb)
    
SENSnative5E <- rbind(SENS_native_clim,SENS_native_clim_hb)
  head(SENSnative5E)
    SENSnative5E <- melt(SENSnative5E,id.vars=c("preds","species","range"))
    colnames(SENSnative5E) <- c("preds","species","range","algorithm","SENS")
    SENSnative5E$threshold <- rep(5,nrow(SENSnative5E))
    head(SENSnative5E)
    
    SENSnative <- rbind(SENSnative2.5E,SENSnative5E)
    SENSnative$SENS <- replace(SENSnative$SENS, SENSnative$SENS ==0,0.001) 
    SENSnative$SENS <- replace(SENSnative$SENS, SENSnative$SENS ==1,0.999) 

SENSnative$threshold <- as.factor(SENSnative$threshold )        
m3 <- glmmTMB(SENS ~ (preds+threshold+algorithm)^3 + (1|species), SENSnative, family=beta_family(link="logit"))
    m3 <- glmmTMB(SENS ~ (preds+threshold+algorithm)^2 + (1|species), SENSnative, family=beta_family(link="logit"))
    m3 <- glmmTMB(SENS ~ preds+threshold+algorithm+preds:algorithm + (1|species), SENSnative, family=beta_family(link="logit"))
    car::Anova(m3) 
    plot(allEffects(m3))
    allEffects(m3)
    emmeans (m3,  ~ preds|algorithm,type="response")
    pairs(emmeans (m3,  ~ preds|algorithm,type="response"))
    
    
    notFNE <- SENSnative[!SENSnative$algorithm=="fne",] 
    #CLIM
    notFNEclim <- notFNE[notFNE$preds=="clim",]
      notFNEclimBART <- notFNEclim[notFNEclim$algorithm=="sens_bart",]
        mean(notFNEclimBART$SENS)
        sd(notFNEclimBART$SENS)
      notFNEclimGLM <- notFNEclim[notFNEclim$algorithm=="sens_glm",]
        mean(notFNEclimGLM$SENS)
        sd(notFNEclimGLM$SENS)
    #CBH
    notFNEchb <- notFNE[notFNE$preds=="clim_hb",]
      notFNEchbBART <- notFNEchb[notFNEchb$algorithm=="sens_bart",]
        mean(notFNEchbBART$SENS)
        sd(notFNEchbBART$SENS)
    notFNEchbGLM <- notFNEchb[notFNEchb$algorithm=="sens_glm",]
        mean(notFNEchbGLM$SENS)
        sd(notFNEchbGLM$SENS)
    
    
    
    
