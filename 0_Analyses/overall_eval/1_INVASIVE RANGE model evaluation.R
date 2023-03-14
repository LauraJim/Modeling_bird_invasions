# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#libraries
library(reshape2)
library(lme4)
library(lmerTest)
library(effects)
library(glmmTMB)
  
###########################################################################################################################################################################
###READ IN MODELS WITH CLIMATE BIOTIC AND HABITAT VARIABLES################################################################################################################
###########################################################################################################################################################################

bias_clamping_5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/5E/", "eval_metrics_eur_full.csv"))
bias_clamping_5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/5E/", "eval_metrics_eur_mop.csv"))
bias_clamping_5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/5E/", "eval_metrics_inv_full.csv"))
bias_clamping_5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/5E/", "eval_metrics_inv_mop.csv"))

bias_clamping_2.5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/2.5E/", "eval_metrics_eur_full.csv"))
bias_clamping_2.5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/2.5E/", "eval_metrics_eur_mop.csv"))
bias_clamping_2.5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/2.5E/", "eval_metrics_inv_full.csv"))
bias_clamping_2.5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/CLAMPING/2.5E/", "eval_metrics_inv_mop.csv"))

bias_extrapol_5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/5E/", "eval_metrics_eur_full.csv"))
bias_extrapol_5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/5E/", "eval_metrics_eur_mop.csv"))
bias_extrapol_5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/5E/", "eval_metrics_inv_full.csv"))
bias_extrapol_5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/5E/", "eval_metrics_inv_mop.csv"))

bias_extrapol_2.5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/2.5E/", "eval_metrics_eur_full.csv"))
bias_extrapol_2.5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/2.5E/", "eval_metrics_eur_mop.csv"))
bias_extrapol_2.5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/2.5E/", "eval_metrics_inv_full.csv"))
bias_extrapol_2.5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS_BIOTIC_HABITAT/EXTRAPOL/2.5E/", "eval_metrics_inv_mop.csv"))

invaded <- c(rep('whole_europe',40),rep('dispersal',40))
invaded <- c(invaded,invaded,invaded,invaded)



eval_data <- rbind(bias_clamping_5E_eur_full,bias_clamping_5E_eur_mop,bias_clamping_5E_inv_full,bias_clamping_5E_inv_mop,
                   bias_clamping_2.5E_eur_full,bias_clamping_2.5E_eur_mop,bias_clamping_2.5E_inv_full,bias_clamping_2.5E_inv_mop,
                   bias_extrapol_5E_eur_full,bias_extrapol_5E_eur_mop,bias_extrapol_5E_inv_full,bias_extrapol_5E_inv_mop,
                   bias_extrapol_2.5E_eur_full,bias_extrapol_2.5E_eur_mop,bias_extrapol_2.5E_inv_full,bias_extrapol_2.5E_inv_mop)
                   
eval_data <- eval_data[c(1:4,11:13,14:16,17:21)] 
  eval_data$threshold <- as.factor(eval_data$threshold)
  eval_data$invaded <- invaded
  unique(eval_data$BIAS_BIOTIC_HABITAT)
  eval_data$preds <- rep("climate_biot_hab",nrow(eval_data))
  unique(eval_data$extrapol)
  unique(eval_data$threshold)
  unique(eval_data$mop)
  unique(eval_data$invaded)
  eval_data_cbh <- eval_data
  
###########################################################################################################################################################################
###READ IN MODELS WITH CLIMATE ONLYIC AND HABITAT VARIABLES################################################################################################################
###########################################################################################################################################################################
bias_clamping_5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/5E/", "eval_metrics_eur_full.csv"))
bias_clamping_5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/5E/", "eval_metrics_eur_mop.csv"))
bias_clamping_5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/5E/", "eval_metrics_inv_full.csv"))
bias_clamping_5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/5E/", "eval_metrics_inv_mop.csv"))
  
bias_clamping_2.5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/", "eval_metrics_eur_full.csv"))
bias_clamping_2.5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/", "eval_metrics_eur_mop.csv"))
bias_clamping_2.5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/", "eval_metrics_inv_full.csv"))
bias_clamping_2.5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/CLAMPING/2.5E/", "eval_metrics_inv_mop.csv"))
  
bias_extrapol_5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/5E/", "eval_metrics_eur_full.csv"))
bias_extrapol_5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/5E/", "eval_metrics_eur_mop.csv"))
bias_extrapol_5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/5E/", "eval_metrics_inv_full.csv"))
bias_extrapol_5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/5E/", "eval_metrics_inv_mop.csv"))

bias_extrapol_2.5E_eur_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/2.5E/", "eval_metrics_eur_full.csv"))
bias_extrapol_2.5E_eur_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/2.5E/", "eval_metrics_eur_mop.csv"))
bias_extrapol_2.5E_inv_full <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/2.5E/", "eval_metrics_inv_full.csv"))
bias_extrapol_2.5E_inv_mop <- read.csv(paste0("../1_PA_modelling/eval_metrics/INVASIVE/BIAS/EXTRAPOL/2.5E/", "eval_metrics_inv_mop.csv"))
  
invaded <- c(rep('whole_europe',40),rep('dispersal',40))
invaded <- c(invaded,invaded,invaded,invaded)
  
 eval_data <- rbind(bias_clamping_5E_eur_full,bias_clamping_5E_eur_mop,bias_clamping_5E_inv_full,bias_clamping_5E_inv_mop,
                     bias_clamping_2.5E_eur_full,bias_clamping_2.5E_eur_mop,bias_clamping_2.5E_inv_full,bias_clamping_2.5E_inv_mop,
                     bias_extrapol_5E_eur_full,bias_extrapol_5E_eur_mop,bias_extrapol_5E_inv_full,bias_extrapol_5E_inv_mop,
                     bias_extrapol_2.5E_eur_full,bias_extrapol_2.5E_eur_mop,bias_extrapol_2.5E_inv_full,bias_extrapol_2.5E_inv_mop)
  
  eval_data <- eval_data[c(1:4,11:13,14:16,17:21)] 
  eval_data$threshold <- as.factor(eval_data$threshold)
  eval_data$invaded <- invaded
  unique(eval_data$BIAS)
  eval_data$preds <- rep("climate",nrow(eval_data))
  unique(eval_data$extrapol)
  unique(eval_data$threshold)
  unique(eval_data$mop)
  unique(eval_data$invaded)
  colnames(eval_data)
  
  eval_data_c <- eval_data
  
#merge datasets
  colnames(eval_data_cbh)
    colnames(eval_data_cbh)[11] <- "BIAS"
  colnames(eval_data_c)
    colnames(eval_data_c)[11] <- "BIAS"
  
eval_data <- rbind(eval_data_cbh,eval_data_c)
  colnames(eval_data)
  
  
AUCratio <- melt(eval_data[c(1:4,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))
sens <- melt(eval_data[c(1,5:7,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))
spec <- melt(eval_data[c(1,8:10,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))

#####################################################################################################################################
###ANALYSES BASED ON CLIMATE BIOTIC AND HABITAT######################################################################################
#####################################################################################################################################
 
#
###AUC ratio##########################################################################################
#
colnames(AUCratio) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","AUC")
head(AUCratio)

m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^5 + (1|species),data=AUCratio)
  anova(m_AUC) 

m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^4 + (1|species),data=AUCratio)
  anova(m_AUC) 

m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^3 + (1|species),data=AUCratio)
  anova(m_AUC)  

m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^2 + (1|species),data=AUCratio)
  anova(m_AUC)  
  
m_AUC <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded + preds 
              + extrapol:preds
              + algorithm:extrapol
              + algorithm:mop
              + algorithm:preds + (1|species),data=AUCratio)
  anova(m_AUC)
  p3 <- plot(allEffects(m_AUC),select =3)
  p4 <- plot(allEffects(m_AUC),select =4)
  p5 <- plot(allEffects(m_AUC),select =5)
  p6 <- plot(allEffects(m_AUC),select =6)
  library(gridExtra)
  grid.arrange(p3,p4,p5,p6, nrow=2,  ncol=2)
  allEffects(m_AUC)
  
  emmeans(m_AUC,  ~ preds|algorithm)
    pairs(emmeans(m_AUC,  ~ preds|algorithm))
  
  #get raw summaries
  #BART
  mean(eval_data_c$AUCratio_bart)
  sd(eval_data_c$AUCratio_bart)
  
  mean(eval_data_cbh$AUCratio_bart,na.rm=TRUE)
  sd(eval_data_cbh$AUCratio_bart,na.rm=TRUE)
  
  #GLM
  mean(eval_data_c$AUCratio_glm)
  sd(eval_data_c$AUCratio_glm)
  
  mean(eval_data_cbh$AUCratio_glm,na.rm=TRUE)
  sd(eval_data_cbh$AUCratio_glm,na.rm=TRUE)
  
  #FNE
  mean(eval_data_c$AUCratio_fne,na.rm=TRUE)
  sd(eval_data_c$AUCratio_fne,na.rm=TRUE)
  
  mean(eval_data_cbh$AUCratio_fne,na.rm=TRUE)
  sd(eval_data_cbh$AUCratio_fne,na.rm=TRUE)
  
#
###SENSITIVITY##########################################################################################
#
colnames(sens) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","sens")
  head(sens)
  sens$sens <- replace(sens$sens, sens$sens ==0,0.001) 
  sens$sens <- replace(sens$sens, sens$sens ==1,0.999) 
  
m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded + preds)^5+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded + preds)^4+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)
  
m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded + preds)^3+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm +invaded + preds)^3+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm +invaded + preds)^2+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)    
  
m_sens <- glmmTMB(sens ~ extrapol +  threshold + algorithm +invaded + preds + algorithm:preds + (1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)    
  plot(allEffects(m_sens))
  allEffects(m_sens)

#
###SPECIFICITY##########################################################################################
#
colnames(spec) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","spec")
head(spec)
  spec$spec <- replace(spec$spec, spec$spec ==0,0.001) 
  spec$spec <- replace(spec$spec, spec$spec ==1,0.999) 
  
m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^5+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)
  
m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^4+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)
  
m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded + preds)^3+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)  

m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded + preds )^2+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)     
  
m_spec <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded + preds + algorithm:preds+ (1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)    
  plot(allEffects(m_spec))
  allEffects(m_spec)

  
#####################################################################################################################################
###ANALYSES BASED ON CLIMATE ONLY####################################################################################################
#####################################################################################################################################
 
AUCratio <- melt(eval_data_c[c(1:4,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))
sens <- melt(eval_data_c[c(1,5:7,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))
spec <- melt(eval_data_c[c(1,8:10,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))  

#
###AUC ratio##########################################################################################
#
  colnames(AUCratio) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","AUC")
  head(AUCratio)
  
  m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded)^5 + (1|species),data=AUCratio)
  anova(m_AUC) 
  
  m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded)^4 + (1|species),data=AUCratio)
  anova(m_AUC) 
  
  m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded )^3 + (1|species),data=AUCratio)
  anova(m_AUC)  
  
  m_AUC <- lmer(AUC ~ (extrapol +  threshold + algorithm + mop + invaded)^2 + (1|species),data=AUCratio)
  anova(m_AUC)  
  
  m_AUC <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded 
                + mop:algorithm
                + extrapol:algorithm
                + (1|species),data=AUCratio)
  anova(m_AUC)
  plot(allEffects(m_AUC))
  allEffects(m_AUC)
  
  #get p-values
  pairs(emmeans(m_AUC,  ~ mop|algorithm))
  
  pairs(emmeans(m_AUC,  ~ extrapol|algorithm))
  
  
  #get some raw summaries
  
  clim_mop_yes <- eval_data_c[eval_data_c$mop=="yes_mop",]
  clim_mop_no <- eval_data_c[eval_data_c$mop=="no",]
  #GLM
  mean(clim_mop_yes$AUCratio_glm)
  sd(clim_mop_yes$AUCratio_glm)
  mean(clim_mop_no$AUCratio_glm,na.rm=TRUE)
  sd(clim_mop_no$AUCratio_glm,na.rm=TRUE)
  #BART
  mean(clim_mop_yes$AUCratio_bart)
  sd(clim_mop_yes$AUCratio_bart)
  mean(clim_mop_no$AUCratio_bart,na.rm=TRUE)
  sd(clim_mop_no$AUCratio_bart,na.rm=TRUE)
  #FNE
  mean(clim_mop_yes$AUCratio_fne,na.rm=TRUE)
  sd(clim_mop_yes$AUCratio_fne,na.rm=TRUE)
  mean(clim_mop_no$AUCratio_fne,na.rm=TRUE)
  sd(clim_mop_no$AUCratio_fne,na.rm=TRUE)
  
  clim_extrapol_yes <- eval_data_c[eval_data_c$extrapol=="yes",]
  clim_extrapol_no <- eval_data_c[eval_data_c$extrapol=="clamping",]
  #GLM
  mean(clim_extrapol_yes$AUCratio_glm)
  sd(clim_extrapol_yes$AUCratio_glm)
  mean(clim_extrapol_no$AUCratio_glm,na.rm=TRUE)
  sd(clim_extrapol_no$AUCratio_glm,na.rm=TRUE)
  #BART
  mean(clim_extrapol_yes$AUCratio_bart)
  sd(clim_extrapol_yes$AUCratio_bart)
  mean(clim_extrapol_no$AUCratio_bart,na.rm=TRUE)
  sd(clim_extrapol_no$AUCratio_bart,na.rm=TRUE)
  #FNE
  mean(clim_extrapol_yes$AUCratio_fne,na.rm=TRUE)
  sd(clim_extrapol_yes$AUCratio_fne,na.rm=TRUE)
  mean(clim_extrapol_no$AUCratio_fne,na.rm=TRUE)
  sd(clim_extrapol_no$AUCratio_fne,na.rm=TRUE)
  
#
###SENSITIVITY##########################################################################################
#
colnames(sens) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","sens")
head(sens)
  sens$sens <- replace(sens$sens, sens$sens ==0,0.001) 
  sens$sens <- replace(sens$sens, sens$sens ==1,0.999) 
  
  m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded)^5+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
  m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded)^4+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)
  
  m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop +invaded)^3+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
  m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm +mop +invaded)^3+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)  
  
  m_sens <- glmmTMB(sens ~ (extrapol +  threshold + algorithm +invaded+mop)^2+(1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)    
  
  m_sens <- glmmTMB(sens ~ extrapol +  threshold + algorithm +invaded+mop + (1|species), sens, family=beta_family(link="logit"))
  car::Anova(m_sens)    
  plot(allEffects(m_sens))
  allEffects(m_sens)

  #comparisons
  pairs(emmeans(m_sens,  ~ threshold))
  pairs(emmeans(m_sens,  ~ algorithm))
  
  
  #get some raw data
  sens_glm <- sens[sens$algorithm=="sens_glm",]
    
  sens_bart <- sens[sens$algorithm=="sens_bart",]
  sens_fne <- sens[sens$algorithm=="sens_fne",]
  
  #glm
  mean(sens_glm[sens_glm$threshold==2.5,]$sens,na.rm=TRUE)
  sd(sens_glm[sens_glm$threshold==2.5,]$sens,na.rm=TRUE)
  
  mean(sens_glm[sens_glm$threshold==5,]$sens,na.rm=TRUE)
  sd(sens_glm[sens_glm$threshold==5,]$sens,na.rm=TRUE)
  
  #fne
  mean(sens_fne[sens_fne$threshold==2.5,]$sens,na.rm=TRUE)
  sd(sens_fne[sens_fne$threshold==2.5,]$sens,na.rm=TRUE)
  
  mean(sens_fne[sens_fne$threshold==5,]$sens,na.rm=TRUE)
  sd(sens_fne[sens_fne$threshold==5,]$sens,na.rm=TRUE)
  
  #bart
  mean(sens_bart[sens_bart$threshold==2.5,]$sens,na.rm=TRUE)
  sd(sens_bart[sens_bart$threshold==2.5,]$sens,na.rm=TRUE)
  
  mean(sens_bart[sens_bart$threshold==5,]$sens,na.rm=TRUE)
  sd(sens_bart[sens_bart$threshold==5,]$sens,na.rm=TRUE)
  
#
###SPECIFICITY##########################################################################################
#
  colnames(spec) <- c("species","bias","preds","extrapol","threshold","mop","invaded","algorithm","spec")
  head(spec)
  spec$spec <- replace(spec$spec, spec$spec ==0,0.001) 
  spec$spec <- replace(spec$spec, spec$spec ==1,0.999) 
  
  m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded)^5+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)
  
  m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded)^4+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)
  
  m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded)^3+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)  
  
  m_spec <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded)^2+(1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)     
  
  m_spec <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded + (1|species), spec, family=beta_family(link="logit"))
  car::Anova(m_spec)    
  plot(allEffects(m_spec))
  allEffects(m_spec)
  
  #comparisons
  pairs(emmeans(m_spec,  ~ threshold))
  pairs(emmeans(m_spec,  ~ algorithm))
  
  
  #get some raw data
  spec_glm <- spec[spec$algorithm=="spec_glm",]
  spec_bart <- spec[spec$algorithm=="spec_bart",]
  spec_fne <- spec[spec$algorithm=="spec_fne",]
  
  #glm
  mean(spec_glm[spec_glm$threshold==2.5,]$spec,na.rm=TRUE)
  sd(spec_glm[spec_glm$threshold==2.5,]$spec,na.rm=TRUE)
  
  mean(spec_glm[spec_glm$threshold==5,]$spec,na.rm=TRUE)
  sd(spec_glm[spec_glm$threshold==5,]$spec,na.rm=TRUE)
  
  #fne
  mean(spec_fne[spec_fne$threshold==2.5,]$spec,na.rm=TRUE)
  sd(spec_fne[spec_fne$threshold==2.5,]$spec,na.rm=TRUE)
  
  mean(spec_fne[spec_fne$threshold==5,]$spec,na.rm=TRUE)
  sd(spec_fne[spec_fne$threshold==5,]$spec,na.rm=TRUE)
  
  #bart
  mean(spec_bart[spec_bart$threshold==2.5,]$spec,na.rm=TRUE)
  sd(spec_bart[spec_bart$threshold==2.5,]$spec,na.rm=TRUE)
  
  mean(spec_bart[spec_bart$threshold==5,]$spec,na.rm=TRUE)
  sd(spec_bart[spec_bart$threshold==5,]$spec,na.rm=TRUE)
  
########################################################################################################################  
######TEST FOR INFLUENCE OF SAMPLE SIZE, NICHE TRUNCATION AND NICHE DYNAMICS ON INVASIVE RANGE EVALUATION STATISTICS####
########################################################################################################################
  
###READ IN REQUIRED DATA FOR TESTING  
  
#read in sample size data
  ss <- read.csv(paste0("../0_Nf_modeling/", "allspecies_samplesizes_v3.csv"))
  colnames(ss)[1] <- "species"
  
#read niche dynamics metrics
  nd <- read.table(paste0("../2_broenni/eu/", "niche_dynamics_EU_75perc.txt"),sep="\t",h=T)
  
#read niche truncation index
  nti <- read.table(paste0("../3_NTI/", "NTI.txt"),sep="\t",h=T)  

head(AUCratio)  
head(ss)
head(nd)
head(nti)
  

#############################################################################################################################################
###3################################# influence of native range sample size##################################################################
############################################################################################################################################# 

#AUC
eval_AUCratio2 <- merge(AUCratio,ss,by="species")
  #full model (=reduced model)
  m_AUCratio2 <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded  + N.rarefied +  (1|species),data=eval_AUCratio2)
  m_AUCratio2 <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded  + N.rarefied
                      + mop:algorithm
                      + extrapol:algorithm
                      + algorithm:N.rarefied
                      + (1|species),data=eval_AUCratio2)
  anova(m_AUCratio2)
  plot(allEffects(m_AUCratio2))
  
  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_fne",]
  m_fne.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=fne.alone)
  anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_glm",]
  m_glm.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.rarefied  + (1|species),data=glm.alone)
  anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_bart",]
  m_bart.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=bart.alone)
  anova(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#SENSITIVITY
eval_sens2 <- merge(sens,ss,by="species")
  #full model (=reduced model)
  m_sens2 <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop + invaded  + N.rarefied)^2 +  (1|species),data=eval_sens2,family=beta_family(link="logit"))
  car::Anova(m_sens2)
  m_sens2 <- glmmTMB(sens ~ extrapol +  threshold + algorithm + mop + invaded  + N.rarefied + algorithm:N.rarefied
                      + (1|species),data=eval_sens2,family=beta_family(link="logit"))
  car::Anova(m_sens2)
  plot(allEffects(m_sens2))  

  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_sens2[eval_sens2$algorithm =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_sens2[eval_sens2$algorithm =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_sens2[eval_sens2$algorithm =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#SPECIFICITY
  eval_spec2 <- merge(spec,ss,by="species")
  #full model (=reduced model)
  m_spec2 <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded  + N.rarefied)^2 +  (1|species),data=eval_spec2,family=beta_family(link="logit"))
  car::Anova(m_spec2)
  m_spec2 <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded  + N.rarefied
                     +algorithm:N.rarefied
                     + (1|species),data=eval_spec2,family=beta_family(link="logit"))
  car::Anova(m_spec2)
  plot(allEffects(m_spec2))  
  
  #explore interaction: sample size for each SDM separately
  #fne alone
  fne.alone <- eval_spec2[eval_spec2$algorithm =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_spec2[eval_spec2$algorithm =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_spec2[eval_spec2$algorithm =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.rarefied + (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  plot(allEffects(m_bart.alone))  
 
#############################################################################################################################################
###4################################# influence of invasive range sample size################################################################
############################################################################################################################################# 
  
#AUC
  eval_AUCratio2 <- merge(AUCratio,ss,by="species")
  #full model (=reduced model)
  m_AUCratio2 <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded  + N.inv.rarefied
                      + mop:algorithm
                      + extrapol:algorithm
                      + (1|species),data=eval_AUCratio2)
  anova(m_AUCratio2)
  plot(allEffects(m_AUCratio2))
  
  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_fne",]
  m_fne.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied +N.inv.rarefied + (1|species),data=fne.alone)
  anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_glm",]
  m_glm.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied  + (1|species),data=glm.alone)
  anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_AUCratio2[eval_AUCratio2$algorithm =="AUCratio_bart",]
  m_bart.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied  + (1|species),data=bart.alone)
  anova(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#SENSITIVITY
  eval_sens2 <- merge(sens,ss,by="species")
  #full model (=reduced model)
  m_sens2 <- glmmTMB(sens ~ (extrapol +  threshold + algorithm + mop + invaded  + N.inv.rarefied)^2 +  (1|species),data=eval_sens2,family=beta_family(link="logit"))
  car::Anova(m_sens2)
  m_sens2 <- glmmTMB(sens ~ extrapol +  threshold + algorithm + mop + invaded  + N.inv.rarefied
                     + (1|species),data=eval_sens2,family=beta_family(link="logit"))
  car::Anova(m_sens2)
  plot(allEffects(m_sens2))  
  
  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_sens2[eval_sens2$algorithm =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_sens2[eval_sens2$algorithm =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_sens2[eval_sens2$algorithm =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#SPECIFICITY
  eval_spec2 <- merge(spec,ss,by="species")
  #full model (=reduced model)
  m_spec2 <- glmmTMB(spec ~ (extrapol +  threshold + algorithm + mop + invaded  + N.inv.rarefied)^2 +  (1|species),data=eval_spec2,family=beta_family(link="logit"))
  car::Anova(m_spec2)
  m_spec2 <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded  + N.inv.rarefied
                     +algorithm:N.inv.rarefied
                     + (1|species),data=eval_spec2,family=beta_family(link="logit"))
  car::Anova(m_spec2)
  plot(allEffects(m_spec2))  
  
  #explore interaction: sample size for each SDM separately
  #fne alone
  fne.alone <- eval_spec2[eval_spec2$algorithm =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_spec2[eval_spec2$algorithm =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_spec2[eval_spec2$algorithm =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + N.inv.rarefied + (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  plot(allEffects(m_bart.alone))  
  
  
#############################################################################################################################################
###5#########################################################influence of NTI################################################################
#############################################################################################################################################   
#AUC
  AUCratio_nti <- merge(AUCratio,nti,by="species") 
  head(AUCratio_nti)
  #full model
  m_AUCratio_nti <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded  + NTI
                      + mop:algorithm
                      + extrapol:algorithm
                      + algorithm:NTI
                      + (1|species),data=AUCratio_nti)
  anova(m_AUCratio_nti)
  plot(allEffects(m_AUCratio_nti))
  
  #fne alone
  fne.alone <- AUCratio_nti[AUCratio_nti$algorithm =="AUCratio_fne",]
  m_fne.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=fne.alone)
  anova(m_fne.alone)
  
  #glm alone
  glm.alone <- AUCratio_nti[AUCratio_nti$algorithm =="AUCratio_glm",]
  m_glm.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=glm.alone)
  anova(m_glm.alone)
  
  #bart alone
  bart.alone <- AUCratio_nti[AUCratio_nti$algorithm =="AUCratio_bart",]
  m_bart.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=bart.alone)
  anova(m_bart.alone)   
  
#SENSITIVITY
eval_nti <- merge(sens,nti,by="species")  
  #full model
  m_eval_nti <- glmmTMB(sens ~ extrapol +  threshold + algorithm + mop + invaded  + NTI
                       + mop:algorithm
                       + extrapol:algorithm
                       + algorithm:NTI
                       + (1|species),data=eval_nti,family=beta_family(link="logit"))
  car::Anova(m_eval_nti)
  plot(allEffects(m_eval_nti))

  #fne alone
  fne.alone <- eval_nti[eval_nti$algorithm =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  
  #glm alone
  glm.alone <- eval_nti[eval_nti$algorithm =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  
  #bart alone
  bart.alone <- eval_nti[eval_nti$algorithm =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + NTI +(1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  
#SPECIFICITY
eval_nti <- merge(spec,nti,by="species")  
  #full model
  m_eval_nti <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded  + NTI
                        + mop:algorithm
                        + extrapol:algorithm
                        + algorithm:NTI
                        + (1|species),data=eval_nti,family=beta_family(link="logit"))
  car::Anova(m_eval_nti)
  plot(allEffects(m_eval_nti))
  
  #fne alone
  fne.alone <- eval_nti[eval_nti$algorithm =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  
  #glm alone
  glm.alone <- eval_nti[eval_nti$algorithm =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + NTI + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_nti[eval_nti$algorithm =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + NTI +(1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#############################################################################################################################################
###6#########################################################influence of NICHE EXPANSION####################################################
#############################################################################################################################################

#AUC  
AUCratio_ne <- merge(AUCratio,nd,by="species")
  head(AUCratio_ne)
  #full model
  m_AUCratio_ne <- lmer(AUC ~ extrapol +  threshold + algorithm + mop + invaded  + niche.expansion
                         + mop:algorithm
                         + extrapol:algorithm
                         + algorithm:niche.expansion
                         + (1|species),data=AUCratio_ne)
  anova(m_AUCratio_ne)
  plot(allEffects(m_AUCratio_ne))
  
  #fne alone
  fne.alone <- AUCratio_ne[AUCratio_ne$algorithm =="AUCratio_fne",]
  m_fne.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=fne.alone)
  anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- AUCratio_ne[AUCratio_ne$algorithm =="AUCratio_glm",]
  m_glm.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + niche.expansion +  (1|species),data=glm.alone)
  anova(m_glm.alone)
  summary(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- AUCratio_ne[AUCratio_ne$algorithm =="AUCratio_bart",]
  m_bart.alone <- lmer(AUC ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=bart.alone)
  anova(m_bart.alone)  
  
#SENSITIVITY
eval_ne <- merge(sens,nd,by="species")  
  #full model
  m_eval_ne <- glmmTMB(sens ~ extrapol +  threshold + algorithm + mop + invaded  + niche.expansion
                       + mop:algorithm
                       + extrapol:algorithm
                       + algorithm:niche.expansion
                       + (1|species),data=eval_ne,family=beta_family(link="logit"))
  car::Anova(m_eval_ne)
  plot(allEffects(m_eval_ne))
  
  #fne alone
  fne.alone <- eval_ne[eval_ne$algorithm =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_ne[eval_ne$algorithm =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  summary(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_ne[eval_ne$algorithm =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ extrapol +  threshold + mop + invaded  + niche.expansion +(1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  summary(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
#SPECIFICITY
  eval_ne <- merge(spec,nd,by="species")  
  #full model
  m_eval_ne <- glmmTMB(spec ~ extrapol +  threshold + algorithm + mop + invaded  + niche.expansion
                       + mop:algorithm
                       + extrapol:algorithm
                       + algorithm:niche.expansion
                       + (1|species),data=eval_ne,family=beta_family(link="logit"))
  car::Anova(m_eval_ne)
  plot(allEffects(m_eval_ne))
  
  #fne alone
  fne.alone <- eval_ne[eval_ne$algorithm =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=fne.alone,family=beta_family(link="logit"))
  summary(m_fne.alone)
  car::Anova(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_ne[eval_ne$algorithm =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + niche.expansion + (1|species),data=glm.alone,family=beta_family(link="logit"))
  summary(m_glm.alone)
  car::Anova(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_ne[eval_ne$algorithm =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ extrapol +  threshold + mop + invaded  + niche.expansion +(1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  summary(m_bart.alone)
  plot(allEffects(m_bart.alone))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  