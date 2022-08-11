# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# LOAD PACKAGES ####
library(lme4)
library(reshape2)
library(lmerTest)
library(sjPlot)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(effects)
library(car)

###################################################################################################################################
###READ IN DATA####################################################################################################################
###################################################################################################################################

#read in sample size data
ss <- read.csv(paste0("../../Modeling_bird_invasions/Nf_modeling/", "allspecies_samplesizes_v3.csv"))
colnames(ss)[1] <- "species"

#READ IN 2.5 data

#read in SDM model evaluation statistics
eur_full <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/2.5E/DS_eval_metrics/", "eval_metrics_eur_full.csv"))
  eur_full$bg <- rep("eur_full",nrow(eur_full))
eur_mess <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/2.5E/DS_eval_metrics/", "eval_metrics_eur_mess.csv"))
  eur_mess$bg <- rep("eur_mess",nrow(eur_mess))
inv_full <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/2.5E/DS_eval_metrics/", "eval_metrics_inv_full.csv"))
  inv_full$bg <- rep("inv_full",nrow(inv_full))
inv_mess <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/2.5E/DS_eval_metrics/", "eval_metrics_inv_mess.csv"))
  inv_mess$bg <- rep("inv_mess",nrow(inv_mess))
  
eval_data <- rbind(eur_full,eur_mess,inv_full, inv_mess)
  #head(eval_data)

  eval_data_2.5 <- eval_data
  eval_data_2.5$E <- rep(2.5,nrow(eval_data_2.5))

#READ IN 5 data

#read in SDM model evaluation statistics
eur_full <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_eur_full.csv"))
  eur_full$bg <- rep("eur_full",nrow(eur_full))
eur_mess <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_eur_mess.csv"))
  eur_mess$bg <- rep("eur_mess",nrow(eur_mess))
inv_full <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_inv_full.csv"))
  inv_full$bg <- rep("inv_full",nrow(inv_full))
inv_mess <- read.csv(paste0("../../Modeling_bird_invasions/PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_inv_mess.csv"))
  inv_mess$bg <- rep("inv_mess",nrow(inv_mess))

eval_data <- rbind(eur_full,eur_mess,inv_full, inv_mess)
  #head(eval_data)

  eval_data_5 <- eval_data
  eval_data_5$E <- rep(5,nrow(eval_data_5))

#all data
eval_all <- rbind(eval_data_2.5,eval_data_5)
  head(eval_all)
  eval_all$E <- as.factor(eval_all$E)

###################################################################################################################################
###AUC RATIO#######################################################################################################################
###################################################################################################################################  
eval_AUCratio <- melt(eval_all[c(1:4,17:18)],id.vars = c("bg","species","E"), variable.name = "sdm_method")
  head(eval_AUCratio)  
  colnames(eval_AUCratio)[5] <- "AUCratio"
  #full model
  compare <- lmer(AUCratio ~ sdm_method+bg+E+sdm_method:E + bg:E + (1|species),eval_AUCratio)
  anova(compare)
  #reduced model
  compare <- lmer(AUCratio ~ sdm_method+bg+E + (1|species),eval_AUCratio)
  anova(compare)

###################################################################################################################################
###SENSITIVITY#####################################################################################################################
###################################################################################################################################  
eval_sens <- melt(eval_all[c(1,11:13,17:18)],id.vars = c("bg","species","E"), variable.name = "sdm_method")
  head(eval_sens)
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==0,0.001) 
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==1,0.999) 
  colnames(eval_sens)[5] <- "sens"
  #full model
  m_sens <- glmmTMB(sens ~ sdm_method+bg+E+sdm_method:E+bg:E+(1|species), eval_sens, family=beta_family(link="logit"))
  car::Anova(m_sens)
  #reduced model
  m_sens <- glmmTMB(sens ~ sdm_method+bg+E+(1|species), eval_sens, family=beta_family(link="logit"))
  car::Anova(m_sens)
  
###################################################################################################################################
###SPECIFICITY#####################################################################################################################
################################################################################################################################### 
eval_spec <- melt(eval_all[c(1,14:16,17:18)],id.vars = c("bg","species","E"), variable.name = "sdm_method")
    head(eval_spec)    
    eval_spec$value <- replace(eval_spec$value, eval_spec$value ==0,0.001) 
    eval_spec$value <- replace(eval_spec$value, eval_spec$value ==1,0.999) 
    colnames(eval_spec)[5] <- "spec"
    #full modem
    m_spec <- glmmTMB(spec ~ sdm_method+bg+E+sdm_method:E+bg:E+(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec)
    #reduced model
    m_spec <- glmmTMB(spec ~ sdm_method+bg+E+(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec)