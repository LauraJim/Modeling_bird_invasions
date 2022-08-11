# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# LOAD PACKAGES ####
library(raster)
library(reshape2)
library(lmerTest)
library(sjPlot)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(effects)
library(dplyr)
library(ggsignif)

###################################################################################################################################
###READ IN DATA####################################################################################################################
###################################################################################################################################
#EXTRAPOLATION

#read in SDM model evaluation statistics
eur_full <- read.csv(paste0("../PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_eur_full.csv"))
  eur_full$bg <- rep("eur_full",nrow(eur_full))
eur_mess <- read.csv(paste0("../PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_eur_mess.csv"))
  eur_mess$bg <- rep("eur_mess",nrow(eur_mess))
inv_full <- read.csv(paste0("../PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_inv_full.csv"))
  inv_full$bg <- rep("inv_full",nrow(inv_full))
inv_mess <- read.csv(paste0("../PA_modelling/DS/5E/DS_eval_metrics/", "eval_metrics_inv_mess.csv"))
  inv_mess$bg <- rep("inv_mess",nrow(inv_mess))
  
eval_data_extrapol <- rbind(eur_full,eur_mess,inv_full, inv_mess)
  
#CLAMPING

#read in SDM model evaluation statistics
eur_full <- read.csv(paste0("../PA_modelling/DS/5ECL/DS_eval_metrics/", "eval_metrics_eur_full.csv"))
eur_full$bg <- rep("eur_full",nrow(eur_full))
eur_mess <- read.csv(paste0("../PA_modelling/DS/5ECL/DS_eval_metrics/", "eval_metrics_eur_mess.csv"))
eur_mess$bg <- rep("eur_mess",nrow(eur_mess))
inv_full <- read.csv(paste0("../PA_modelling/DS/5ECL/DS_eval_metrics/", "eval_metrics_inv_full.csv"))
inv_full$bg <- rep("inv_full",nrow(inv_full))
inv_mess <- read.csv(paste0("../PA_modelling/DS/5ECL/DS_eval_metrics/", "eval_metrics_inv_mess.csv"))
inv_mess$bg <- rep("inv_mess",nrow(inv_mess))

eval_data_clamping <- rbind(eur_full,eur_mess,inv_full, inv_mess)  
  
  
###CHECK FOR DIFFERENCES
head(eval_data_extrapol)
head(eval_data_clamping)
  
eval_data_extrapol$extrapol_method <- rep("extra",nrow(eval_data_extrapol))  
eval_data_clamping$extrapol_method <- rep("clamping",nrow(eval_data_extrapol))

eval_data <- rbind(eval_data_clamping,eval_data_extrapol)

###################################################################################################################################
###AUC ratio (SDM only)############################################################################################################
###################################################################################################################################  
#
###1: select AUC ratio data only
#
eval_AUCratio <- melt(eval_data[c(1:4,17:18)],id.vars = c("bg","species","extrapol_method"), variable.name = "sdm_method")
  head(eval_AUCratio)
  summary(eval_AUCratio$value)
  mean(eval_AUCratio[eval_AUCratio$extrapol_method=="clamping",]$value,na.rm=TRUE)
  mean(eval_AUCratio[eval_AUCratio$extrapol_method=="extra",]$value,na.rm=TRUE)
#
###2: influence of SDM algorithm and background and clamping
#  
  #values contained between (above) 0 and 2, so probably OK to use a simple linear model...
  #full model
  m_AUCratio <- lmer(value ~ sdm_method + bg + extrapol_method+
                       sdm_method:bg + sdm_method:extrapol_method + bg:extrapol_method+
                       sdm_method:bg:extrapol_method+
                       sdm_method:bg + (1|species),data=eval_AUCratio)
    anova(m_AUCratio)
    shapiro.test(residuals(m_AUCratio)) 
    
  #reduced model
  m_AUCratio <- lmer(value ~ sdm_method + bg + extrapol_method+
                          + (1|species),data=eval_AUCratio)
  anova(m_AUCratio)

  m_AUCratio_plot<- plot_model(m_AUCratio, type = "eff",terms=c("sdm_method","bg","extrapol_method"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
  m_AUCratio_plot
  ggsave(filename="./5E/m_AUCratio_plot.png",
         plot=m_AUCratio_plot,
         device='png',
         width=297,
         height=210,
         units="mm")
  emmeans(m_AUCratio, specs = pairwise ~ extrapol_method)
  

###################################################################################################################################  
###SENSITIVITY#####################################################################################################################
###################################################################################################################################  
#
###1: select sensitivity data only
#
eval_sens <- melt(eval_data[c(1,11:13,17:18)],id.vars = c("bg","species","extrapol_method"), variable.name = "sdm_method")
  head(eval_sens)
  #values contained between 0 and 1 so beta regression
  #glmmTMB can accomodate beta regression with random effects
  #data fudge to avoid actual 0 and 1 results
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==0,0.001) 
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==1,0.999) 
  colnames(eval_sens)[5] <- "sens"
  
  mean(eval_sens[eval_sens$extrapol_method=="clamping",]$sens,na.rm=TRUE)
    sd(eval_sens[eval_sens$extrapol_method=="clamping",]$sens,na.rm=TRUE)
  mean(eval_sens[eval_sens$extrapol_method=="extra",]$sens,na.rm=TRUE)
    sd(eval_sens[eval_sens$extrapol_method=="extra",]$sens,na.rm=TRUE)

#
###2: influence of SDM algorithm and background
#   
  #full model
  m_sens <- glmmTMB(sens ~ sdm_method + bg + extrapol_method+
                      sdm_method:bg + sdm_method:extrapol_method + bg:extrapol_method+
                      sdm_method:bg:extrapol_method
                    +(1|species), eval_sens, family=beta_family(link="logit"))
    car::Anova(m_sens)
    
  #reduced model
  m_sens <- glmmTMB(sens ~ sdm_method + bg + extrapol_method+
                        +(1|species), eval_sens, family=beta_family(link="logit"))
    car::Anova(m_sens)
    
    m_sens_plot<- plot_model(m_sens, type = "eff",terms=c("sdm_method","bg","extrapol_method"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
    m_sens_plot
    ggsave(filename="./5E/m_sens_plot.png",
         plot=m_sens_plot,
         device='png',
         width=297,
         height=210,
         units="mm")
    emmeans(m_sens, specs = pairwise ~ bg|extrapol_method)
  
###################################################################################################################################    
###SPECIFICITY#####################################################################################################################
################################################################################################################################### 
#
###1: select specificity data only
# 
eval_spec <- melt(eval_data[c(1,14:16,17:18)],id.vars = c("bg","species","extrapol_method"), variable.name = "sdm_method")
  head(eval_spec)
  #values contained between 0 and 1 so beta regression
  #glmmTMB can accommodate beta regression with random effects
  #data fudge to avoid actual 0 and 1 results
  eval_spec$value <- replace(eval_spec$value, eval_spec$value ==0,0.001) 
  eval_spec$value <- replace(eval_spec$value, eval_spec$value ==1,0.999) 
  colnames(eval_spec)[5] <- "spec"
  
  mean(eval_spec[eval_spec$extrapol_method=="clamping",]$spec,na.rm=TRUE)
    sd(eval_spec[eval_spec$extrapol_method=="clamping",]$spec,na.rm=TRUE)
  mean(eval_spec[eval_spec$extrapol_method=="extra",]$spec,na.rm=TRUE)
    sd(eval_spec[eval_spec$extrapol_method=="extra",]$spec,na.rm=TRUE)

#
###2: influence of SDM algorithm and background
#  
  #full model
    m_spec <- glmmTMB(spec ~ sdm_method + bg + extrapol_method+
                        sdm_method:bg + sdm_method:extrapol_method + bg:extrapol_method+
                        sdm_method:bg:extrapol_method
                      +(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec)
  
    #reduced model
    m_spec <- glmmTMB(spec ~ sdm_method + bg + extrapol_method+
                                           +(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec) 
    m_spec_plot<- plot_model(m_spec, type = "eff",terms=c("sdm_method","bg","extrapol_method"))+theme_sjplot2()+ylab("specitivity")+ylim(0,1)
    m_spec_plot
    ggsave(filename="./5E/m_spec_plot.png",
           plot=m_spec_plot,
           device='png',
           width=297,
           height=210,
           units="mm")
    emmeans(m_spec, specs = pairwise ~ bg|extrapol_method)
    
