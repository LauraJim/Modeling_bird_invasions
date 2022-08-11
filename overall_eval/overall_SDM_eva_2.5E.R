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
library(dplyr)

dir.create("../overall_eval/2.5E") #directory to store data

###################################################################################################################################
###READ IN DATA####################################################################################################################
###################################################################################################################################

#read in sample size data
ss <- read.csv(paste0("../../Modeling_bird_invasions/Nf_modeling/", "allspecies_samplesizes_v3.csv"))
colnames(ss)[1] <- "species"

#read niche dynamics metrics
nd <- read.csv(paste0("../../Modeling_bird_invasions/broenni/", "niche_dynamics_EU_75perc.csv"))

#read niche truncation index
nti <- read.csv(paste0("../../Modeling_bird_invasions/NTI/", "nti.csv"))


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
# head(eval_data)

###################################################################################################################################
###AUC ratio (SDM only)############################################################################################################
###################################################################################################################################  
#
###1: select AUC ratio data only
#
eval_AUCratio <- melt(eval_data[c(1:4,ncol(eval_data))],id.vars = c("bg","species"), variable.name = "sdm_method")
  head(eval_AUCratio)
  summary(eval_AUCratio$value)

#
###2: influence of SDM algorithm and background
#  
  #values contained between (above) 0 and 2, so probably OK to use a simple linear model...
  #full model
  m_AUCratio <- lmer(value ~ sdm_method + bg + sdm_method:bg + (1|species),data=eval_AUCratio)
    anova(m_AUCratio)
    shapiro.test(residuals(m_AUCratio)) #not perfect, but W~0.95
  #reduced model
  m_AUCratio <- lmer(value ~ sdm_method + bg + (1|species),data=eval_AUCratio)
    shapiro.test(residuals(m_AUCratio)) #not perfect, but W~0.95 so close enough for me
    anova(m_AUCratio)
  m_AUCratio_plot<- plot_model(m_AUCratio, type = "eff",terms=c("sdm_method","bg"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
  m_AUCratio_plot
  ggsave(filename="./2.5E/m_AUCratio_plot.png",
         plot=m_AUCratio_plot,
         device='png',
         width=297,
         height=210,
         units="mm")
  emmeans(m_AUCratio, specs = pairwise ~ sdm_method)
  
#
###3: influence of invasive range sample size
#  
eval_AUCratio2 <- merge(eval_AUCratio,ss,by="species")
  #full model (=reduced model)
  m_AUCratio2 <- lmer(value ~ sdm_method + bg + N.inv.rarefied + sdm_method:N.inv.rarefied + (1|species),data=eval_AUCratio2)
  anova(m_AUCratio2)
    m_AUCratio_plot2<- plot_model(m_AUCratio2, type = "eff",terms=c("sdm_method","N.inv.rarefied"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
    m_AUCratio_plot2
    ggsave(filename="./2.5E/m_AUCratio_plot_invasive_sample_size.png",
         plot=m_AUCratio_plot2,
         device='png',
         width=297,
         height=210,
         units="mm")
    emmeans(m_AUCratio2, specs = pairwise ~ N.inv.rarefied|sdm_method)
  
  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_fne",]
  m_fne.alone <- lmer(value ~ N.inv.rarefied+bg+ (1|species),data=fne.alone)
  anova(m_fne.alone)
  summary(m_fne.alone)
  
  #glm alone
  glm.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_glm",]
  m_glm.alone <- lmer(value ~ N.inv.rarefied+bg+ (1|species),data=glm.alone)
  anova(m_glm.alone)
  
  #bart alone
  bart.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_bart",]
  m_bart.alone <- lmer(value ~ N.inv.rarefied+bg+ (1|species),data=bart.alone)
  anova(m_bart.alone)
  
#
##4: influence of native range sample size
#   
  #full model
  m_AUCratio2 <- lmer(value ~ sdm_method + bg + N.rarefied + sdm_method:N.rarefied + (1|species),data=eval_AUCratio2)  
    anova(m_AUCratio2)
  #reduced model  
  m_AUCratio2 <- lmer(value ~ sdm_method + bg + N.rarefied + (1|species),data=eval_AUCratio2)  
    anova(m_AUCratio2) 
    m_AUCratio_plot2<- plot_model(m_AUCratio2, type = "eff",terms=c("sdm_method","N.rarefied"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
    m_AUCratio_plot2
  
  #explore interaction: sample size for each SDM seperatly
  #fne alone
  fne.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_fne",]
  m_fne.alone <- lmer(value ~ N.rarefied+bg+ (1|species),data=fne.alone)
  anova(m_fne.alone)  
    
  #glm alone
  glm.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_glm",]
  m_glm.alone <- lmer(value ~ N.rarefied+bg+ (1|species),data=glm.alone)
  anova(m_glm.alone)
    
  #bart alone
  bart.alone <- eval_AUCratio2[eval_AUCratio2$sdm_method =="AUCratio_bart",]
  m_bart.alone <- lmer(value ~ N.rarefied+bg+ (1|species),data=bart.alone)
  anova(m_bart.alone)
    
#
##5: influence of niche dynamics
#     
eval_AUCratio3 <- merge(eval_AUCratio,nd,by="species")
  m_AUCratio3 <- lmer(value ~ sdm_method + niche.expansion + sdm_method:niche.expansion + bg + (1|species),data=eval_AUCratio3)
  anova(m_AUCratio3)  
  m_AUCratio_plot3<- plot_model(m_AUCratio3, type = "eff",terms=c("sdm_method","niche.expansion"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
    
  #fne alone
  fne.alone <- eval_AUCratio3[eval_AUCratio3$sdm_method =="AUCratio_fne",]
  m_fne.alone <- lmer(value ~ niche.expansion+bg+ (1|species),data=fne.alone)
  anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_AUCratio3[eval_AUCratio3$sdm_method =="AUCratio_glm",]
  m_glm.alone <- lmer(value ~ niche.expansion+bg+ (1|species),data=glm.alone)
  anova(m_glm.alone)
  summary(m_glm.alone)
  plot(allEffects(m_glm.alone))
  
  #bart alone
  bart.alone <- eval_AUCratio3[eval_AUCratio3$sdm_method =="AUCratio_bart",]
  m_bart.alone <- lmer(value ~ niche.expansion+bg+ (1|species),data=bart.alone)
  anova(m_bart.alone)
  
#
##6: influence of niche truncation index
#  
eval_AUCratio4 <- merge(eval_AUCratio,nti,by="species")  
  #full model
  m_AUCratio4 <- lmer(value ~ sdm_method + NTI + sdm_method:NTI + bg + (1|species),data=eval_AUCratio4)
  anova(m_AUCratio4)   
  plot_model(m_AUCratio4, type = "eff",terms=c("sdm_method","NTI"))+theme_sjplot2()+ylab("AUCratio")+ylim(1,2)
  
  #fne alone
  fne.alone <- eval_AUCratio4[eval_AUCratio4$sdm_method =="AUCratio_fne",]
  m_fne.alone <- lmer(value ~ NTI+bg+ (1|species),data=fne.alone)
  anova(m_fne.alone)

  #glm alone
  glm.alone <- eval_AUCratio4[eval_AUCratio4$sdm_method =="AUCratio_glm",]
  m_glm.alone <- lmer(value ~ NTI+bg+ (1|species),data=glm.alone)
  anova(m_glm.alone)

  #bart alone
  bart.alone <- eval_AUCratio4[eval_AUCratio4$sdm_method =="AUCratio_bart",]
  m_bart.alone <- lmer(value ~ NTI+bg+ (1|species),data=bart.alone)
  anova(m_bart.alone) 

###################################################################################################################################  
###SENSITIVITY#####################################################################################################################
###################################################################################################################################  
#
###1: select sensitivity data only
#
eval_sens <- melt(eval_data[c(1,11:13,ncol(eval_data))],id.vars = c("bg","species"), variable.name = "sdm_method")
  head(eval_sens)
  #values contained between 0 and 1 so beta regression
  #glmmTMB can accomodate beta regression with random effects
  #data fudge to avoid actual 0 and 1 results
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==0,0.001) 
  eval_sens$value <- replace(eval_sens$value, eval_sens$value ==1,0.999) 
  colnames(eval_sens)[4] <- "sens"

#
###2: influence of SDM algorithm and background
#   
  #full model
  m_sens <- glmmTMB(sens ~ sdm_method+bg+sdm_method:bg+(1|species), eval_sens, family=beta_family(link="logit"))
    car::Anova(m_sens)
  #reduced model
  m_sens <- glmmTMB(sens ~ sdm_method+bg+(1|species), eval_sens, family=beta_family(link="logit"))
    car::Anova(m_sens) 
  m_sens_plot<- plot_model(m_sens, type = "eff",terms=c("sdm_method","bg"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
    m_sens_plot
    ggsave(filename="./2.5E/m_sens_plot.png",
         plot=m_sens_plot,
         device='png',
         width=297,
         height=210,
         units="mm")
    emmeans(m_sens, specs = pairwise ~ bg|sdm_method)
  
#
###3: influence of invasive range sample size
#  
eval_sens2 <- merge(eval_sens,ss,by="species")
  m_eval_sens2 <- glmmTMB(sens ~ sdm_method + N.inv.rarefied + sdm_method:N.inv.rarefied + bg+ (1|species),eval_sens2, family=beta_family(link="logit"))
    car::Anova(m_eval_sens2)
  m_eval_sens2 <- glmmTMB(sens ~ sdm_method + N.inv.rarefied + bg+ (1|species),eval_sens2, family=beta_family(link="logit"))
    car::Anova(m_eval_sens2)
    m_eval_sens_plot2<- plot_model(m_eval_sens2, type = "eff",terms=c("sdm_method","N.inv.rarefied"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
    m_eval_sens_plot2
    ggsave(filename="./2.5E/m_sens_plot_sample_size.png",
         plot= m_eval_sens_plot2,
         device='png',
         width=297,
         height=210,
         units="mm")  
 
#
##4: influence of native range sample size
#   
  #full model
  m_eval_sens2 <- glmmTMB(sens ~ sdm_method + N.rarefied + sdm_method:N.rarefied + bg+ (1|species),eval_sens2, family=beta_family(link="logit"))
  car::Anova(m_eval_sens2)
  
  #fne alone
  fne.alone <- eval_sens2[eval_sens2$sdm_method =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ N.rarefied+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot_model(m_fne.alone, type = "eff",terms=c("N.rarefied"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
  
  
  #glm alone
  glm.alone <- eval_sens2[eval_sens2$sdm_method =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ N.rarefied+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)

  #bart alone
  bart.alone <- eval_sens2[eval_sens2$sdm_method =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ N.rarefied+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)  
  
#
##5: influence of niche dynamics
# 
eval_sens3 <- merge(eval_sens,nd,by="species")  
  #full model
  m_eval_sens3 <- glmmTMB(sens ~ sdm_method + niche.expansion + sdm_method:niche.expansion + bg+ (1|species),eval_sens3, family=beta_family(link="logit"))
  car::Anova(m_eval_sens3)
  
  #fne alone
  fne.alone <- eval_sens3[eval_sens3$sdm_method =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ niche.expansion+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  summary(m_fne.alone)
  plot(allEffects(m_fne.alone))
  
  #glm alone
  glm.alone <- eval_sens3[eval_sens3$sdm_method =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ niche.expansion+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  
  #bart alone
  bart.alone <- eval_sens3[eval_sens3$sdm_method =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ niche.expansion+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  summary(m_bart.alone)
  
#
##6: influence of niche truncation index
#  
eval_sens4 <- merge(eval_sens,nti,by="species")  
  #full model
  m_eval_sens4 <- glmmTMB(sens ~ sdm_method + NTI + sdm_method:NTI + bg+ (1|species),eval_sens4, family=beta_family(link="logit"))
  car::Anova(m_eval_sens4)  
  plot_model(m_eval_sens4, type = "eff",terms=c("sdm_method","NTI"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
  
  #fne alone
  fne.alone <- eval_sens4[eval_sens4$sdm_method =="sens_fne",]
  m_fne.alone <- glmmTMB(sens ~ NTI+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)

  #glm alone
  glm.alone <- eval_sens4[eval_sens4$sdm_method =="sens_glm",]
  m_glm.alone <- glmmTMB(sens ~ NTI+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  
  #bart alone
  bart.alone <- eval_sens4[eval_sens4$sdm_method =="sens_bart",]
  m_bart.alone <- glmmTMB(sens ~ NTI+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  
#
###7A: Mechanistic model sensitivity  
#
  #remove mess as not available/relevant for NM (and not different anyhow)
  eval_sens <- eval_sens[!eval_sens$bg =="eur_mess",]
  eval_sens <- eval_sens[!eval_sens$bg =="inv_mess",]
  unique(eval_sens$bg)
  head(eval_sens)
  
  #read in NicheMapper eval statististic###
  f_nm_eval <- read.csv(paste0("../../Modeling_bird_invasions/NicheMapper/", "NicheMapperStats.csv"))
  head(f_nm_eval) 
    #some housekeeping to align terminology
    f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="europe","eur_full")  
    f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="dispersal","inv_full")
  
    f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="species","NM_species")
    f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="intraspecific","NM_intra")
  
    nm_temp <- data.frame(f_nm_eval$background,f_nm_eval$species,f_nm_eval$level,f_nm_eval$SENS)
    colnames(nm_temp) <- c("bg","species","sdm_method","sens")
    head(nm_temp)
    nm_temp$sens <- replace(nm_temp$sens, nm_temp$sens ==0,0.001) 
    nm_temp$sens <- replace(nm_temp$sens, nm_temp$sens ==1,0.999) 
    
    eval_sens_all <- rbind(eval_sens,nm_temp)
    head(eval_sens_all)
    unique(eval_sens_all$sdm_method)
    unique(eval_sens_all$bg)
  
  #full model
  m_eval_sens_all <- glmmTMB(sens ~ sdm_method + bg + sdm_method:bg+ (1|species),eval_sens_all, family=beta_family(link="logit"))
    car::Anova(m_eval_sens_all) 
  m_eval_sens_all <- glmmTMB(sens ~ sdm_method + bg + (1|species),eval_sens_all, family=beta_family(link="logit"))
    car::Anova(m_eval_sens_all)
  m_eval_sens_all_f <- glmmTMB(sens ~ sdm_method + (1|species),eval_sens_all, family=beta_family(link="logit"))
    car::Anova(m_eval_sens_all_f)
    summary(m_eval_sens_all_f)
    emmeans(m_eval_sens_all_f, specs = pairwise ~ sdm_method)
    #plot
    m_eval_sens_all_plot<- plot_model(m_eval_sens_all, type = "eff",terms=c("sdm_method","bg"))+theme_sjplot2()+ylab("sensitivity")+ylim(0,1)
    m_eval_sens_all_plot
    #raw data summarized
    eval_sens_all %>% group_by(sdm_method) %>% summarise(mean = mean(sens),sd=sd(sens))
    #plot of raw data (select one background only as results are identical)  
    eval_sens_all <- eval_sens_all[eval_sens_all$bg =="eur_full",]
    sens_cloud_plot <- ggplot(eval_sens_all, aes(x = sdm_method, y = sens)) + ggdist::stat_halfeye(adjust = .5, width = 1, justification = -.2, .width = 0, point_colour = NA) + 
      geom_boxplot(width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
      coord_cartesian(xlim = c(1.2, NA))+theme_bw()
      sens_cloud_plot <- sens_cloud_plot + scale_x_discrete(limit = c("sens_glm", "sens_bart", "sens_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
      sens_cloud_plot <- sens_cloud_plot + xlab("") + ylab("SENSITIVITY")
      sens_cloud_plot <- sens_cloud_plot + theme(text = element_text(size=25)) 
      sens_cloud_plot
      ggsave(filename="./2.5E/sens_cloud_plot.png",
           plot=sens_cloud_plot,
           device='png',
           width=420,
           height=210,
           units="mm")
#
##7B: influence of invasive range sample size
#   
  nm_sens <- merge(nm_temp,ss,by="species")     
  #full model (select one background only as results are identical)
  nm_sens <- nm_sens[nm_sens$bg =="eur_full",]
  m_nm_sens <- glmmTMB(sens ~ sdm_method + N.inv.rarefied + sdm_method:N.inv.rarefied +(1|species),nm_sens, family=beta_family(link="logit"))
  car::Anova(m_nm_sens)
  #reduced model
  m_nm_sens <- glmmTMB(sens ~ sdm_method + N.inv.rarefied +(1|species),nm_sens, family=beta_family(link="logit"))
  car::Anova(m_nm_sens)

#
##7C: influence of niche dynamics
#   
  nm_sens <- merge(nm_temp,nd,by="species")                                  
  #full model (select one background only as results are identical)
  nm_sens <- nm_sens[nm_sens$bg =="eur_full",]
  m_nm_sens <- glmmTMB(sens ~ sdm_method + niche.expansion + sdm_method:niche.expansion +(1|species),nm_sens, family=beta_family(link="logit"))
  car::Anova(m_nm_sens)
  #reduced model
  m_nm_sens <- glmmTMB(sens ~ sdm_method + niche.expansion  +(1|species),nm_sens, family=beta_family(link="logit"))
  car::Anova(m_nm_sens)
                           
###################################################################################################################################    
###SPECIFICITY#####################################################################################################################
################################################################################################################################### 
#
###1: select specificity data only
# 
eval_spec <- melt(eval_data[c(1,14:16,ncol(eval_data))],id.vars = c("bg","species"), variable.name = "sdm_method")
  head(eval_spec)
  #values contained between 0 and 1 so beta regression
  #glmmTMB can accommodate beta regression with random effects
  #data fudge to avoid actual 0 and 1 results
  eval_spec$value <- replace(eval_spec$value, eval_spec$value ==0,0.001) 
  eval_spec$value <- replace(eval_spec$value, eval_spec$value ==1,0.999) 
  colnames(eval_spec)[4] <- "spec"

#
###2: influence of SDM algorithm and background
#  
  #full model
  m_spec <- glmmTMB(spec ~ sdm_method+bg+sdm_method:bg+(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec)
  #reduced model
  m_spec <- glmmTMB(spec ~ sdm_method+bg+(1|species), eval_spec, family=beta_family(link="logit"))
    car::Anova(m_spec) 
    m_spec_plot<- plot_model(m_spec, type = "eff",terms=c("sdm_method","bg"))+theme_sjplot2()+ylab("specitivity")+ylim(0,1)
    m_spec_plot
    ggsave(filename="./2.5E/m_spec_plot.png",
           plot=m_spec_plot,
           device='png',
           width=297,
           height=210,
           units="mm")
    emmeans(m_spec, specs = pairwise ~ bg|sdm_method)
    
#
###3: influence of invasive range sample size
#  
eval_spec2 <- merge(eval_spec,ss,by="species")
    m_eval_spec2 <- glmmTMB(spec ~ sdm_method + N.inv.rarefied + sdm_method:N.inv.rarefied + bg+  + (1|species),eval_spec2, family=beta_family(link="logit"))
      car::Anova(m_eval_spec2)
    m_eval_spec2 <- glmmTMB(spec ~ sdm_method + N.inv.rarefied + bg+  + (1|species),eval_spec2, family=beta_family(link="logit"))
      car::Anova(m_eval_spec2)
      m_eval_spec_plot2<- plot_model(m_eval_spec2, type = "eff",terms=c("sdm_method","N.inv.rarefied"))+theme_sjplot2()+ylab("specitivity")+ylim(0,1)
      m_eval_spec_plot2
    
      ggsave(filename="./2.5E/m_spec_plot_sample_size.png",
           plot= m_eval_spec_plot2,
           device='png',
           width=297,
           height=210,
           units="mm")  
    
#
##4: influence of native range sample size
#   
#full model
  m_eval_spec2 <- glmmTMB(spec ~ sdm_method + N.rarefied + sdm_method:N.rarefied + bg+ (1|species),eval_spec2, family=beta_family(link="logit"))
  car::Anova(m_eval_spec2)   
  
  #fne alone
  fne.alone <- eval_spec2[eval_spec2$sdm_method =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ N.inv.rarefied+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone) 

  #glm alone
  glm.alone <- eval_spec2[eval_spec2$sdm_method =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ N.inv.rarefied+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)

  #bart alone
  bart.alone <- eval_spec2[eval_spec2$sdm_method =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ N.inv.rarefied+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  
#
##5: influence of niche dynamics
# 
eval_spec3 <- merge(eval_spec,nd,by="species")
  #full model
  m_eval_spec3 <- glmmTMB(spec ~ sdm_method + niche.expansion + sdm_method:niche.expansion + bg+  + (1|species),eval_spec3, family=beta_family(link="logit"))
  car::Anova(m_eval_spec3)

  #fne alone
  fne.alone <- eval_spec3[eval_spec3$sdm_method =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ niche.expansion+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone) 
  summary(m_fne.alone)
  
  #glm alone
  glm.alone <- eval_spec3[eval_spec3$sdm_method =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ niche.expansion+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)
  summary(m_glm.alone)
  
  #bart alone
  bart.alone <- eval_spec3[eval_spec3$sdm_method =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ niche.expansion+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  summary(m_bart.alone)
  
#
##6: influence of niche truncation index
#   
eval_spec4 <- merge(eval_spec,nti,by="species") 
  #full model
  m_eval_spec4 <- glmmTMB(spec ~ sdm_method + NTI + sdm_method:NTI + bg+  + (1|species),eval_spec4, family=beta_family(link="logit"))
  car::Anova(m_eval_spec4)
  plot_model(m_eval_spec4, type = "eff",terms=c("sdm_method","NTI"))+theme_sjplot2()+ylab("specificity")+ylim(0,1)
  
  #fne alone
  fne.alone <- eval_spec4[eval_spec4$sdm_method =="spec_fne",]
  m_fne.alone <- glmmTMB(spec ~ NTI+bg+ (1|species),data=fne.alone,family=beta_family(link="logit"))
  car::Anova(m_fne.alone)
  
  #glm alone
  glm.alone <- eval_spec4[eval_spec4$sdm_method =="spec_glm",]
  m_glm.alone <- glmmTMB(spec ~ NTI+bg+ (1|species),data=glm.alone,family=beta_family(link="logit"))
  car::Anova(m_glm.alone)

  #bart alone
  bart.alone <- eval_spec4[eval_spec4$sdm_method =="spec_bart",]
  m_bart.alone <- glmmTMB(spec ~ NTI+bg+ (1|species),data=bart.alone,family=beta_family(link="logit"))
  car::Anova(m_bart.alone)
  
#
###7A: Mechanistic model specificity  
#  
  #remove mess as not available/relevant for NM (and not different anyhow)
  eval_spec <- eval_spec[!eval_spec$bg =="eur_mess",]
  eval_spec <- eval_spec[!eval_spec$bg =="inv_mess",]
  unique(eval_spec$bg)
  head(eval_spec)
  #read in NicheMapper eval statististic###
  f_nm_eval <- read.table("D:/Ststrubb/OneDrive - UGent/Projects/MC_paper/NatComm/Review/NicheMapperStats.txt",sep="\t",h=T)
  head(f_nm_eval) 
  #some housekeeping to align terminology
    f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="europe","eur_full")  
    f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="dispersal","inv_full")
    f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="species","NM_species")
    f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="intraspecific","NM_intra")
  
    nm_temp <- data.frame(f_nm_eval$background,f_nm_eval$species,f_nm_eval$level,f_nm_eval$SPEC)
    colnames(nm_temp) <- c("bg","species","sdm_method","spec")
    head(nm_temp)
    nm_temp$spec <- replace(nm_temp$spec, nm_temp$spec ==0,0.001) 
    nm_temp$spec <- replace(nm_temp$spec, nm_temp$spec ==1,0.999) 
    
    eval_spec_all <- rbind(eval_spec,nm_temp)
    head(eval_spec_all)
    unique(eval_spec_all$sdm_method)
    unique(eval_spec_all$bg)
    
  #full model
  m_eval_spec_all <- glmmTMB(spec ~ sdm_method + bg + sdm_method:bg+ (1|species),eval_spec_all, family=beta_family(link="logit"))
    car::Anova(m_eval_spec_all) 
  #reduced model    
  m_eval_spec_all <- glmmTMB(spec ~ sdm_method + bg + (1|species),eval_spec_all, family=beta_family(link="logit"))
    car::Anova(m_eval_spec_all)
    summary(m_eval_spec_all)
    emmeans(m_eval_spec_all, specs = pairwise ~ sdm_method|bg)
    #plot
    m_eval_spec_all_plot<- plot_model(m_eval_spec_all, type = "eff",terms=c("sdm_method","bg"))+theme_sjplot2()+ylab("specitivity")+ylim(0,1)
    m_eval_spec_all_plot
    #raw data summarized
    eval_spec_all %>% group_by(sdm_method) %>% summarize(mean=mean(spec),sd=sd(spec))
    #plots of raw data USING 'eur_full' background
    eval_spec_all_eur_full <- eval_spec_all[eval_spec_all$bg =="eur_full",]
    spec_cloud_plot <- ggplot(eval_spec_all_eur_full, aes(x = sdm_method, y = spec)) + ggdist::stat_halfeye(adjust = .5, width =0.8, justification = -.2, .width = 0, point_colour = NA) + 
      geom_boxplot(width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
      coord_cartesian(xlim = c(1.2, NA))+theme_bw()
      spec_cloud_plot <- spec_cloud_plot + scale_x_discrete(limit = c("spec_glm", "spec_bart", "spec_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
      spec_cloud_plot <- spec_cloud_plot + xlab("") + ylab("SPECITIVITY (EUR)")
      spec_cloud_plot <- spec_cloud_plot + theme(text = element_text(size=25)) 
      spec_cloud_plot
      ggsave(filename="./2.5E/spec_cloud_plot_eur_full.png",
           plot=spec_cloud_plot,
           device='png',
           width=420,
           height=210,
           units="mm")
    #plots of raw data USING 'inv_full' background    units="mm")
    eval_spec_all_inv_full <- eval_spec_all[eval_spec_all$bg =="inv_full",]
    spec_cloud_plot <- ggplot(eval_spec_all_inv_full, aes(x = sdm_method, y = spec)) + ggdist::stat_halfeye(adjust = .5, width =0.8, justification = -.2, .width = 0, point_colour = NA) + 
      geom_boxplot(width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
      coord_cartesian(xlim = c(1.2, NA))+theme_bw()
      spec_cloud_plot <- spec_cloud_plot + scale_x_discrete(limit = c("spec_glm", "spec_bart", "spec_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
      spec_cloud_plot <- spec_cloud_plot + xlab("") + ylab("SPECITIVITY (INV)")
      spec_cloud_plot <- spec_cloud_plot + theme(text = element_text(size=25)) 
      spec_cloud_plot
      ggsave(filename="./2.5E/spec_cloud_plot_inv_full.png",
           plot=spec_cloud_plot,
           device='png',
           width=420,
           height=210,
           units="mm")
#
##7B: influence of invasive range sample size
#   
  nm_spec <- merge(nm_temp,ss,by="species")     
  #full model
  m_nm_spec <- glmmTMB(spec ~ sdm_method + N.inv.rarefied + sdm_method:N.inv.rarefied + bg + (1|species),nm_spec, family=beta_family(link="logit"))
  car::Anova(m_nm_spec)
  #reduced model
  m_nm_spec <- glmmTMB(spec ~ sdm_method + N.inv.rarefied + bg + (1|species),nm_spec, family=beta_family(link="logit"))
  car::Anova(m_nm_spec)
  
#
##7C: influence of niche dynamics
#   
  nm_spec <- merge(nm_temp,nd,by="species")                                  
  #full model
  m_nm_spec <- glmmTMB(spec ~ sdm_method + niche.expansion + sdm_method:niche.expansion + bg + (1|species),nm_spec, family=beta_family(link="logit"))
  car::Anova(m_nm_spec)
  plot_model(m_nm_spec, type = "eff",terms=c("sdm_method","niche.expansion"))+theme_sjplot2()+ylab("specificity")+ylim(0,1)
      
  #species-level only
  nm_spec_species <- nm_spec[nm_spec$sdm_method=="NM_species",]    
  m_nm_spec_species <- glmmTMB(spec ~ niche.expansion + bg + (1|species),nm_spec_species, family=beta_family(link="logit"))
  car::Anova(m_nm_spec_species)
  
  #intraspecific-level only
  nm_spec_intra <- nm_spec[nm_spec$sdm_method=="NM_intra",]    
  m_nm_spec_intra <- glmmTMB(spec ~ niche.expansion + bg + (1|species),nm_spec_intra, family=beta_family(link="logit"))
  car::Anova(m_nm_spec_intra)
  
##
###NUMBER OF SUCCESFULL MODELS###############################################################################################################  
##  

#CRITERION SDM: AUCratio > 1, significant, varying omission rates  
#CRITERION NM: varying omission rates 

#omision for NichMapper              
f_nm_eval$omis_NM <- 1-f_nm_eval$SENS
  head(f_nm_eval) 

#CALCULATE NUMBER OF SUCCESFUL MODELS  
model_perf <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(model_perf) <- c("species","bg","succes","allowed_invas_omis","method")
  
  for.omission <- seq(0,0.2,0.01)
  #SDMs
  for (i in 1:length(for.omission)){
  allowed.invas.omission <- for.omission[i]  
    glm_OK <- ifelse(eval_data$AUCratio_glm > 1 & eval_data$pROCpval_glm <= 0.05 & eval_data$omis_glm <= allowed.invas.omission,1,0)
    bart_OK <- ifelse(eval_data$AUCratio_bart > 1 & eval_data$pROCpval_bart <= 0.05 & eval_data$omis_bart <= allowed.invas.omission,1,0) 
    fne_OK <- ifelse(eval_data$AUCratio_fne > 1 & eval_data$pROCpval_fne <= 0.05 & eval_data$omis_fne <= allowed.invas.omission,1,0)  
      t_glm <- data.frame(eval_data$species,eval_data$bg,glm_OK,rep(for.omission[i],nrow(eval_data)),rep("glm",nrow(eval_data)))
        colnames(t_glm) <- c("species","bg","succes","allowed_invas_omis","method")
      t_bart <- data.frame(eval_data$species,eval_data$bg,bart_OK,rep(for.omission[i],nrow(eval_data)),rep("bart",nrow(eval_data)))
        colnames(t_bart) <- c("species","bg","succes","allowed_invas_omis","method")
      t_fne <- data.frame(eval_data$species,eval_data$bg,fne_OK,rep(for.omission[i],nrow(eval_data)),rep("fne",nrow(eval_data)))
        colnames(t_fne) <- c("species","bg","succes","allowed_invas_omis","method")
    temp <- rbind(t_glm,t_bart,t_fne)
  model_perf <- rbind(model_perf,temp)
  }
  
  #NICHEMAPPER
  model_perf_NM <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(model_perf_NM) <- c("species","level","bg","succes","allowed_invas_omis","method")
      for (i in 1:length(for.omission)){
        allowed.invas.omission <- for.omission[i]  
        f_nm_eval$NM_OK <- ifelse(f_nm_eval$omis_NM <= allowed.invas.omission,1,0)
        t_NM <- data.frame(f_nm_eval$species,f_nm_eval$level, f_nm_eval$background,
                     f_nm_eval$NM_OK,rep(for.omission[i],nrow(f_nm_eval)),rep("NM",nrow(f_nm_eval)))
        colnames(t_NM) <- c("species","level","bg","succes","allowed_invas_omis","method")
    model_perf_NM <- rbind(model_perf_NM,t_NM)
    model_perf_NM 
  }
  

#PLOT1: EUROPE_FULL (eur_full)  
model_perf_eur_full <- model_perf[model_perf$bg=="eur_full",]
  unique(model_perf_eur_full$bg)
  model_perf_eur_full <- group_by(model_perf_eur_full, allowed_invas_omis, method)
  model_perf_eur_full <- summarise(model_perf_eur_full, species.succes = sum(succes,na.rm=TRUE))
model_perf_NM_eur_full <- model_perf_NM[model_perf_NM$bg=="eur_full",]
  unique(model_perf_NM_eur_full$bg)
  model_perf_NM_eur_full <- group_by(model_perf_NM_eur_full, allowed_invas_omis, level)
  model_perf_NM_eur_full <- summarise(model_perf_NM_eur_full, species.succes = sum(succes,na.rm=TRUE))
  colnames(model_perf_NM_eur_full)[2] <- "method" 
final_eur_full <- rbind(model_perf_eur_full,model_perf_NM_eur_full) 
  p_eur_full <- ggplot(final_eur_full,aes(x=allowed_invas_omis,y=species.succes,colour=method))+geom_point(size=2)+geom_line()+theme_sjplot2()+ylim(0,20)
  p_eur_full <- p_eur_full +geom_vline(xintercept = 0.05)
  p_eur_full <- p_eur_full + xlab("error rate in predicting invasive occurreces") + ylab("number of succesfull species")
  p_eur_full
  ggsave(filename="./2.5E/succesfull_models_2.5E.png",
         plot=  p_eur_full,
         device='png',
         width=297,
         height=210,
         units="mm")


##graphs look exactly the same as for all backgrounds
##AUCratio always above 1 and significant except for always myimon
##sensitivity=omission rate is not influenced by background
  
  
#PLOT2: INV_FULL (inv_full) 
#model_perf_inv_full <- model_perf[model_perf$bg=="inv_full",]
#  unique(model_perf_inv_full$bg)
#  model_perf_inv_full <- group_by(model_perf_inv_full, allowed_invas_omis, method)
#  model_perf_inv_full <- summarise(model_perf_inv_full, species.succes = sum(succes,na.rm=TRUE))
#model_perf_NM_inv_full <- model_perf_NM[model_perf_NM$bg=="europe",]
#  unique(model_perf_NM_inv_full$bg)
#  model_perf_NM_inv_full <- group_by(model_perf_NM_inv_full, allowed_invas_omis, level)
#  model_perf_NM_inv_full <- summarise(model_perf_NM_inv_full, species.succes = sum(succes,na.rm=TRUE))
#  colnames(model_perf_NM_inv_full)[2] <- "method" 
#final_inv_full <- rbind(model_perf_inv_full,model_perf_NM_inv_full) 
#  p_inv_full <- ggplot(final_inv_full,aes(x=allowed_invas_omis,y=species.succes,colour=method))+geom_point(size=2)+geom_line()+theme_sjplot2()
#  p_inv_full  