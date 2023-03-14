# set working directory to this script's location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#libraries
library(reshape2)
library(lme4)
library(lmerTest)
library(effects)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(ggsignif)

###########################################################################################################################################################################
###READ IN MODELS WITH CLIMATE ONLY################################################################################################################
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
  eval_data_c <- subset(eval_data_c, select = -c(AUCratio_glm,AUCratio_bart,AUCratio_fne))  #remove AUC metrids
  eval_data_c <- eval_data_c[eval_data_c$mop =="no",]
  
#sens <- melt(eval_data[c(1,5:7,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))
  #spec <- melt(eval_data[c(1,8:10,11:16)],id.vars=c("species","BIAS","preds","extrapol","threshold","mop","invaded"))

###########################################################################################################################################################################
###READ IN NICHEMAPPER MODEL RESULS########################################################################################################################################
###########################################################################################################################################################################
#read in NicheMapper eval statististic###
  f_nm_eval <- read.csv(paste0("../4_NicheMapper/", "NicheMapperStats.csv"))
  head(f_nm_eval) 
#some housekeeping to align terminology
  f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="europe","eur_full")  
  f_nm_eval$background<-replace(f_nm_eval$background, f_nm_eval$background =="dispersal","inv_full")
  
  f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="species","NM_species")
  f_nm_eval$level<-replace(f_nm_eval$level, f_nm_eval$level =="intraspecific","NM_intra")
  
  nm_temp <- data.frame(f_nm_eval$background,f_nm_eval$species,f_nm_eval$level,f_nm_eval$SENS,f_nm_eval$SPEC)
  colnames(nm_temp) <- c("bg","species","sdm_method","sens","spec")
  head(nm_temp)
  nm_temp$sens <- replace(nm_temp$sens, nm_temp$sens ==0,0.001) 
  nm_temp$sens <- replace(nm_temp$sens, nm_temp$sens ==1,0.999)
  nm_temp$spec <- replace(nm_temp$spec, nm_temp$spec ==0,0.001) 
  nm_temp$spec <- replace(nm_temp$spec, nm_temp$spec ==1,0.999)
  
  nm_temp$invaded <- nm_temp$bg
  nm_temp <- nm_temp[c(2:6)]
  nm_temp$invaded[nm_temp$invaded == 'eur_full'] <- 'whole_europe'
  nm_temp$invaded[nm_temp$invaded == 'inv_full'] <- 'dispersal'
  colnames(nm_temp)[2] <- "algorithm"
  
###########################################################################################################################################################################
###cOMPARE MECHANISTIC AND CORRELATIVE MODELS#############################################################################################################################
#######################################################################################################################################################################
 
###   
#COMP 1: whole-europe background for NM and SDM, extrapolation only SDM, 2.5E SDM
###
  sdm_comp1 <- eval_data_c[eval_data_c$invaded =="whole_europe",]
  sdm_comp1 <- sdm_comp1[sdm_comp1$extrapol =="yes",]
  sdm_comp1 <- sdm_comp1[sdm_comp1$threshold ==5,]
  sdm_comp1 <- sdm_comp1[c(1:7,13)]
    t1 <- melt(sdm_comp1[c(1:4,8)],id=c('species','invaded'))
      colnames(t1)[4] <- "sens"
    t2 <- melt(sdm_comp1[c(1,5:8)],id=c('species','invaded'))
      colnames(t2)[4] <- "spec"
    t1$spec <- t2$spec
    colnames(t1)[3] <- 'algorithm'
    sdm_comp1 <- t1[c(1,3:5,2)]
    sdm_comp1$algorithm <- t2$variable
    head(sdm_comp1)
    sdm_comp1$sens <- replace(sdm_comp1$sens, sdm_comp1$sens ==0,0.001) 
    sdm_comp1$sens <- replace(sdm_comp1$sens, sdm_comp1$sens ==1,0.999)
    sdm_comp1$spec <- replace(sdm_comp1$spec, sdm_comp1$spec ==0,0.001) 
    sdm_comp1$spec <- replace(sdm_comp1$spec, sdm_comp1$spec ==1,0.999)
  nm_comp1 <- nm_temp[nm_temp$invaded =="whole_europe",]
  head(nm_comp1)
  data_comp1 <- rbind(nm_comp1,sdm_comp1)

  #model
  m_spec <- glmmTMB(spec ~ algorithm + (1|species),data_comp1, family=beta_family(link="logit"))
  m_spec <- lmer(spec ~ algorithm + (1|species),data_comp1)
    car::Anova(m_spec)
    pairs(emmeans(m_spec,  ~ algorithm))
    
    plot(allEffects(m_spec))
  
    #manual jittering of upper and lower data points
    upper <- data_comp1[data_comp1$spec >0.99,]
    upper.jitter <- runif(nrow(upper),min=0,max=0.05)
    upper$spec <- upper$spec - upper.jitter
    edit.upper <- 1-max(upper$spec)
    upper$spec <- upper$spec+edit.upper
    lower <- data_comp1[data_comp1$spec <0.01,]
    lower.jitter <- runif(nrow(lower),min=0,max=0.05)
    lower$spec <- lower$spec + lower.jitter
    edit.lower <- min(lower$spec)
    lower$spec <- lower$spec-edit.lower
    other <- data_comp1[!data_comp1$spec <0.01,]
    other <- other[!other$spec >0.99,]
    data_comp1 <- rbind(lower, other,upper)
    nrow(data_comp1)
    
data_comp1$modelclass <- ifelse(grepl("spec_",data_comp1$algorithm),'corr','mech')
  spec_cloud_plot <- ggplot(data_comp1, aes(x = algorithm, y = spec,fill=modelclass)) + scale_fill_manual(values=c("#E0E0E0", "#A0A0A0"))+
      geom_boxplot(fatten=NULL,width = .12, outlier.color = NA,lwd=0.2)+
      ggdist::stat_halfeye(adjust = .5, width =1.1, justification = -0.1, .width = 0, point_colour = NA) + 
      ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01,color='black') + 
      coord_cartesian(xlim = c(1.2, NA),ylim=c(0,1.1))+theme_classic()+ stat_summary(fun=mean, geom="point", shape=20, size=1.5, color="black", fill="black")
    spec_cloud_plot <- spec_cloud_plot + 
      scale_x_discrete(limit = c("spec_glm", "spec_bart", "spec_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NM(sp)","NM(intra)"))
    spec_cloud_plot <- spec_cloud_plot + xlab("") + ylab("absences succes rate")
    spec_cloud_plot <- spec_cloud_plot + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
    spec_cloud_plot <- spec_cloud_plot + theme(axis.text.x = element_text(size = 7))
    spec_cloud_plot <- spec_cloud_plot + theme(legend.position = "none") 
    spec_cloud_plot 
    spec_cloud_plot2 <- spec_cloud_plot + geom_signif(comparisons = list(c("spec_fne", "NM_intra")),y_position = 1.025,annotation = "",tip_length = 0.005,size=0.1)
    spec_cloud_plot2 <- spec_cloud_plot2 + geom_signif(comparisons = list(c("spec_bart", "NM_species")),y_position = 1.05,annotation = "",tip_length = 0.005,size=0.1)
    spec_cloud_plot2 <- spec_cloud_plot2 +  scale_y_continuous(breaks=seq(0,1,0.2))
    spec_cloud_plot2 <- spec_cloud_plot2 + annotate("text", label = "B", x = 0.78, y = 1.1,size=3)
    spec_cloud_plot2
    ggsave(filename="./spec_cloud_plot.png",
           plot=spec_cloud_plot2,
           device='png',
           width=89,
           height=69,
           units="mm")
  
###   
#COMP 2: whole-europe background for NM and SDM, extrapolation only SDM, 2.5E SDM
###
sdm_comp2 <- eval_data_c[eval_data_c$invaded =="whole_europe",]
  sdm_comp2 <- sdm_comp2[sdm_comp2$extrapol =="yes",]
  sdm_comp2 <- sdm_comp2[sdm_comp2$threshold ==2.5,]
  sdm_comp2 <- sdm_comp2[c(1:7,13)]
  t1 <- melt(sdm_comp2[c(1:4,8)],id=c('species','invaded'))
  colnames(t1)[4] <- "sens"
  t2 <- melt(sdm_comp2[c(1,5:8)],id=c('species','invaded'))
  colnames(t2)[4] <- "spec"
  t1$spec <- t2$spec
  colnames(t1)[3] <- 'algorithm'
  sdm_comp2 <- t1[c(1,3:5,2)]
  head(sdm_comp2)
  sdm_comp2$sens <- replace(sdm_comp2$sens, sdm_comp2$sens ==0,0.001) 
  sdm_comp2$sens <- replace(sdm_comp2$sens, sdm_comp2$sens ==1,0.999)
  sdm_comp2$spec <- replace(sdm_comp2$spec, sdm_comp2$spec ==0,0.001) 
  sdm_comp2$spec <- replace(sdm_comp2$spec, sdm_comp2$spec ==1,0.999)
  nm_comp2 <- nm_temp[nm_temp$invaded =="whole_europe",]
  head(nm_comp2)
  data_comp2 <- rbind(nm_comp2,sdm_comp2)
  
  #model
  m_spec <- glmmTMB(spec ~ algorithm + (1|species),data_comp2, family=beta_family(link="logit"))
  m_spec <- lmer(spec ~ algorithm + (1|species),data_comp2)
  car::Anova(m_spec)
  pairs(emmeans(m_spec,  ~ algorithm))
  
  plot(allEffects(m_spec))
  
  sens_cloud_plot <- ggplot(data_comp2, aes(x = algorithm, y = sens)) + ggdist::stat_halfeye(adjust = .5, width =0.8, justification = -.2, .width = 0, point_colour = NA) + 
    geom_boxplot(fatten=NULL,width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
    coord_cartesian(xlim = c(1.2, NA))+theme_bw() + stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")
  sens_cloud_plot <- sens_cloud_plot + scale_x_discrete(limit = c("sens_glm", "sens_bart", "sens_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
  sens_cloud_plot <- sens_cloud_plot + xlab("") + ylab("SENSITIVITY")
  sens_cloud_plot <- sens_cloud_plot + theme(text = element_text(size=15)) 
  sens_cloud_plot
  
###   
#COMP 3: 'dispersal since introduction' background for NM and SDM, extrapolation only SDM, 2.5E SDM
###
sdm_comp3 <- eval_data_c[eval_data_c$invaded =="dispersal",]
  sdm_comp3 <- sdm_comp3[sdm_comp3$extrapol =="yes",]
  sdm_comp3 <- sdm_comp3[sdm_comp3$threshold ==5,]
  sdm_comp3 <- sdm_comp3[c(1:7,13)]
  t1 <- melt(sdm_comp3[c(1:4,8)],id=c('species','invaded'))
  colnames(t1)[4] <- "sens"
  t2 <- melt(sdm_comp3[c(1,5:8)],id=c('species','invaded'))
  colnames(t2)[4] <- "spec"
  t1$spec <- t2$spec
  colnames(t1)[3] <- 'algorithm'
  sdm_comp3 <- t1[c(1,3:5,2)]
  head(sdm_comp3)
  sdm_comp3$sens <- replace(sdm_comp3$sens, sdm_comp3$sens ==0,0.001) 
  sdm_comp3$sens <- replace(sdm_comp3$sens, sdm_comp3$sens ==1,0.999)
  sdm_comp3$spec <- replace(sdm_comp3$spec, sdm_comp3$spec ==0,0.001) 
  sdm_comp3$spec <- replace(sdm_comp3$spec, sdm_comp3$spec ==1,0.999)
  nm_comp3 <- nm_temp[nm_temp$invaded =="whole_europe",]
  head(nm_comp3)
  data_comp3 <- rbind(nm_comp3,sdm_comp3)
  
  #model
  m_spec <- glmmTMB(spec ~ algorithm + (1|species),data_comp3, family=beta_family(link="logit"))
  m_spec <- lmer(spec ~ algorithm + (1|species),data_comp3)
  car::Anova(m_spec)
  pairs(emmeans(m_spec,  ~ algorithm))
  
  plot(allEffects(m_spec))
  
  sens_cloud_plot <- ggplot(data_comp3, aes(x = algorithm, y = sens)) + ggdist::stat_halfeye(adjust = .5, width =0.8, justification = -.2, .width = 0, point_colour = NA) + 
    geom_boxplot(fatten=NULL,width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
    coord_cartesian(xlim = c(1.2, NA))+theme_bw() + stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")
  sens_cloud_plot <- sens_cloud_plot + scale_x_discrete(limit = c("sens_glm", "sens_bart", "sens_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
  sens_cloud_plot <- sens_cloud_plot + xlab("") + ylab("SENSITIVITY")
  sens_cloud_plot <- sens_cloud_plot + theme(text = element_text(size=15)) 
  sens_cloud_plot
  
###   
#COMP 4: 'dispersal since introduction' background for NM and SDM, extrapolation only SDM, 2.5E SDM
###
sdm_comp4 <- eval_data_c[eval_data_c$invaded =="dispersal",]
  sdm_comp4 <- sdm_comp4[sdm_comp4$extrapol =="yes",]
  sdm_comp4 <- sdm_comp4[sdm_comp4$threshold ==2.5,]
  sdm_comp4 <- sdm_comp4[c(1:7,13)]
  t1 <- melt(sdm_comp4[c(1:4,8)],id=c('species','invaded'))
  colnames(t1)[4] <- "sens"
  t2 <- melt(sdm_comp4[c(1,5:8)],id=c('species','invaded'))
  colnames(t2)[4] <- "spec"
  t1$spec <- t2$spec
  colnames(t1)[3] <- 'algorithm'
  sdm_comp4 <- t1[c(1,3:5,2)]
  head(sdm_comp4)
  sdm_comp4$sens <- replace(sdm_comp4$sens, sdm_comp4$sens ==0,0.001) 
  sdm_comp4$sens <- replace(sdm_comp4$sens, sdm_comp4$sens ==1,0.999)
  sdm_comp4$spec <- replace(sdm_comp4$spec, sdm_comp4$spec ==0,0.001) 
  sdm_comp4$spec <- replace(sdm_comp4$spec, sdm_comp4$spec ==1,0.999)
  nm_comp4 <- nm_temp[nm_temp$invaded =="whole_europe",]
  head(nm_comp4)
  data_comp4 <- rbind(nm_comp4,sdm_comp4)
  
  #model
  m_spec <- glmmTMB(spec ~ algorithm + (1|species),data_comp4, family=beta_family(link="logit"))
  m_spec <- lmer(spec ~ algorithm + (1|species),data_comp4)
  car::Anova(m_spec)
  pairs(emmeans(m_spec,  ~ algorithm))
  
  plot(allEffects(m_spec))
  
  sens_cloud_plot <- ggplot(data_comp4, aes(x = algorithm, y = sens)) + ggdist::stat_halfeye(adjust = .5, width =0.8, justification = -.2, .width = 0, point_colour = NA) + 
    geom_boxplot(fatten=NULL,width = .12, outlier.color = NA) + ggdist::stat_dots(side = "left", justification = 1.1,binwidth = .01) + 
    coord_cartesian(xlim = c(1.2, NA))+theme_bw() + stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")
  sens_cloud_plot <- sens_cloud_plot + scale_x_discrete(limit = c("sens_glm", "sens_bart", "sens_fne","NM_species","NM_intra"),labels = c("GLM","BART","FNE","NicheMapper(species)","NicheMapper(intra)"))
  sens_cloud_plot <- sens_cloud_plot + xlab("") + ylab("SENSITIVITY")
  sens_cloud_plot <- sens_cloud_plot + theme(text = element_text(size=15)) 
  sens_cloud_plot