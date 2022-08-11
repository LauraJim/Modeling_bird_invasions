library(lme4)
library(ggplot2)
library(ggpubr)


data.niche_th_sp <- read.table("E:/SDM/broenni/niche_results_th_sp.txt",sep="\t",h=T)
  head(data.niche_th_sp)
  
data.niche_th_env <- read.table("E:/SDM/broenni/niche_results_th_env.txt",sep="\t",h=T)
  head(data.niche_th_env)
  data.niche_th_env$equ[is.na(data.niche_th_env$equ)] <- 1
  data.niche_th_env$stability.index[is.na(data.niche_th_env$stability.index)] <- 0
  data.niche_th_env$expansion.index[is.na(data.niche_th_env$expansion.index)] <- 1
  data.niche_th_env$unfilling.index[is.na(data.niche_th_env$unfilling.index)] <- 1

###restrict to 0 -0.25range###  
  data.niche_th_env <- data.niche_th_env[!data.niche_th_env$threshold>0.25,]  

#############    
###Figures###
#############
  
#lines and points  
niche.dyn.exp <- ggplot(data.niche_th_env, aes(x=threshold, y=expansion.index, group=species)) +
    geom_line(aes(color=species))+
    geom_point(aes(color=species),size=3)+
    geom_point(aes(shape=species),size=3)+scale_shape_manual(values=seq(0,20)) +
    xlab("\n ") + ylab("niche expansion index") 
    niche.dyn.exp
  ggsave("niche.dyn.exp.png", plot=niche.dyn.exp, width=14, height=14, units="cm", dpi=1200)

niche.dyn.unf <- ggplot(data.niche_th_env, aes(x=threshold, y=unfilling.index, group=species)) +
    geom_line(aes(color=species))+
    geom_point(aes(color=species),size=3)+
    geom_point(aes(shape=species),size=3)+scale_shape_manual(values=seq(0,20))+
    xlab("\n ") + ylab("niche unfilling index") 
  niche.dyn.unf
  ggsave("niche.dyn.unf.png", plot=niche.dyn.unf, width=14, height=14, units="cm", dpi=1200)
  
niche.dyn.stab <- ggplot(data.niche_th_env, aes(x=threshold, y=stability.index, group=species)) +
    geom_line(aes(color=species))+
    geom_point(aes(color=species),size=3)+
    geom_point(aes(shape=species),size=3)+scale_shape_manual(values=seq(0,20))+
    xlab("\n ") + ylab("niche stability index")
  niche.dyn.stab
  ggsave("niche.dyn.stab.png", plot=niche.dyn.stab, width=14, height=14, units="cm", dpi=1200 )
  
#violin plots
niche.dyn.exp_vio <- ggplot(data.niche_th_env, aes(x=as.factor(threshold), y=expansion.index)) + geom_violin()  
  niche.dyn.exp_vio <- niche.dyn.exp_vio + 
                      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
                      stat_summary(fun.y=median, geom="point", size=5, color="red")+
                      xlab("\n ") + ylab("niche expansion index") 
  niche.dyn.exp_vio
  ggsave("niche.dyn.exp_vio.png", plot=niche.dyn.exp_vio, width=14, height=14, units="cm", dpi=1200)

niche.dyn.unf_vio <- ggplot(data.niche_th_env, aes(x=as.factor(threshold), y=unfilling.index)) + geom_violin()  
  niche.dyn.unf_vio <- niche.dyn.unf_vio + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
    stat_summary(fun.y=median, geom="point", size=5, color="red")+
    xlab("\n percentile climate density overlap between ranges") + ylab("niche unfilling index")
  niche.dyn.unf_vio
  ggsave("niche.dyn.unf_vio.png", plot=niche.dyn.unf_vio, width=14, height=14, units="cm", dpi=1200)
  
niche.dyn.stab_vio <- ggplot(data.niche_th_env, aes(x=as.factor(threshold), y=stability.index)) + geom_violin()  
  niche.dyn.stab_vio <- niche.dyn.stab_vio + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
    stat_summary(fun.y=median, geom="point", size=5, color="red") + 
    xlab("\n ") + ylab("niche stability index")
  niche.dyn.stab_vio
  ggsave("niche.dyn.stab_vio.png", plot=niche.dyn.stab_vio, width=14, height=14, units="cm", dpi=1200)
  
  
  
  niche.dyn.exp
  niche.dyn.unf
  niche.dyn.stab
  
  niche.dyn.exp_vio
  niche.dyn.unf_vio
  niche.dyn.stab_vio
  
  

  test <- ggarrange(niche.dyn.exp, NULL,niche.dyn.unf,NULL, niche.dyn.stab,niche.dyn.exp_vio, NULL,niche.dyn.unf_vio,NULL, niche.dyn.stab_vio, 
            labels = c("A","", "B","", "C","D","","E","","F"), widths=c(1,0.1,1,0.1,1,1,0.1,1,0.1,1),
            ncol = 5, nrow = 2,common.legend = TRUE, legend="right")
  
  setwd("E:/SDM/broenni")
  ggsave(filename="Fig_Revision_X1.png",
         plot=test,
         device='tiff',
         width=297,
         height=210,
         units="mm")
         
  
###############
###SUMMARIES###
###############
library(dplyr)
data.niche_th_env.75 <- data.niche_th_env[data.niche_th_env$threshold==0.25,]  
  data.niche_th_env.75
  
niche.summary <- data.niche_th_env.75 %>%summarise(niche_overlap_F.D.mean = mean(niche_overlap_F.D),
                                       niche_overlap_F.D.sd = sd(niche_overlap_F.D),
                                       niche_overlap_F.T.mean = mean(niche_overlap_T.D),
                                       niche_overlap_F.T.sd = sd(niche_overlap_T.D),
                                       stability.index.mean = mean(stability.index),
                                       stability.index.sd = sd(stability.index),
                                       expansion.index.mean = mean(expansion.index),
                                       expansion.index.sd = sd(expansion.index),
                                       unfilling.index.sd = sd(unfilling.index),
                                       unfilling.index.mean = mean(unfilling.index))
  niche.summary$sim1 <- sum(ifelse(data.niche_th_env.75$sim1>0.05,0,1))
  niche.summary$sim2 <- sum(ifelse(data.niche_th_env.75$sim2>0.05,0,1))
  niche.summary$equ <- sum(ifelse(data.niche_th_env.75$equ>0.05,0,1))
  niche.summary <- data.frame(t(niche.summary))
  colnames(niche.summary) <- "value"
  niche.summary
  write.table(niche.summary,"niche.summary.txt",sep="\t",row.names=TRUE)
  
niche.summary <- data.niche_th_env.90 %>%summarise(niche_overlap_F.D.median = median(niche_overlap_F.D),
                                                     niche_overlap_F.D.sd = sd(niche_overlap_F.D),
                                                     niche_overlap_F.T.median = median(niche_overlap_T.D),
                                                     niche_overlap_F.T.sd = sd(niche_overlap_T.D),
                                                     stability.index.median = median(stability.index),
                                                     stability.index.sd = sd(stability.index),
                                                     expansion.index.median = median(expansion.index),
                                                     expansion.index.sd = sd(expansion.index),
                                                     unfilling.index.sd = sd(unfilling.index),
                                                     unfilling.index.median = median(unfilling.index))
  niche.summary$sim1 <- sum(ifelse(data.niche_th_env.90$sim1>0.05,0,1))
  niche.summary$sim2 <- sum(ifelse(data.niche_th_env.90$sim2>0.05,0,1))
  niche.summary$equ <- sum(ifelse(data.niche_th_env.90$equ>0.05,0,1))
  niche.summary <- data.frame(t(niche.summary))
  colnames(niche.summary) <- "value"
  niche.summary 
                                  
                                  
  
  