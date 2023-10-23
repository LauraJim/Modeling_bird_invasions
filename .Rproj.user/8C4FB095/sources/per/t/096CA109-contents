# Read functions
source(".\\Nf_modeling\\fit_wn_maha_model.R")

# Packages
library(scales)

# Read species IDs
sp.id <- read.csv("./Nf_modeling/allspecies_samplesizes_v3.csv",header=T)[,c(1,4,6)]

# Summary table with estimated parameters
mle.summary <- matrix(0,nrow = nrow(sp.id),ncol = 8)
colnames(mle.summary) <- c("speciesID", "N", "mu1", "mu2", "sigma11",
                           "sigma12", "sigma22", "def.pos")
# Confidence levels of ellipses in the figures
lvs <- c(0.25,0.5,0.75,0.95)

# Set colorpalette: 
colpal <- c(alpha("grey70",0.7), alpha("gold3",0.7), "purple3", "grey10",
            "brown")
# colpal[1] = random sample in accessible area (native range)
# colpal[2] = random sample in invaded region
# colpal[3] = fitted model for the fundamental niche
# colpal[4] = occurrences in accessible area
# colpal[5] = occurrences in invaded region


# BEFORE RUNNING: change file paths!!!
# See lines: 37,40,45,49,72,107,130,139
#c(12,15,18)
for (j in 1:nrow(sp.id)) {
  # 1) Sample sizes
  n.nat <- sp.id[j,2]
  n.inv <- sp.id[j,3]
  
  # 2) Read occurrence data
  # native range, used to fit the models
  f.occ.nat <- paste0("./Nf_modeling/occurrences-v3/",sp.id[j,1],"_native.csv")
  sp.occ.nat <- read.csv(f.occ.nat,header=T)[,-1]
  # invasive range, used to evaluate the models
  f.occ.inv <- paste0("./Nf_modeling/occurrences-v3/",sp.id[j,1],"_invasive.csv")
  sp.occ.inv <- read.csv(f.occ.inv,header=T)[,-1]
  
  # 3) Read tables with random samples in the native and invaded areas
  # native range, used to fit the models
  f.rs.nat <- paste0("./Nf_modeling/accessible-areas-v3/",sp.id[j,1],
                     "_native_range.csv")
  sp.rs.nat <- read.csv(f.rs.nat,header=T)[,-1]
  # invasive range, used to evaluate the models
  f.rs.inv <- paste0("./Nf_modeling/accessible-areas-v3/",sp.id[j,1],
                     "_invaded_region.csv")
  sp.rs.inv <- read.csv(f.rs.inv,header=T)[,-1]
  
  # 4) Apply functions to estimate parameters
  mle <- fitNiche(E.occ = sp.occ.nat[,3:4], E.samM = sp.rs.nat[,3:4])
  
  # 5)
  if(mle$dp == 1){ # Sigma is positive definite
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
    png(paste0("./Nf_modeling/Results-v3/",sp.id[j,1],"_modelfit.png"),
        width = 1800, height = 1800, res = 300, pointsize = 8)
    # x11()
    # background points from accessible area, using invaded regions to set plot limits
    plot(rbind(sp.rs.nat[,3:4],sp.rs.inv[,3:4]),col=colpal[1],pch=1,
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
  }
  else{
    # 5.1) save Mahalanobis model in the summary table
    mle.summary[j,] <- c(sp.id[j,1], n.nat, mle$maha.mu,
                         as.vector(mle$maha.sigma)[-2],mle$dp)
    # 5.2) Define ellipses using these estimates
    ellis <- list()
    for(i in 1:length(lvs)){
      ellis[[i]] <- ellipse::ellipse(x=mle$maha.sigma, centre=as.numeric(mle$maha.mu), level=lvs[i])
    }
    
    # Warning message
    print(paste("Warning!","Estimated matrix is not positive definite for species:",
                sp.id[j,1]))
    
    # 5.3) PLOT
    # plot will be saved as .png
    png(paste0("./Nf_modeling/Results-v3/",sp.id[j,1],"_modelfit.png"),
        width = 1800, height = 1800, res = 300, pointsize = 8)
    # x11()
    # background points from accessible area
    plot(rbind(sp.rs.nat[,3:4],sp.rs.inv[,3:4]),col=colpal[1],pch=1, xlab="PC1", ylab="PC2", cex.lab=1.3,
         main="Environmental space")
    # add points from invaded region
    points(sp.rs.inv[,3:4],col=colpal[2],pch=1) 
    # ellipse wn
    for(x in 1:length(lvs)){lines(ellis[[x]],col=colpal[3],lwd=2)} 
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
  }
  
  # 7) SAVE estimated parameters for all the species
  if(j==nrow(sp.id))
    write.csv(mle.summary,"./Nf_modeling/Results-v3/mle_allspecies_v3.csv",row.names = F)
  #print(mle.summary)
}


### END ####
# Laura Jimenez
# November, 2021
