# Functions: rs.inE, negloglike
# Laura Jimenez
# First version: February 2020
# Last version: June 2021

# Description: ----------------------------------------------------------
# Maximum likelihood approach of the fundamental niche estimation problem using 
# a weighted distribution where the weights represent the availability of 
# environmental combinations inside M. This approach contains three functions.

## Parameters rs.inE:
# region = a shapefile of the study area (polygon)
# N = the sample size
# Estck = a rasterstack that contains at least two layers with environmental data

## Parameters negloglike
# guess = a vector of length 5 when d=2, it contains the mu and A values as elements
# sam1 = matrix containing the original sample of environmental combinations that 
#       correspond to presences
# sam2 = matrix containing a second random sample of environmental combinations 
#       which come from the area of study (M)

# rasterstack, shapefile and occurrence points

# Calling packages


# FUNCTIONS -------------------------------------------------------------

# Negative log-likelihood function for theta=(mu,A)
## guess -- is a vector of length 5 when d=2, it contains the mu and A values as elements
## sam1 -- matrix containing the original sample of environmental combinations that correspond to presences
## sam2 -- matrix containing a second random sample of environmental combinations which come from the area of study (M)
negloglike <- function(guess,sam1,sam2){
  # define the parameters of interest from the guess parameter
  mu <- guess[1:2]
  A <- matrix(c(guess[3],guess[4],guess[4],guess[5]),nrow=2,ncol=2)
  # original sample size: number of presence points in the sample
  n <- nrow(sam1)
  # function that calculates quadratic terms, use inverse of matrix
  quad <- function(xi) { ax<-as.matrix(xi - mu); t(ax) %*% A %*% ax }
  q1 <- apply(sam1, 1, quad) # quadratic terms of presence points
  q2 <- apply(sam2, 1, quad) # quadratic terms of M points
  # negative log-likelihood value
  S <- 0.5*sum(q1) + n*log(sum(exp(-0.5*q2)))
  return(S)
}

# maximum likelihood
fitNiche <- function(E.occ, E.samM) {
  # calculate mu
  mu.ini <- colMeans(E.occ)
  # calculate A (covariance)
  Sig.ini <- cov(E.occ)
  # invert matrix sig.ini
  A.ini <- chol2inv(chol(Sig.ini))
  # whole vector of inicial values
  vals.ini <- c(mu.ini, A.ini[1,1], A.ini[1,2], A.ini[2,2])#c(mu.ini,A.ini[1,1],A.ini[1,2],A.ini[2,2])
  # fix the values of the samples used to evaluate the neg-log-likelihood
  like.fn <- function(theta){ negloglike(theta, sam1= E.occ, sam2= E.samM) } 
  find.mle <- optim(par=vals.ini, fn=like.fn, method="Nelder-Mead")
  mle <- find.mle$par
  mle.mu <- mle[1:2]
  mle.A <- matrix(c(mle[3:4],mle[4:5]),nrow=2,ncol=2)
  mle.Sig <- tryCatch(expr={chol2inv(chol(mle.A))}, error= function(e){NULL})

  # change column names for mle.Sig
  if(!is.null(mle.Sig)){
  colnames(mle.Sig) <- colnames(Sig.ini)
  rownames(mle.Sig) <- rownames(Sig.ini)
  }
    
  # wn = weighted normal distribution
  return(list(wn.mu = mle.mu, wn.sigma = mle.Sig, maha.mu = mu.ini, maha.sigma = Sig.ini))
}


## Read libraries -------------- 
# library(sp)
library(raster)
library(rgdal)
library(rgeos)

# End