# Laura Jimenez and Carola Franzen
# First version: June 2020
# Last version: June 2021
# Project resulting ellipses back into G-space

# Description: -----------------------------------
# The function niche.G projects ellipses that define suitable environments for a 
# species on a map as potential niches. The regions in the geographical space
# are colored by different degrees of suitability.

## Parameters:
# Estck = a raster stack with more than two climatic layers
# mu = the mean of the columns that contain environmental data, such as 
#       temperature and precipitation 
# Sigma = the covariance of the environmental data linked with a species' 
#         occurrence

## Output:
# The function will produce a geographical map that represents areas that have 
# suitable environmental conditions for a species. Those potential ecological 
# niche regions are colored by different degrees of suitability. The map will
# automatically be saved as a tiff and a asci file.


# the function's code: niche.G --------------------------------------

niche.G <- function(Estck, mu, Sigma) {
  # calculate suitability for each cell
  sui.fun <- function(cell){
    X <- as.vector(cell, mode = "numeric")
    sui.ind <- exp(-mahalanobis(x= X, center= mu, cov= Sigma)/2)
    return(sui.ind)
  }
  # apply this function to the whole raster layer
  suit.rast <- calc(Estck,fun=sui.fun)
  
  return(suit.rast)
}

  
## read libraries -------------------

library(raster)


# End