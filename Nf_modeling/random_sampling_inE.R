# Function: rs.inE
# Laura Jimenez
# Last version: June 2021
# Sampling points for a geographical area

# Description: -----------------------------------------------------------
# The function "rs.inE" takes a random sample of environmental combinations 
# inside the study area. The environmental data from a rasterstack is clipped by
# a shapefile that delimits the area before the sample is taken.

## Parameters:
# region = a shapefile of the study area (polygon)
# N = the sample size
# Estck = a rasterstack that contains at least two layers with environmental data

## Output:
# A matrix with two or more columns of samples that contain environmental data

# Function's code: -------------------------------------------------------------

# Get a random sample of points inside the given polygon
# and extract their environmental values
rs.inE <- function(region,N,Estck){
  # crop and mask the environmental layers with the polygon
  clip <- mask(crop(Estck,region),region)
  # get rid of cells with NA values = indices
  ind <- which(!is.na(clip[[1]][]))
  # get a random sample of indices
  sam <- sample(ind,N,replace = T)
  # choose the points corresponding to the selected indices
  pnts <- clip[][sam,]
  return(pnts)
}


## read libraries --------------
library(raster)
library(rgdal)


# End