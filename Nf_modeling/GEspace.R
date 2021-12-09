# Function "GE.space"
# Authors: Laura Jiménez and Carola Franzen
# 2021-05-25

# Description:"GE.space" ------------------
# This function creates a plot with two two graphs in ggplot. One is the 
# geographical space, the other is the environmental space of a species' 
# occurrence. 

## Parameters: 
# bckgrnd = contains random points within a limited geographical space that 
#           serve as a background for a map
# add.poly = a polygon of a geographical area in which a species occurs
# GE.occ = a matrix that contains the occurrence points of a species and the 
#       environmental data of the find location
# save.p = if set, the graph will be saved, if not, it will not be saved

## Output:
# A graph with two plots. The left plot shows a map with the geographical spread
# of the species, the right plot shows a the environmental space. There are two
# different geographical outputs: 1) with background points, 2) with a polygon.


# Function's code: GE.space -------------------

GE.space <- function(bckgrnd, GE.occ, add.poly = NULL, save.p = NULL) {
  
  # rgb is color; alpha is a factor of transparency
  # Mcol <- rgb(0.35,0,0.2,alpha = 0.6)
  pal <- c("springgreen3", "tomato")
  
  # if there are too many background points the size of the points is reduced
  if(nrow(bckgrnd)<1000)
    pch.shape = 20
  else
    pch.shape = 42
  
  ## if no polygon is added, the map is created with random background points
  if(is.null(add.poly)) {   
    
    ## Geographic Space:
    
    # create new data-frame with combined coordinates of random background points
    # and occurrence. 
    bckgrnd1 <- data.frame(longitude = bckgrnd[, 1], latitude = bckgrnd[, 2])
    occ1 <- data.frame(longitude = GE.occ[, 1], latitude = GE.occ[, 2])
    data <- cbind(rbind(bckgrnd1[,1:2],occ1[,1:2]),
                  # an extra column is added to differentiate bckgrnd and GE.occ 
                  # (1 for bckgrnd, 2 for GE.occ)
                  c(rep(1,nrow(bckgrnd1)),rep(2,nrow(occ1))))
    # rename columns
    data2 <- data.frame(Longitude = data[, 1], Latitude = data[, 2], 
                        Type = data[,3])
    
    # create plot for G-Space
    p1 <- ggplot(data2, aes(x = Longitude, y = Latitude, color = factor(Type), 
                            shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(values= c("1"=pal[1], "2"= pal[2]), guide = FALSE)
    
    
    ## Environmental Space:
    
    # create new data-frame with combined environmental data of random background 
    # points and occurrence.
    bckgrnd3 <- data.frame(PC1 = bckgrnd[, 3], PC2 = bckgrnd[, 4])
    occ3 <- data.frame(PC1 = GE.occ[, 3], PC2 = GE.occ [, 4])
    data3 <- cbind(rbind(bckgrnd3[,1:2],occ3[,1:2]),
                   # an extra column is added to differentiate bckgrnd and GE.occ 
                   # (1 for bckgrnd, 2 for GE.occ)
                   c(rep(1,nrow(bckgrnd3)),rep(2,nrow(occ3))))
    # rename columns
    data4 <- data.frame(PC1 = data3[, 1], PC2 = data3[, 2], 
                        Type = data3[,3])
    
    # create plot for E-Space
    p2 <- ggplot(data4, aes(x = PC1, y = PC2, 
                            color = factor(Type), shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(name= "Data",
                         labels= c("Background", "Presence"),
                         values= c("1"=pal[1], "2"= pal[2])) +
      theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the 
            # left side, value of 1 to the right, for y, value of 0 puts it to 
            # the bottom, # value of 1 puts it to the top
            legend.justification = c("left", "top"))
    
    # save as plots as a file if save.p is added to the function
    if(!is.null(save.p)) {
      ggarrange(p1, p2, ncol = 2, nrow = 1)
      ggsave(save.p,  width = 24, height = 12, units = "cm",
             dpi = 600, pointsize = 6)
      
    }
    # open plots in R if save.p is not added to the function, no file is saved
    else{
      x11()
      ggarrange(p1, p2, ncol = 2, nrow = 1)
    }
  }
 
  ## if a polygon is added, create a map with the polygon
  else { 

    
    ## Geographical space:
    
    # create new data-frame with the coordinates of the species' occurrence
    occ5 <- data.frame(Longitude = GE.occ[, 1], Latitude = GE.occ[, 2])
    # transform polygon into WGS 84
    M <- spTransform(add.poly, CRS("+proj=longlat +datum=WGS84"))
    # This function turns a map into a data frame that can more easily be 
    # plotted with ggplot2.
    M <- fortify(M)
    # takes id that is a "character" and converts it to a number
    M$id = as.numeric(M$id)
    
    ## Extent of map
    # Sets the boundaries to which the worldmap is cut
    Mext <- extent(add.poly)
    # worldmap is added from the rnaturalearth packages
    world <- ne_countries(scale = "medium", returnclass = "sf")
    class(world)
    
    # create plot for G-Space
    p3 <- ggplot(data = world) +
      geom_sf( ) +
      theme_bw() +
      # alpha -> transparency of polygon, 0.1 = high, 0.5 = medium transparency
      geom_map(map=M, data=M, aes(map_id=id), color= pal[1], fill = pal[1], 
               alpha = 0.4) +
      geom_point(data = occ5, aes(x = Longitude, y = Latitude), size = 2, 
                 shape = 23, fill = pal[2]) +
      #  annotation_scale(location = "br", width_hint = 0.2) +
      #  annotation_north_arrow(location = "br", which_north = "true", 
      #                         pad_x = unit(0.75, "in"), 
      #                         pad_y = unit(0.5, "in"),
      #                         style = north_arrow_fancy_orienteering) +
      coord_sf(xlim = c(Mext[1] - 10 ,Mext[2] + 10), 
               ylim = c(Mext[3] - 10,Mext[4] + 10), expand = FALSE) +
      theme(panel.grid.major = element_line(color = "white"),
            panel.background = element_rect(fill = "aliceblue"))
    
    
    ## Environmental Space:
    # same code as the environmental space above
    bckgrnd3 <- data.frame(PC1 = bckgrnd[, 3], PC2 = bckgrnd[, 4])
    occ3 <- data.frame(PC1 = GE.occ[, 3], PC2 = GE.occ[, 4])
    data3 <- cbind(rbind(bckgrnd3[,1:2],occ3[,1:2]),
                   c(rep(1,nrow(bckgrnd3)),rep(2,nrow(occ3))))
    data4 <- data.frame(PC1 = data3[, 1], PC2 = data3[, 2], 
                        Type = data3[,3])
    
    p4 <- ggplot(data4, aes(x = PC1, y = PC2, 
                            color = factor(Type), shape = factor(Type))) +
      geom_point() +
      scale_shape_manual(values=c(pch.shape,19), guide = FALSE) +
      scale_color_manual(name= "Data",
                         labels= c("Background", "Presence"),
                         values= c("1"=pal[1], "2"= pal[2])) +
      theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the 
            # left side, value of 1 to the right, for y, value of 0 puts it to 
            # the bottom, # value of 1 puts it to the top
            legend.justification = c("left", "top"))
    
    # save as plots as a file if save.p is added to the function
    if(!is.null(save.p)) {
      ggarrange(p3, p4, ncol = 2, nrow = 1)
      ggsave(save.p,  width = 24, height = 12, units = "cm",
             dpi = 600, pointsize = 6)
    }
    
    # open plots in R if save.p is not added to the function, no file is saved
    else{
      x11()
      ggarrange(p3, p4, ncol = 2, nrow = 1)
    }
    
  }
  
}

## read libraries -----------

library(rgdal)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)
# package "sf" needs to be installed
# package "rgeos" needs to be installed


# End