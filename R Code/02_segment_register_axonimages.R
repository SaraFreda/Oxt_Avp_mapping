 ##libraries to run the analysis:

library(wholebrain) #used for registration of brain slices
library(OpenImageR) #used to upload image files (usually .tif file) and run edge detection 
library(rio) #used to convert .R files to .csv files 
library(dbscan) #used for avp projection mapping, post-processing

#script: quantify_axon_pixels.R will modify wholebrain package to register segmented axons
###see below for workflow to modify the existing package for the session

## read image into R workspace, might take some time because it is a large image
filename <- file.choose(new=FALSE)

image<-readImage(filename)


edges.image<-edge_detection(image,method='Sobel',conv_mode='full',approx='TRUE')

## if you want to view the image and detected axons (edges) you can use the following code

imageShow(edges.image)


## manually convert image to binary based on pixel gradient, consistent for whole sample

edges.image[which(edges.image<0.17)]=0
edges.image[which(edges.image>0.17)]=1

## organize pixel gradient/intensities into data frame 

data=data.frame(xvals=c(1:length(which(edges.image==1))), yvals=c(1:length(which(edges.image==1))), region=c(1:length(which(edges.image==1))))

data=data.frame(xvals=rep(0,length(which(edges.image==1))), yvals=rep(0,length(which(edges.image==1))))

edges.idx=which(edges.image==1,arr.ind=TRUE)

data$xvals=edges.idx[,2]
data$yvals=edges.idx[,1]

#if you want to plot and visualize the detected pixels in scatter plot format
quartz()
plot(data$xvals, data$yvals, ylim=rev(range(data$yvals)), pch='.', col="black", asp=1)

#post-processing for avp data, only needed for images with high background (or very low axonal signal):
coords <- data.frame(x = data$xvals, y = data$yvals)
clustering <- dbscan(coords, eps = 1.95, minPts = 6)
clustered_data <- coords[clustering$cluster != 0, ]

###########################################



#to visualize detected axons, post-processing
quartz()
plot(clustered_data, xlim=range(data$xvals), ylim=rev(range(data$yvals)),pch=".", col="black",asp=1)

data = clustered_data
colnames(data) <- c("xvals", "yvals")

# first, need to segment and register brain using wholebrain code

#make sure you have the correct working directory:
setwd('/Users/snfreda/Documents/R')

#segmentation code, only need to capture brain outline and determine Bregma coordinate i.e. do not need to find signal or use for thresholding
seg <- segment(filename, 2,0.25, filter=seg$filter)
browseURL("https://mouse.brain-map.org/experiment/siv?id=100142143&imageId=102162486&imageType=atlas&initImage=atlas&showSubImage=y&contrast=0.5,0.5,0,255,4")

#then register brain slice
quartz()
regi<-registration(filename, coordinate=1.00, filter=seg$filter)

#if need to add or change correspondance points for slice registration, use the following code:
regi<-add.corrpoints(regi,10)
regi<-change.corrpoints(regi, 11:18)
regi<- registration(filename, coordinate=1.00, filter=seg$filter, correspondance=regi)

#now that registration is complete, inspect registration and add pixels for each registered region
### USE	trace(get.cell.ids, edit=TRUE) AND REPLACE CODE WITH 'quantify_axon_pixels.R' at the START of your R session (i.e. should only need to do this once per R session)

scale.factor2 <- mean(dim(regi$transformationgrid$mx)/c(regi$transformationgrid$height, 
        regi$transformationgrid$width))
        
trace(get.cell.ids, edit=TRUE) ###this will modify the existing wholebrain package to quantify axon pixels for the current R session

###RUN Wholebrain code to find axon pixels in registered brain regions, depending on axonal density this function may take a few seconds
dataset<-get.cell.ids(regi,seg,forward.warp=FALSE)

################## to manually check results #####################

numPaths <- regi$atlas$numRegions
outlines <- regi$atlas$outlines
scale.factor <- mean(dim(regi$transformationgrid$mx)/c(regi$transformationgrid$height, 
        regi$transformationgrid$width))
xMax <- max(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
        na.rm = TRUE) * (1/scale.factor)
xMin <- min(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
        na.rm = TRUE) * (1/scale.factor)
yMax <- max(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
        na.rm = TRUE) * (1/scale.factor)
yMin <- min(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
        na.rm = TRUE) * (1/scale.factor)
quartz()
plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
        xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", regi$coordinate, 
            "mm"), font.main = 1)

#plot(data$xvals, data$yvals, ylim = c(yMax, yMin), xlim = c(xMin, xMax),  pch='.', col="darkmagenta", asp=1)
#points(data$xvals, data$yvals, ylim = c(yMax, yMin), xlim = c(xMin, xMax), pch=20, cex=0.2, col="darkmagenta",asp=1)            
            
  lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xrT/scale.factor, outlines[[x]]$yrT/scale.factor, 
            border = "lightskyblue2", lwd=0.4)
    })
    lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xlT/scale.factor, outlines[[x]]$ylT/scale.factor, 
            border = "lightskyblue2", lwd=0.4)
    })
    
 #### if you need to delete any outliers such as pixels found outside the brain slice or in ventricles
dataset<- dataset[!dataset$id==0,]
dataset<- dataset[!dataset$id==129,]
dataset<- dataset[!dataset$id==81,]

#if you need to delete rows and re-number:
dataset<- dataset[-c(1:3),]     #this will delete rows 1:3
row.names(dataset)<-c(1:nrow(dataset))	

##save segmentation, registation, and quantification (dataset) files

setwd('/Users/snfreda/Desktop/Avp Projection analysis/Quantification/Revisions')
save(dataset,file='Avp_F2503_File 3 Slice 1_axons.RData')
save(seg,regi, file='Avp_F2503_File 1 Slice 10_axons_Seg-Reg.RData')

convert('Avp_F2503_File 3 Slice 1_axons.RData', 'Avp_F2503_File 3 Slice 1_axons.csv')






################ plotting density contour plots for a single slice ############


###load library
library(tidyverse)
##convert x and y values from "original" data set so that it can be plotted in the registered ML,DV axis 
##first obtain forward warps:
regi<-get.forward.warpRCPP(regi)


##option 1: limit KDE (kernel density estimation) to axon pixel data limits within each slice (may result in "cut-off" contours)

#then calculate "new" x and y values to transform onto new grid system
scale.factor <- mean(dim(regi$transformationgrid$mx)/c(regi$transformationgrid$height, 
        regi$transformationgrid$width))
index <- round(scale.factor * cbind(data$yvals, data$xvals))  #may generate an error if axon pixels are outside registration lines
newX <- regi$transformationgrid$mxF[index]/scale.factor
newY <- regi$transformationgrid$myF[index]/scale.factor

#tomecoord <- stereotactic.coordinates(newX, newY, regi, 
            #inverse = FALSE)
# xt<- tomecoord$x
# yt <- tomecoord$y        


#data <- data.frame(newX = newX,newY = newY)
#save(data,file='Oxt_F234_File 1 Slice 0_axonpixels.RData')

#now initialize data for contour plots
data.plot<-ggplot(data, aes(x=newX, y=newY))

#establish x and y axis limits:
xMax <- max(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
        na.rm = TRUE) * (1/scale.factor)
xMin <- min(c(regi$transformationgrid$mx, regi$transformationgrid$mxF), 
        na.rm = TRUE) * (1/scale.factor)
yMax <- max(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
        na.rm = TRUE) * (1/scale.factor)
yMin <- min(c(regi$transformationgrid$my, regi$transformationgrid$myF), 
        na.rm = TRUE) * (1/scale.factor)
                  
numPaths <- regi$atlas$numRegions
outlines <- regi$atlas$outlines        
        
# next, call quartz window and initialize plot with title ***must initialize plot first in order for lapply funct to work
quartz()
plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
        xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", regi$coordinate, 
            "mm"), font.main = 1)   

#plot contours 

data.plot + geom_density_2d_filled(contour_var="ndensity") + theme(legend.position="none") + 
coord_fixed(ratio=1, ylim=c(yMax,yMin),xlim=c(xMin,xMax))

           
#data.plot + geom_density_2d_filled(contour_var="ndensity") + ylim(rev(range(newY)))  + theme(legend.position="none") + 
#coord_fixed(ratio=1, ylim=c(yMax,yMin),xlim=c(xMin,xMax))

#then overlay outlines
 lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xr/scale.factor, outlines[[x]]$yr/scale.factor, 
            border = "black", col=regi$atlas$col2[x])
            })
    lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xl/scale.factor, outlines[[x]]$yl/scale.factor, 
            border = "black", col=regi$atlas$col2[x])
    })           
                
########option 2: 
### perform KDE using registered atlas outlines as boundry for contour estimation (may stretch data estimation slightly)      

#prepare data for plot contours
xrange <- c(xMin, xMax)
yrange <- c(yMax, yMin)

# Perform KDE (kernel density estimate) across the whole brain atlas/boundaries
kde <- kde2d(
  x = newX,
  y = newY,
  n = 300,  # grid resolution
  lims = c(xrange[1], xrange[2], yrange[1], yrange[2])
)

# Convert to data frame for ggplot
kde.df <- with(kde, expand.grid(x = x, y = y))
kde.df$z <- as.vector(kde$z)

# next, call quartz window and initialize plot with title ***must initialize plot first in order for lapply funct to work
quartz()
plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
        xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", regi$coordinate, 
            "mm"), font.main = 1)
            
#plot contours with atlas outlines

ggplot(kde.df, aes(x, y, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level))) +
  coord_fixed(ratio=1, xlim=xrange, ylim=yrange, expand=FALSE) +
  theme(legend.position="none")

#then overlay outlines (these atlas outlines from CCF)
 #lapply(1:numPaths, function(x) {
        #polygon(outlines[[x]]$xr/scale.factor, outlines[[x]]$yr/scale.factor, 
            #border = "black", col=regi$atlas$col2[x])
            #})
    #lapply(1:numPaths, function(x) {
        #polygon(outlines[[x]]$xl/scale.factor, outlines[[x]]$yl/scale.factor, 
            #border = "black", col=regi$atlas$col2[x])
    #})
  
 #outlines from raw registered image 
  lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xrT/scale.factor, outlines[[x]]$yrT/scale.factor, 
            border = "black")
    })
    lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xlT/scale.factor, outlines[[x]]$ylT/scale.factor, 
            border = "black")
    })
    
   
  
  

############## plotting axon pixels and brain slice, transformed to allen brain atlas grid  ########

#first, you need to load two .R data files saved under "Documents/R" folder 
	#1. EPSatlas.RData
	#2. AtlasIndex.RData
	
load("/Users/snfreda/Documents/R/EPSatlas.RData")
load("/Users/snfreda/Documents/R/atlasIndex.RData")
 
numPaths <- regi$atlas$numRegions
outlines <- regi$atlas$outlines
regi<-get.forward.warpRCPP(regi)


scale.factor <- mean(dim(regi$transformationgrid$mx)/c(regi$transformationgrid$height, 
        regi$transformationgrid$width))
index <- round(scale.factor * cbind(data$yvals, data$xvals))  #may generate an error if axon pixels are outside registration lines
newX <- regi$transformationgrid$mxF[index]/scale.factor
newY <- regi$transformationgrid$myF[index]/scale.factor
 
k <- which(abs(regi$coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(regi$coordinate - atlasIndex$mm.from.bregma[1:132])))
plate <- atlasIndex$plate.id[which(abs(regi$coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(regi$coordinate - 
                  atlasIndex$mm.from.bregma)))]
plate.info <- EPSatlas$plate.info[[k]]              
#scale.factor = 0.9579832
#numPaths <- EPSatlas$plates[[k]]@summary@numPaths
#bregmaY = 200
#bregmaX = 5640
#xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2                  

plate.info$style<-"#f2f2f2"
#this code changes the fiber tracts to be a dark gray
idx.fibers<-which(plate.info$structure_id==1009)
plate.info$style[idx.fibers]="#d7d7d7"
idx.ventricles<-which(plate.info$structure_id==129)
plate.info$style[idx.ventricles]="#000000"

 quartz()
 plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
        xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", regi$coordinate, 
            "mm"), font.main = 1)
       
   lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xr/scale.factor, outlines[[x]]$yr/scale.factor, 
            border = "black", col=plate.info$style[x], lwd= 0.5)
            })
    lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xl/scale.factor, outlines[[x]]$yl/scale.factor, 
            border = "black", col=plate.info$style[x], lwd=0.5)
    })
    points(newX, newY, ylim = c(yMax, yMin), xlim = c(xMin, xMax), pch=20, cex=0.14, col="darkmagenta",asp=1)


#col="darkmagenta" for avp axons
#col = chartreuse4 for oxt axons


###### just to plot regi outlines #####
regi<-get.forward.warpRCPP(regi)
numPaths <- regi$atlas$numRegions
outlines <- regi$atlas$outlines


quartz()
plot(c(xMin, xMax), c(yMin, yMax), ylim = c(yMax, yMin), 
        xlim = c(xMin, xMax), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", regi$coordinate, 
            "mm"), font.main = 1)
 lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xr/scale.factor, outlines[[x]]$yr/scale.factor, 
            border = "black", lwd= 0.5)
            })
    lapply(1:numPaths, function(x) {
        polygon(outlines[[x]]$xl/scale.factor, outlines[[x]]$yl/scale.factor, 
            border = "black", lwd=0.5)
    })           
            
            

#col=green4, for Oxt proj

    
