function (registration, segmentation, forward.warp = FALSE) 
{
    if ("plane" %in% names(regi)) {
        if (regi$plane == "sagittal") {
            EPSatlas <- SAGITTALatlas
            atlasIndex <- atlasIndex[atlasIndex$plane == "sagittal", 
                ]
            plane <- "sagittal"
        }
        else {
            atlasIndex <- atlasIndex[atlasIndex$plane == "coronal", 
                ]
            plane <- "coronal"
        }
    }
    else {
        atlasIndex <- atlasIndex[atlasIndex$plane == "coronal", 
            ]
        plane <- "coronal"
    }
    coordinate <- regi$coordinate
    k <- which(abs(coordinate - atlasIndex$mm.from.bregma) == 
        min(abs(coordinate - atlasIndex$mm.from.bregma)))
    plate.info <- EPSatlas$plate.info[[k]]
    scale.factor <- mean(c(dim(regi$transformationgrid$mx)[1]/regi$transformationgrid$height, 
        dim(regi$transformationgrid$mx)[2]/regi$transformationgrid$width))
    outlines <- regi$atlas$outlines
    namelist <- as.numeric(as.vector(plate.info$structure_id))
    colorlist <- plate.info$style
    #neuronregion <- rep(0, length(seg$soma$x))
    #neuroncolor <- rep("#000000", length(seg$soma$x))
    #right.hemisphere <- rep(NA, length(seg$soma$x))
    #for (i in 1:length(outlines)) {
        #temp <- point.in.polygon(seg$soma$x, seg$soma$y, 
            #c(outlines[[i]]$xlT/scale.factor), c(outlines[[i]]$ylT/scale.factor))
        #neuronregion[which(temp == 1)] <- namelist[i]
        #neuroncolor[which(temp == 1)] <- as.character(colorlist[i])
        #right.hemisphere[which(temp == 1)] <- FALSE
        #temp <- point.in.polygon(seg$soma$x, seg$soma$y, 
            #c(outlines[[i]]$xrT/scale.factor), c(outlines[[i]]$yrT/scale.factor))
        #neuronregion[which(temp == 1)] <- namelist[i]
        #neuroncolor[which(temp == 1)] <- as.character(colorlist[i])
        #right.hemisphere[which(temp == 1)] <- TRUE
    #}
    #seg$soma$id <- neuronregion
    #segmentation$soma$color <- neuroncolor
    #segmentation$soma$right.hemisphere <- right.hemisphere
    #dataset <- data.frame(seg$soma[c(1:4, 8:10)])
    
##############################################################################    
 neuronregion <- rep(0,nrow(data))
 neuroncolor <- rep("#000000", nrow(data))     
 right.hemisphere <- rep(NA, nrow(data))
        
 for (i in 1:length(outlines)) {
 	temp <- point.in.polygon(data$xvals, data$yvals, 
            c(outlines[[i]]$xlT/scale.factor2), c(outlines[[i]]$ylT/scale.factor2))
             neuronregion[which(temp == 1)] <- namelist[i]
        neuroncolor[which(temp == 1)] <- as.character(colorlist[i])
        right.hemisphere[which(temp == 1)] <- FALSE
            #x.bound<-c(outlines[[i]]$xr/scale.factor2)
            #y.bound<-c(outlines[[i]]$yr/scale.factor2)
          #rarea[which(temp ==1)]<- polyarea(outlines[i]$xr/scale.factor2, outlines[i]$yr/scale.factor2)

   
   	temp <- point.in.polygon(data$xvals, data$yvals, 
            c(outlines[[i]]$xrT/scale.factor2), c(outlines[[i]]$yrT/scale.factor2))
        neuronregion[which(temp == 1)] <- namelist[i]
        neuroncolor[which(temp == 1)] <- as.character(colorlist[i])
        right.hemisphere[which(temp == 1)] <- TRUE
		#x.bound<-c(outlines[[i]]$xr/scale.factor2)
        #y.bound<-c(outlines[[i]]$yr/scale.factor2)
		#larea[which(temp==1)] <- polyarea(x.bound, y.bound)
 	
 	}

data$region <- neuronregion
data$color <- neuroncolor
data$right.hemisphere <- right.hemisphere
temp.dataset <- data.frame(data[c(1:5)])

unique.regions<-unique(temp.dataset$region)

id<-c(1:length(unique.regions))
total.pixels<-c(1:length(unique.regions))
color<-c(1:length(unique.regions))


for (j in 1:length(unique.regions)) {
	idx<-which(temp.dataset[,3]==unique.regions[[j]])
	num.pixels<-length(idx)
	id[j]<-unique.regions[j]
	total.pixels[j]<-num.pixels		
}
dataset<-data.frame(cbind(id, total.pixels))


dataset$acronym <- rep(NA, length(dataset$id))
    class(dataset$acronym) <- "character"
    dataset$acronym[dataset$id > 0] <- as.character(acronym.from.id(dataset$id[dataset$id > 
        0]))
    dataset$name <- rep(NA, length(dataset$id))
    class(dataset$name) <- "character"
    dataset$name[dataset$id > 0] <- as.character(name.from.id(dataset$id[dataset$id > 
        0]))
    imagename <- substr(basename(registration$outputfile), nchar("Registration_") + 
        1, nchar(basename(registration$outputfile)))
    dataset$image <- rep(imagename, nrow(dataset))
    dataset <- data.frame(animal = rep(NA, nrow(dataset)), AP = rep(registration$coordinate, 
        nrow(dataset)), dataset)    
    
    
    
    
    if (forward.warp == TRUE) {
        if (!(length(registration$transformationgrid$mxF) > 0)) {
            registration <- get.forward.warpRCPP(registration)
        }
        index <- round(scale.factor * cbind(dataset$y, dataset$x))
        index[index == 0] <- 1
        if (length(which(index[, 1] > dim(registration$transformationgrid$mxF)[1]))) {
            index[which(index[, 1] > dim(registration$transformationgrid$mxF)[1]), 
                1] <- dim(registration$transformationgrid$mxF)[1]
        }
        if (length(which(index[, 2] > dim(registration$transformationgrid$mxF)[2]))) {
            index[which(index[, 2] > dim(registration$transformationgrid$mxF)[2]), 
                2] <- dim(registration$transformationgrid$mxF)[2]
        }
        somaX <- registration$transformationgrid$mxF[index]/scale.factor
        somaY <- registration$transformationgrid$myF[index]/scale.factor
        tomecoord <- stereotactic.coordinates(somaX, somaY, registration, 
            inverse = FALSE)
        dataset$ML <- tomecoord$x
        dataset$DV <- tomecoord$y
    }
    #dataset$acronym <- rep(NA, length(dataset$id))
    #class(dataset$acronym) <- "character"
    #dataset$acronym[dataset$id > 0] <- as.character(acronym.from.id(dataset$id[dataset$id > 
       # 0]))
    #dataset$name <- rep(NA, length(dataset$id))
    #class(dataset$name) <- "character"
    #dataset$name[dataset$id > 0] <- as.character(name.from.id(dataset$id[dataset$id > 
       # 0]))
    #imagename <- substr(basename(registration$outputfile), nchar("Registration_") + 
       # 1, nchar(basename(registration$outputfile)))
    #dataset$image <- rep(imagename, nrow(dataset))
    #dataset <- data.frame(animal = rep(NA, nrow(dataset)), AP = rep(registration$coordinate, 
        #nrow(dataset)), dataset)
    #dataset$color <- as.character(dataset$color)
    if (plane == "sagittal") {
        dataset$AP <- dataset$ML
        dataset$ML <- registration$coordinate
    }
    return(dataset)
}
