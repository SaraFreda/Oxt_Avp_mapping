library(wholebrain)
library(rio)

setwd('')

filename <- file.choose(new=FALSE)

#always adjust intensity first, then brain outline... slight bug in wholebrain package for MACs, prefers to adjust outline first


seg <- segment(filename, 2,0.25, filter=seg$filter)

#resize to 120
#brain outline ~3200

# the following code just shows the properties saved in the "seg" variable... plots the location of each detected cell in a xy plot 
names(seg)
names(seg$soma)
plot(seg$soma$x, seg$soma$y, ylim=rev(range(seg$soma$y)), asp=1)

# below is an example of how to hard code the filter for registration... might be useful for future, but can use seg$filter for registration
#myfilter<-structure(list(alim = c(1, 5), threshold.range = c(110, 300), eccentricity = 420L, Max = 255, Min = 0, brain.threshold = 20, resize = 0.04, blur = 4, downsample = 0.25), .Names = c("alim", "threshold.range", "eccentricity", "Max", "Min", "brain.threshold", "resize", "blur", "downsample"))

#registration

########### BE VERY CAREFUL WITH RESIZE PARAMETER ON SEGMENTATION, SHOULD BE a range of ~140-400 ######

#####if you get error transformationgrid$mx[index]: out of  bounds... then you need to change re-size parameters ####

###open webpage of Allen Brain atlas 2008 coronal reference, to use for Bregmma value in registration code
browseURL("https://mouse.brain-map.org/experiment/siv?id=100142143&imageId=102162486&imageType=atlas&initImage=atlas&showSubImage=y&contrast=0.5,0.5,0,255,4")

quartz()
regi<-registration(filename, coordinate=-1.15, filter=seg$filter)


## registration will output two different graphs- orange and purple...
##orange= reference atlas
##purple= your registered image

##if you need to change or add correspondance points to fix registration: use the following code (inputs= your registration object "regi" and number of points you want to add)
#first click the point you want to add on your purple image, then click corresponding area on orange atlas

regi<-add.corrpoints(regi,8)
regi<-change.corrpoints(regi, 1:32)

#re-run registration with added points

regi<- registration(filename, coordinate=-1.15, filter=seg$filter, correspondance=regi)




#inspect segmentation and registration results once you are done adjusting registration to fit brain slice outline


#inspect your registration... this will provide final table of the data with XY coordinates of your detected cells 
dataset<- inspect.registration(regi,seg, forward.warps=TRUE, draw.trans.grid=FALSE)

#if you need to delete any outliers such as cells found outside the brain outline or in ventricles/aqueducts
dataset<- dataset[!dataset$id==0,]
dataset<- dataset[!dataset$id==1009,]
dataset<- dataset[!dataset$id==129,]

#set save directory
#save segmented cells into dataframe named "dataset"
#save segmentation/registration information for that slice
setwd('') #save directory

save(dataset,file='Avp-Cre_F2507_File 3 Slice #7_mcherry.RData')
save(seg,regi, file='Avp-Cre_F2507_File 3 Slice #7_mcherry_Seg-reg.RData')

#use rio package to convert dataframe variable into a .csv file that can be saved locally
convert('Avp-Cre_F2507_File 3 Slice #7_mcherry.RData', 'Avp-Cre_F2507_File 3 Slice #7_mcherry.csv')


# some plotting functions to inspect results:
schematic.plot(dataset, title=TRUE, scale.bar=TRUE, mm.grid=FALSE, pch=21, col=gray(0.1))


# for plotting glassbrain (3d rendering) and downsampling to make the file slightly smaller
#1. load glassbrain data from wholebrain library folder (can use lib.Path() function to find where your data is installed, should be C:/library/3.6/wholebrain/data.... or something close to that)

#2. once file is loaded you will have two variables in your workspace VOLUME and VOLUMESMALL
#3. VOLUME is the high resolution variable you will want to downsample in the following way:
#   VOLUME$v1<-VOLUME$v1[seq(1, nrow(VOLUME$v1), 3),]   ***this will downsample to factor of 3
# repeat the above step, for VOLUME$v2, and VOLUME$v3
#then you can plot using (drawScene.rgl(list(VOLUME)))