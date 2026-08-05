
#setwd('')



dir_path= getwd()


library(dplyr)
#library(purrr)
library(tidyverse)
library(RColorBrewer)
library(wholebrain)

subfolders <- list.dirs(path=dir_path, full.names=TRUE, recursive=FALSE)

#remove any subdirectories you don't want quantified

# compile all .R data files containing "dataset" variable with identified cells

folders <- filtered_subfolders

#initialize empty data frames to store compiled datasets
compiled_data <- data.frame()

#function to load R data files from each subfolder
load_and_compile_r_data_files <- function(folder_path,compiled_data) {
	
	r_data_path <- file.path(folder_path,"R data")
	r_data_files <- list.files(path= r_data_path, full.names=TRUE)
	
	for (file in r_data_files){
		print(paste("Loading file:", file))
		load(file)
		if (exists("dataset")){
			compiled_data <- rbind(compiled_data,dataset)
		} else {
			warning(paste("File",file,"does not contain a variable named 'dataset'"))
		}
	}
	
	return(compiled_data)
}

#function to load R data files in subfolders and change the image name to include file string name (basename)
change_dataset_image_name <- function(folder_path) {
	r_data_path <- file.path(folder_path, "R data")
	r_data_files <- list.files(path = r_data_path, full.names=TRUE)
	
	for (file in r_data_files){
		load(file)
		if (exists("dataset")){
			file_name <- tools::file_path_sans_ext(basename(file))
			dataset$image <- file_name
			save(dataset, file=file)
			print(paste("Modified and saved file:",file))
		} else {
			warning(paste("File",file,"does not contain a variable named 'dataset'"))
		}
		
	}
}



#function to load "dataset" data frame in subfolders and resave it as .csv file
resave_csv_file <- function(folder_path) {
	all_data_folders <- list.dirs(path = folder_path, full.names=TRUE, recursive=TRUE)
	csv_data_path <- grep("csv", all_data_folders, ignore.case=TRUE, value=TRUE)
	r_data_path <- file.path(folder_path, "R data")
	r_data_files <- list.files(path = r_data_path, full.names=TRUE)
	
	for (file in r_data_files) {
		load(file)
		if (exists("dataset")){
			file_name <- tools::file_path_sans_ext(basename(file))
			csv_file_path <- file.path(csv_data_path, paste0(file_name,".csv"))
			write.csv(dataset, file=csv_file_path, row.names=FALSE)
			print(paste("Saved csv file:", csv_file_path))
		} else {
			warning(paste("File:",file, "does not contain a variable named 'dataset'"))
		}
		
	}	
}


######### execute functions to compile data ###########

#iterate through each main folder to modify the dataset$image name
for (folder in folders) {
	change_dataset_image_name(folder)
	}	

#iterate through each main folder to compile data
compiled_data <- data.frame()
for (folder in folders) {
		compiled_data <-load_and_compile_r_data_files(folder, compiled_data)
	}	

#iterate through each main folder, to resave csv files with current R data
for (folder in folders) {
		resave_csv_file(folder)
	}	


#function to get sample names from dataset variable

get_sample_names <- function(df) {
	image_names <- unique(df$image)
	sample_names <- vector("character", length(image_names))
	for (i in seq_along(image_names)) {
		parts <- strsplit(image_names[i],"_")[[1]]
		sample_names[i] <- parts[2]
		}		
	return(unique(sample_names))	
}



#function to return index of either male or female samples
get_sex_index <- function(dataset,sex) {
	opt_sex <- c("female","male")
	criterion <- match.arg(sex, choices= opt_sex)
	
	sample_names <- vector("character", length(dataset$image))
	for (i in seq_along(sample_names)) {
		parts <- strsplit(dataset$image[i],"_")[[1]]
		sample_names[i] <- parts[2]
	}
	
	if (sex == "female") {
		sex_idx <- grep("F", sample_names)
	}
	
	else {
		sex_idx <- grep("M", sample_names)
	}
	
 return(sex_idx)
}


#function to add sample name to each row of dataframe

add_sample_name <- function(dataset) {
	sample_names <- vector("character", length(dataset$image))
	for (i in seq_along(sample_names)) {
		parts <- strsplit(dataset$image[i],"_")[[1]]
		sample_names[i] <- parts[2]
	}
	
	dataset$animal <- sample_names
	return(dataset)	
}

#function to add sex identifier to each row of dataframe
add_sex_identifier <- function(dataset) {
	
	dataset$sex <- ifelse(grepl("M\\d+", dataset$image), "male",
                 ifelse(grepl("F\\d+", dataset$image), "female", NA))
	
	return(dataset)
	
}

#function to split compiled_data/data frame variable based on input parameter

split_compiled_data <- function(df, criterion) {
	valid_criteria <- c("GFP","sex", "sample")
	
	criterion <- match.arg(criterion, choices=valid_criteria)
	
	if (!"image" %in% colnames(df)) {
		stop("The data frame must contain an 'image' column.")
	}
	
	if (criterion == "GFP") {
		gfp_data <- df[grepl("GFP", df$image), ]
		rabies_data <- df[!grepl("GFP", df$image), ]
		
		return(list(GFP=gfp_data, rabies_data=rabies_data))
	}
	
	else if (criterion == "sex") {
		male_idx <- get_sex_index(df,"male")
		female_idx <- get_sex_index(df,"female")
		
		male_data <- df[male_idx,]
		female_data <- df[female_idx,]
		
		return(list(male=male_data, female=female_data))
		
		}
		
	else if (criterion == "sample") {
		unq_samples <- get_sample_names(df)
		sample_data_list <- list()
		
		for (sample in unq_samples) {
			sample_data_list[[sample]] <-df[grepl(paste0("_",sample,"_"),df$image),]	
		}
		return(sample_data_list)
	}
	
}



########## functions to analyze/plot compiled data ###################

### function to plot somas for viral spread
plot_schematic_by_sex <- function(compiled_data, sex) {
	opt_sex <- c("female","male")
	criterion <- match.arg(sex, choices=opt_sex)
	
	
	split_data <- split_compiled_data(compiled_data, "sex")
	female_data <- split_data$female
	male_data <- split_data$male
	
	if (criterion == "female") {
		cell_data <- female_data
	}
	
	else if (criterion == "male") {
		cell_data <- male_data
		
	}

	#k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
	unq.APcoord<-unique(cell_data[,2])
	
	
	#########  code to assign a specific color to each sample ############
	#unq.samples <- unique(cell_data$image)
	#all.sample.names <- sapply(strsplit(unq.samples,"_"), "[[",2)
	#unq.samples <- unique(all.sample.names)
	
	#colors <- viridis(length(unq.samples),option="viridis")
	
	
	#assign colors of each sample in the cell_data variable
	
	
	#for (i in seq_along(unq.samples)) {
		#idx <- grepl(paste0("_",unq.samples[i],"_"),cell_data$image)
		#cell_data$color[idx] <- colors[i]
	#}
	
	############## assign color for all somatic signals ####################
	### all somas will be automatically plotted on the right hemisphere
	
	#cell_data$color <- '#8DC63F' #green color for male samples, oxt
	cell_data$color <- '#005A88' #teal/blue color for female samples, oxt
	
	#cell_data$color <- '#330066' #dark purple for female samples, avp
	#cell_data$color <- '#F08000' #orange for male samples, avp
	
	#cell_data$color <- '#B72467' #magenta/purple for avp samples, all sexes
	#cell_data$color <- '#00ADDC' #light blue color for oxt samples, all sexes

	
	for (j in 1:length(unq.APcoord)) {
	
	idx<-which(cell_data$AP == unq.APcoord[j])
	dataset <- cell_data[idx,]
	dataset$right.hemisphere <- TRUE
	dataset$ML <- abs(dataset$ML)
	schematic.plot(dataset, title=FALSE, scale.bar=TRUE, mm.grid=TRUE, pch=21, col=gray(0.1))
	}	
	
}


#glassbrain from wholebrain package for viral spread
cell_data$right.hemisphere <- TRUE
cell_data$ML <- abs(cell_data$ML)
glassbrain(cell_data, high.res=TRUE, laterality=FALSE)




#code to analyze each sample and determine %total of each parent region (i.e. striatum, pallidum, hypothalamus, etc.)

compiled_data <- compiled_data[!compiled_data$id %in% c(0,129,98,140,81,145), ]    #delete outliers/pixels found in ventricles, aqueducts, etc.
compiled_data <- compiled_data[!compiled_data$id %in% c(512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]  #delete cerebellar regions


calculate_percent_parent_regions <-function(compiled_data) {
	compiled_data_samples <- split_compiled_data(compiled_data,"sample")
	
		isocortex_regions <- c("Frontal","Primary motor","Secondary motor","Somatomotor","Primary somatosensory","Supplemental somatosensory","Gustatory","Visceral area","Auditory","Visual","Anterior cingulate","Prelimbic","Infralimbic","Orbital","Agranular","Retrosplenial","Posterior parietal","Temporal association","Perirhinal area","Ectorhinal","Cerebral cortex")
		
		striatum_regions <- c("Caudoputamen","Nucleus accumbens","Fundus of striatum","Olfactory tubercle","Calleja","Lateral septal nucleus","Septofimbrial","Septohippocampal","Anterior amygdalar","Central amygdalar","Intercalated amygdalar","Medial amygdalar","Bed nucleus of the accessory olfactory","Striatum-like")
		
		pallidum_regions <- c("Globus pallidus","Substantia innominata","Magnocellular","Diagonal band nucleus","Medial septal","Bed nuclei of the stria","Bed nucleus of the anterior","Triangular nucleus of septum")
		
		thalamus_regions <- c("Ventral anterior-lateral","Ventral medial nucleus","Ventral posterolateral","Ventral posteromedial","Subparafascicular","Peripeduncular","Nucleus of reunions","geniculate complex","Parafascicular","habenula","Rhomboid nucleus","Central medial","Paracentral","Perireunensis","Reticular nucleus of the thalamus","Subgeniculate","Mediodorsal nucleus","Intermediodorsal nculeus","reuniens","parataenial","paraventricular nucleus of the thalamus","Anteroventral nucleus","Anteromedial nucleus","Anterodorsal nucleus","Interanteromedial nucleus","Interanterodorsal nucleus","Central lateral nucleus","Paracentral","Subthalamic","Posterior limiting nucleus","Thalamus","nucleus of the thalamus","Suprageniculate")
		
		hypothal_regions <- c("Supraoptic","Nucleus circularis","hypothalamic","Tuberal","Zona incerta","Median eminence","Retrochiasmatic","preoptic","Suprachiasmatic","hypothalamus","Parastrial","premammillary","Dopaminergic A13 group","mammillary","Anteroventral periventricular","Forel","Vascular organ","Subparaventricular zone","Subfornical organ")
		
		midbrain_regions <- c("Midbrain","colliculus","Parabigeminal","Midbrain trigeminal","Substantia nigra","Ventral tegmental","Midbrain reticular","Periaqueductal","Red nucleus","Edinger-Westphal","Cuneiform nuclues","Oculomotor nucleus","Anterior tegmental","Dorsal nucleus raphe","Central linear nucleus raphe","Interfascicular","linear nucleus raphe","Pedunculopontine","Interpeduncular nucleus","Nucleus sagulum","Trochlear nucleus","pretectal","nucleus of the optic tract","Nucleus of Darkschewitsch","Interstitial nucleus of Cajal","Cuneiform nucleus","Lateral terminal nucleus","Nucleus of the posterior commissure","Precommissural nucleus")
		
		hindbrain_regions <- c("Parabrachial","Superior olivary","lateral lemniscus","Dorsal tegmental","Pontine gray","Supratrigeminal","Tegmental reticular nucleus","nucleus of the trigeminal","Motor nucleus of trigeminal","Supragenual","Superior central nucleus raphe","Laterodorsal tegmental","Nucleus incertus","Subceruleus","Nucleus raphe pontis","Area postrema","cochlear","Cuneate","Gracile","spinal nucleus","solitary tract","Pons","Nucleus of the trapezoid body","Pontine reticular","Sublaterodorsal nucleus","Nucleus raphe magnus","Pontine central gray","Koelliker-Fuse","facial motor","Barrington's nucleus","Nucleus raphe pallidus","Medulla","Intermediate reticular nucleus","Locus ceruleus","Gigantocellular","Medial vestibular","Abducens nucleus","Inferior salivatory nucleus","vestibular nucleus","Parvicellular reticular nucleus","Parapyramidal nucleus")
		
		other_regions <- c("fiber tracts","Basic cell groups and regions","Claustrum","Basomedial amygdalar","Posterior amygdalar","Basolateral amygdalar","Lateral amygdalar","Endopiriform")
		
		olfactory_regions <- c("Cortical amygdalar area","olfactory tract","Piriform area","Postpiriform","olfactory bulb","Taenia tecta","Dorsal peduncular","olfactory nucleus","olfactory areas","Piriform-amygdalar")
		
		hippocampus_regions <- c("CA3","CA1","CA2","Postsubiculum","Entorhinal","Parasubiculum","Dentate gyrus","Subiculum","Presubiculum","Induseum griseum","Fasciola cinerea")
		
		all_regions <- c(isocortex_regions,striatum_regions,pallidum_regions,thalamus_regions,hypothal_regions,midbrain_regions,hindbrain_regions,other_regions, hippocampus_regions,olfactory_regions)
		
	
	
	
	norm_parent_regions <- list()
	non_matching_brain_regions <- list()
	for (s in 1:length(compiled_data_samples)) {
		sample_data <- compiled_data_samples[[s]]
		sample_name <- names(compiled_data_samples)[s]
		sample_data$brain_region <- sapply(strsplit(sample_data$name, ","), function(x) x[1])
		sum_total_pixels <- sum(sample_data$total.pixels)
		
		isocortex_idx <- unique(unlist(lapply(isocortex_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		striatum_idx <- unique(unlist(lapply(striatum_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		pallidum_idx <- unique(unlist(lapply(pallidum_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		thalamus_idx <- unique(unlist(lapply(thalamus_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		hypothal_idx <- unique(unlist(lapply(hypothal_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		midbrain_idx <- unique(unlist(lapply(midbrain_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		hindbrain_idx <- unique(unlist(lapply(hindbrain_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		other_idx <- unique(unlist(lapply(other_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		hippocampus_idx <- unique(unlist(lapply(hippocampus_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		olfactory_idx <- unique(unlist(lapply(olfactory_regions, function(region) {which(grepl(region,sample_data$brain_region, ignore.case=TRUE))})))
		
		isocortex_data <- sample_data[isocortex_idx, ]
		striatum_data <- sample_data[striatum_idx, ]
		#need to create exception for striatum so it does not count Basomedial data
		striatum_data <- striatum_data[!striatum_data$id %in% c(334,327), ]
		pallidum_data <- sample_data[pallidum_idx, ]
		thalamus_data <- sample_data[thalamus_idx, ]
		##need to create exception for thalamus, so it does not count hypothalamus data
		thalamus_data <- thalamus_data[!thalamus_data$id %in% c(1097,668,676,684), ]
		hypothal_data <- sample_data[hypothal_idx, ]
		midbrain_data <- sample_data[midbrain_idx, ]
		hindbrain_data <- sample_data[hindbrain_idx, ]
		other_data <- sample_data[other_idx, ]
		hippocampus_data <- sample_data[hippocampus_idx, ]
		##need to create exception for olfactory, so it does not count BAT data
		olfactory_data <- sample_data[olfactory_idx, ]
		olfactory_data <- olfactory_data[!olfactory_data$id %in% c(292), ]
		
		
		isocortex_norm <- sum(isocortex_data$total.pixels)/sum_total_pixels*100
		striatum_norm <- sum(striatum_data$total.pixels)/sum_total_pixels*100
		pallidum_norm <- sum(pallidum_data$total.pixels)/sum_total_pixels*100
		thalamus_norm <- sum(thalamus_data$total.pixels)/sum_total_pixels*100
		hypothal_norm <- sum(hypothal_data$total.pixels)/sum_total_pixels*100
		midbrain_norm <- sum(midbrain_data$total.pixels)/sum_total_pixels*100
		hindbrain_norm <- sum(hindbrain_data$total.pixels)/sum_total_pixels*100
		hippocampus_norm <- sum(hippocampus_data$total.pixels)/sum_total_pixels*100
		olfactory_norm <- sum(olfactory_data$total.pixels)/sum_total_pixels*100
		other_norm <- sum(other_data$total.pixels)/sum_total_pixels*100
		
		normalized_data <- rbind(isocortex_norm, striatum_norm,pallidum_norm,thalamus_norm,hypothal_norm,midbrain_norm,hindbrain_norm,other_norm,hippocampus_norm,olfactory_norm)
		
		non_matching_pattern <- paste(all_regions, collapse = "|")
		non_matching_regions <- sample_data$name[!grepl(non_matching_pattern, sample_data$brain_region, ignore.case = TRUE)]
		norm_parent_regions[[sample_name]] <- normalized_data
		non_matching_brain_regions[[sample_name]] <- non_matching_regions
	}
	
	print(norm_parent_regions)
	print(non_matching_brain_regions)
	return(norm_parent_regions)
	
}



###### code for average heatmap coronal sections, normalized to all total pixels #################


library(wholebrain)
library(tidyverse)
library(viridis)

#first, you need to load two .R data files from wholebrain package:
	#1. EPSatlas.RData
	#2. AtlasIndex.RData
	
load("/Users/snfreda/Documents/R/EPSatlas.RData")
load("/Users/snfreda/Documents/R/atlasIndex.RData")

plot_average_heatmap_coronal_norm_total <- function(APcoord,sex,compiled_data) {
	valid_criteria_sex <- c("female","male", "both")
	
	criterion <- match.arg(criterion, choices=valid_criteria_sex)
	
	if (criterion == "both") {
    projection_data <- compiled_data
  		} 
  	else {
    projection_data <- split_compiled_data(compiled_data, "sex")[[criterion]]
  		}
  		
  	#if (nrow(projection_data) == 0) {
    #stop("No data available for the selected AP coordinate and sex.")
  		#}	

	
	
	####be sure to delete outliers and cerebellar data #####

	projection_data <- projection_data[!projection_data$id %in% c(0,129,98,140,81,145,1009,8), ]    #delete outliers/pixels found in ventricles, aqueducts, etc.
	projection_data <- projection_data[!projection_data$id %in% c(512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]  #delete cerebellar regions
	
	split_projection_data <- split_compiled_data(projection_data,"sample")
	num_samples <- length(names(split_projection_data))
	
	normalized_output_sample <- list()
	for (i in seq_along(split_projection_data)) {
		data <- split_projection_data[[i]]
		sample_name <- names(split_projection_data)[i]
		sum_total_pixels <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		all_norm_pixels = data.frame(brain_region=character(),structure_id = numeric(), norm_pixel = numeric(), sqrt_norm_pixels = numeric(), stringsAsFactors=FALSE)
		
		for (j in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[j], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				sqrt_pix <- sqrt(norm_pix)
				structure_id <- brain_region_data$id[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[j],structure_id=structure_id,norm_pixel=norm_pix,sqrt_norm_pixels=sqrt_pix))
			}
			
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(sqrt_norm_pixels))
		
		normalized_output_sample[[sample_name]] <- sorted_pixels
			
		}
	
	}
	
	#organize data and calculate average normalized pixel for each brain region
	all_brain_regions <- unique(c(unlist(lapply(normalized_output_sample, function(df) df$brain_region))))
	all_structure_ids <- unique(c(unlist(lapply(normalized_output_sample, function(df) df$structure_id))))
	
	combined_data <- do.call(rbind, normalized_output_sample)
	avg_normalized_pixels <- numeric()
	
	for (r in 1:length(all_brain_regions)) {
		data_to_avg <- combined_data[combined_data$brain_region == all_brain_regions[r], ]
		avg_normalized_pixels[r] <- sum(data_to_avg$norm_pixel)/num_samples
		}
	
	#normalized data organized into a data frame for plotting
	
	plate_data_norm <- data.frame(avg_pixels=avg_normalized_pixels,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_norm$animal <- 'NA'
	
	#### if the brain region does not have >0.1 (or 0.05 for avp) norm pixel, it will be removed from the data frame, this will cause the region to be colored black in the heatmap
	
	plate_data_norm <- plate_data_norm[plate_data_norm$avg_pixels > 0.005, ]
	
	######plate information to plot coronal outlines from wholebrain package (need to load two different .R files, see above)
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on pixel normalization
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	acronymlist<-acronym.from.id(structure_id)
	plateData<-data.frame(structure_id=structure_id, color=character(length(acronymlist)), acronym=acronymlist)
	plateData$color=rep('#000000',length(plateData$color))
	
	#generate "heat map" based on avg normalized pixel values
	g1<-ggplot(plate_data_norm, aes(x=animal, y=brain_region, fill=avg_pixels)) + geom_tile() + scale_fill_viridis(option='D',trans="log",name="norm_pixels")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_norm$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
		
	
}





#coronal heatmap that normalizes pixel counts within slice (NOT entire brain)
plot_average_heatmap_coronal_norm_withinslice <- function(APcoord,sex,compiled_data) {
	valid_criteria_sex <- c("female","male", "both")
	
	criterion <- match.arg(sex, choices=valid_criteria_sex)
	slice_projection_data <- compiled_data[compiled_data$AP == APcoord, ]
	
	if (nrow(slice_projection_data) == 0) {
    stop("The selected AP coordinate does not exist in the data. Please choose another AP coordinate.")
  		}
  	
  	if (criterion == "both") {
    projection_data <- slice_projection_data
  		} 
  	else {
    projection_data <- split_compiled_data(slice_projection_data, "sex")[[criterion]]
  		}
  		
  	if (nrow(projection_data) == 0) {
    stop("No data available for the selected AP coordinate and sex.")
  		}	
	
	####be sure to delete outliers and cerebellar data #####

	projection_data <- projection_data[!projection_data$id %in% c(0,129,98,140,81,145,1009,8), ]    #delete outliers/pixels found in ventricles, aqueducts, etc.
	projection_data <- projection_data[!projection_data$id %in% c(512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]  #delete cerebellar regions
	
	split_projection_data <- split_compiled_data(projection_data,"sample")
	num_samples <- length(names(split_projection_data))
	
	normalized_output_sample <- list()
	for (i in seq_along(split_projection_data)) {
		data <- split_projection_data[[i]]
		sample_name <- names(split_projection_data)[i]
		sum_total_pixels <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		all_norm_pixels = data.frame(brain_region=character(),structure_id = numeric(), norm_pixel = numeric(), sqrt_norm_pixels = numeric(), stringsAsFactors=FALSE)
		
		for (j in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[j], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				sqrt_pix <- sqrt(norm_pix)
				structure_id <- brain_region_data$id[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[j],structure_id=structure_id,norm_pixel=norm_pix,sqrt_norm_pixels=sqrt_pix))
			}
			
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(sqrt_norm_pixels))
		
		normalized_output_sample[[sample_name]] <- sorted_pixels
			
		}
	
	}
	
	#organize data and calculate average normalized pixel for each brain region
	all_brain_regions <- unique(c(unlist(lapply(normalized_output_sample, function(df) df$brain_region))))
	all_structure_ids <- unique(c(unlist(lapply(normalized_output_sample, function(df) df$structure_id))))
	
	combined_data <- do.call(rbind, normalized_output_sample)
	avg_normalized_pixels <- numeric()
	
	for (r in 1:length(all_brain_regions)) {
		data_to_avg <- combined_data[combined_data$brain_region == all_brain_regions[r], ]
		avg_normalized_pixels[r] <- sum(data_to_avg$sqrt_norm_pixels)/num_samples
		}
	
	#normalized data organized into a data frame for plotting
	
	plate_data_norm <- data.frame(avg_pixels=avg_normalized_pixels,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_norm$animal <- 'NA'
	
	#### if the brain region does not have >0.1 avg pixel, it will be removed from the data frame, this will cause the region to be colored black in the heatmap
	
	plate_data_norm <- plate_data_norm[plate_data_norm$avg_pixels >= 0.1, ]
	
	######plate information to plot coronal outlines from wholebrain package (need to load two different .R files, see above)
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on pixel normalization
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	acronymlist<-acronym.from.id(structure_id)
	plateData<-data.frame(structure_id=structure_id, color=character(length(acronymlist)), acronym=acronymlist)
	plateData$color=rep('#000000',length(plateData$color))
	
	#generate "heat map" based on avg normalized pixel values
	g1<-ggplot(plate_data_norm, aes(x=animal, y=brain_region, fill=avg_pixels)) + geom_tile() + scale_fill_viridis(option='D',name="norm_pixels")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_norm$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
	
}




######## BAR GRAPH comparison #########


split_relevant_brain_regions <- function(brain_region_vector) {
  brain_regions <- c()
  for (i in seq_along(brain_region_vector)) {
    if (grepl("Globus pallidus", brain_region_vector[i])) {
      brain_region <- gsub(",", "", brain_region_vector[i])
    } else if (grepl("Substantia nigra", brain_region_vector[i])) {
      brain_region <- gsub(",", "", brain_region_vector[i])
    } else {
      brain_region <- strsplit(brain_region_vector[i], ",")[[1]][1]
    }
    brain_regions[i] <- brain_region
  }
  return(brain_regions)  # Add this line to return the result
}

assign_parent_region <- function(brain_region_vector) {
	brain_regions <- sapply(strsplit(brain_region_vector, ","), function(x) x[1])
	
	isocortex_regions <- c("Frontal","Primary motor","Secondary motor","Somatomotor","Primary somatosensory","Supplemental somatosensory","Gustatory","Visceral area","Auditory","Visual","Anterior cingulate","Prelimbic","Infralimbic","Orbital","Agranular","Retrosplenial","Posterior parietal","Temporal association","Perirhinal area","Ectorhinal","Cerebral cortex")
		
		striatum_regions <- c("Caudoputamen","Nucleus accumbens","Fundus of striatum","Olfactory tubercle","Calleja","Lateral septal nucleus","Septofimbrial","Septohippocampal","Anterior amygdalar","Central amygdalar","Intercalated amygdalar","Medial amygdalar","Bed nucleus of the accessory olfactory","Striatum-like")
		
		pallidum_regions <- c("Globus pallidus","Substantia innominata","Magnocellular","Diagonal band nucleus","Medial septal","Bed nuclei of the stria","Bed nucleus of the anterior","Triangular nucleus of septum")
		
		thalamus_regions <- c("Ventral anterior-lateral","Ventral medial nucleus","Ventral posterolateral","Ventral posteromedial","Subparafascicular","Peripeduncular","geniculate complex","Parafascicular","habenula","Rhomboid nucleus","Central medial","Paracentral","Perireunensis","Nucleus of reunions","Reticular nucleus of the thalamus","Subgeniculate","Mediodorsal nucleus","Intermediodorsal nculeus","reuniens","parataenial","paraventricular nucleus of the thalamus","Anteroventral nucleus","Anteromedial nucleus","Anterodorsal nucleus","Interanteromedial nucleus","Interanterodorsal nucleus","Central lateral nucleus","Paracentral","Subthalamic","Posterior limiting nucleus","Thalamus","nucleus of the thalamus","Suprageniculate")
		
		hypothal_regions <- c("Supraoptic","Nucleus circularis","hypothalamic","Tuberal","Zona incerta","Median eminence","Retrochiasmatic","preoptic","Suprachiasmatic","hypothalamus","Parastrial","premammillary","Dopaminergic A13 group","mammillary","Anteroventral periventricular","Forel","Vascular organ","Subparaventricular zone", "Subfornical organ")
		
		midbrain_regions <- c("Midbrain","colliculus","Parabigeminal","Midbrain trigeminal","Substantia nigra","Ventral tegmental","Midbrain reticular","Periaqueductal","Red nucleus","Edinger-Westphal","Cuneiform nuclues","Oculomotor nucleus","Anterior tegmental","Dorsal nucleus raphe","Central linear nucleus raphe","Interfascicular","linear nucleus raphe","Pedunculopontine","Interpeduncular nucleus","Nucleus sagulum","Trochlear nucleus","pretectal","nucleus of the optic tract","Nucleus of Darkschewitsch","Interstitial nucleus of Cajal","Cuneiform nucleus","Lateral terminal nucleus","Nucleus of the posterior commissure","Precommissural nucleus")
		
		hindbrain_regions <- c("Parabrachial","Superior olivary","lateral lemniscus","Dorsal tegmental","Pontine gray","Supratrigeminal","Tegmental reticular nucleus","nucleus of the trigeminal","Motor nucleus of trigeminal","Supragenual","Superior central nucleus raphe","Laterodorsal tegmental","Nucleus incertus","Subceruleus","Nucleus raphe pontis","Area postrema","cochlear","Cuneate","Gracile","spinal nucleus","solitary tract","Pons","Nucleus of the trapezoid body","Pontine reticular","Sublaterodorsal nucleus","Nucleus raphe magnus","Pontine central gray","Koelliker-Fuse","facial motor","Barrington's nucleus","Nucleus raphe pallidus","Medulla","Intermediate reticular nucleus","Locus ceruleus","Gigantocellular","Medial vestibular","Abducens nucleus","Inferior salivatory nucleus","vestibular nucleus","Parvicellular reticular nucleus","Parapyramidal nucleus")
		
		other_regions <- c("fiber tracts","Basic cell groups and regions","Claustrum","Basomedial amygdalar","Posterior amygdalar","Basolateral amygdalar","Lateral amygdalar","Endopiriform")
		
		olfactory_regions <- c("Cortical amygdalar area","olfactory tract","Piriform area","Postpiriform","olfactory bulb","Taenia tecta","Dorsal peduncular","olfactory nucleus","olfactory areas","Piriform-amygdalar")
		
		hippocampus_regions <- c("CA3","CA1","CA2","Postsubiculum","Entorhinal","Parasubiculum","Dentate gyrus","Subiculum","Presubiculum","Induseum griseum","Fasciola cinerea")
		
		all_parent_regions <- rep(NA,length(brain_regions))
		
		isocortex_idx <- unique(unlist(lapply(isocortex_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		striatum_idx <- unique(unlist(lapply(striatum_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		pallidum_idx <- unique(unlist(lapply(pallidum_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		thalamus_idx <- unique(unlist(lapply(thalamus_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		hypothal_idx <- unique(unlist(lapply(hypothal_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		midbrain_idx <- unique(unlist(lapply(midbrain_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		hindbrain_idx <- unique(unlist(lapply(hindbrain_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		other_idx <- unique(unlist(lapply(other_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		hippocampus_idx <- unique(unlist(lapply(hippocampus_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		olfactory_idx <- unique(unlist(lapply(olfactory_regions, function(region) {which(grepl(region,brain_regions, ignore.case=TRUE))})))
		
		all_parent_regions[isocortex_idx] <- "isocortex"
		all_parent_regions[striatum_idx] <- "striatum"
		all_parent_regions[pallidum_idx] <- "pallidum"
		all_parent_regions[thalamus_idx] <- "thalamus"
		all_parent_regions[hypothal_idx] <- "hypothalamus"
		all_parent_regions[midbrain_idx] <- "midbrain"
		all_parent_regions[hindbrain_idx] <- "hindbrain"
		all_parent_regions[other_idx] <- "other"
		all_parent_regions[hippocampus_idx] <- "hippocampus"
		all_parent_regions[olfactory_idx] <- "olfactory"
		
		na_count <- sum(is.na(all_parent_regions))
		print(paste("Number of NA parent regions:",na_count))
		return(all_parent_regions)
	
	
}

plot_bar_projection_sex_comp <- function(compiled_data) {
	
	female_data <- split_compiled_data(compiled_data, "sex")[["female"]]
	male_data <- split_compiled_data(compiled_data, "sex")[["male"]]
	
	female_data_samples <- split_compiled_data(female_data, "sample")
	male_data_samples <- split_compiled_data(male_data, "sample")
	
	female_sample_names <- get_sample_names(female_data)
	male_sample_names <- get_sample_names(male_data)
	
	normalized_female_output <- list()
	normalized_male_output <- list()
	
	
	##### organize female data, based on each sample #####
	
	for (i in seq_along(female_data_samples)) {
		data <- female_data_samples[[i]]
		data <- data[!data$id %in% c(0,129,98,140,81,145,1009,8), ]  # delete outliers in ventricles and fibers in non-descript fiber tracts/boundaries
		data <- data[!data$id %in% c(512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ] #delete cerebellar outliers
		data<- data[!grepl("paraventricular hypothalamic", data$name,ignore.case=TRUE), ] #delete pixels found in PVN
		
		sample_name <- female_sample_names[i]
		sum_total_pixels <- sum(data$total.pixels)
		
		data$parent_region <- assign_parent_region(data$name)
		data$name <- split_relevant_brain_regions(data$name)
		unq_brain_regions <- unique(data$name)
		
		
		all_norm_pixels = data.frame(brain_region=character(),parent_region=character(), norm_pixel = numeric(), stringsAsFactors=FALSE)
		
		for (j in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[j], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				parent_region <-  brain_region_data$parent_region[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[j],parent_region=parent_region,norm_pixel=norm_pix))
				}
			}
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(norm_pixel))
		
		normalized_female_output[[sample_name]] <- sorted_pixels
			
		}
	
	
	##### organize male data, based on each sample #####
	
	
	for (t in seq_along(male_data_samples)) {
		data <- male_data_samples[[t]]
		data <- data[!data$id %in% c(0,129,98,140,81,145,1009,8), ] # delete outliers in ventricles and fibers in non-descript fiber tracts/boundaries
		data <- data[!data$id %in% c(512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ] #delete cerebellar outliers
		data<- data[!grepl("paraventricular hypothalamic", data$name,ignore.case=TRUE), ] #delete pixels found in PVN
		
		sample_name <- male_sample_names[t]
		sum_total_pixels <- sum(data$total.pixels)
		
		data$parent_region <- assign_parent_region(data$name)
		data$name <- split_relevant_brain_regions(data$name)
		unq_brain_regions <- unique(data$name)
		
		
		all_norm_pixels = data.frame(brain_region=character(),parent_region=character(), norm_pixel = numeric(), stringsAsFactors=FALSE)
		
		for (s in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[s], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				parent_region <-  brain_region_data$parent_region[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[s],parent_region=parent_region,norm_pixel=norm_pix))
				}
			}
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(norm_pixel))
		
		normalized_male_output[[sample_name]] <- sorted_pixels
			
		}
	
	
	
	########### reorganize data for plotting #############
	
	combined_data_female <- data.frame()  
  
  for (i in seq_along(normalized_female_output)) {
    df <- normalized_female_output[[i]]
    df$sample <- names(normalized_female_output)[i]  # Add sample name
    combined_data_female <- rbind(combined_data_female, df)  # Combine with the existing data
  }

	combined_data_male <- data.frame() 
  
  for (i in seq_along(normalized_male_output)) {
    df <- normalized_male_output[[i]]
    df$sample <- names(normalized_male_output)[i]  # Add sample name
    combined_data_male <- rbind(combined_data_male, df)  # Combine with the existing data
  }

	combined_data <- rbind(combined_data_female, combined_data_male)
	
	# replace NA or empty parent_region with a placeholder
	combined_data$parent_region[is.na(combined_data$parent_region) | combined_data$parent_region == ""] <- "Unknown"
	
		
	# reshape the data frame so that each sample has a column
	df_wide <- reshape(combined_data, idvar = c("parent_region", "brain_region"), timevar = "sample", direction = "wide")
	df_wide[is.na(df_wide)] <- 0
	
	# calculate the average norm_pixels for each brain_region across samples
	pixel_avg <- rowMeans(df_wide[,-c(1, 2)], na.rm = TRUE)  #ignore the first two columns (parent_region, brain_region)
	
	# sort the data by parent_region and average norm_pixels in descending order
	df_wide$avg_norm_pixels <- pixel_avg
	df_wide <- df_wide[order(df_wide$parent_region, -df_wide$avg_norm_pixels), ]
	
	# select the top 10 brain regions for each parent_region
	top10_indices <- unlist(tapply(1:nrow(df_wide), df_wide$parent_region, function(x) head(x, 10)))
	
	#subset the data to get the top 10 brain regions for each parent_region
	df_top10 <- df_wide[top10_indices, ]
	
	#print numbers in decimal format (instead of scientific notation)
	options(scipen = 999)
	
	print(df_top10)
	return(df_top10)	
	
	
}

##################### code for control analysis ############################
dir_path= getwd()
folders <- list.dirs(path=dir_path, full.names=TRUE, recursive=FALSE)

for (folder in folders) {
	change_dataset_image_name(folder)
	}	

#iterate through each main folder to compile data
compiled_control_data <- data.frame()
for (folder in folders) {
		compiled_control_data <-load_and_compile_r_data_files(folder, compiled_control_data)
	}



#calc_distance <- function(x1,y1,x2,y2) { sqrt((x1-x2)^2 + (y1-y2)^2) }


#control (+) images downsampled 1:16 for anterograde experiments, may need to add dilation to antibody images 
calculate_controls_overlap_oxt <- function(compiled_control_data) {
	
	split_control_data_samples <- split_compiled_data(compiled_control_data, "sample")
	overlapping_cells_calc <- list()
	
	for (sample in names(split_control_data_samples)) {
	
	data <- split_control_data_samples[[sample]]
	#total_eYFP_cells <- nrow(data[grepl("eYFP", data$image),])
	total_eYFP_cells <- nrow(data[grepl("GFP", data$image),])
	
	#eYFP_images <- unique(data[grepl("eYFP", data$image),14])
	eYFP_images <- unique(data[grepl("GFP", data$image),14])
	
	#antibody_images <- unique(data[!grepl("eYFP", data$image),14])
	antibody_images <- unique(data[!grepl("GFP", data$image),14])
		
	all_overlap_data <- data.frame()
	for (eYFP_image in eYFP_images) {
		overlap_data <- data.frame()
		#get corresponding antibody cells
		base_name <- sub("_GFP","",eYFP_image)
		antibody_image <- antibody_images[grepl(base_name, antibody_images)]
		
		
		if (length(antibody_image) > 0) {
			eYFP_cells <- data[data$image == eYFP_image, ]
			antibody_cells <- data[data$image == antibody_image, ]
			eYFP_cells$radius <- sqrt(eYFP_cells$area / pi)*16	   
    		antibody_cells$radius <- sqrt(antibody_cells$area / pi)*16
			
			eYFP_coords <- as.matrix(eYFP_cells[, c("x", "y")])
    			antibody_coords <- as.matrix(antibody_cells[, c("x", "y")])
    			#eYFP_coords <- as.matrix(eYFP_cells[, c("ML", "DV")])
    				#antibody_coords <- as.matrix(antibody_cells[, c("ML", "DV")])
    			combined_coords <- rbind(eYFP_coords, antibody_coords)
    			distance_matrix <- as.matrix(dist(combined_coords))
    			
    			distances <- distance_matrix[1:nrow(eYFP_coords), (nrow(eYFP_coords) + 1):ncol(distance_matrix)]
    			radius_sums <- outer(eYFP_cells$radius, antibody_cells$radius, "+")
    			overlap_idx <- which(distances < radius_sums, arr.ind = TRUE)
    			overlap_data <- rbind(overlap_data, eYFP_cells[overlap_idx[, 1], ])
    			#overlap_data <- overlap_data[!duplicated(overlap_data[, c("x", "y")]), ]
    			overlap_data <- eYFP_cells[unique(overlap_idx[,1]),]
    			all_overlap_data <- rbind(all_overlap_data,overlap_data)
	}
	}
	#overlapping_cells_calc[[sample]] <- nrow(overlap_data)/total_eYFP_cells*100
	overlapping_cells_calc[[sample]] <- nrow(all_overlap_data) / total_eYFP_cells *100		
}
print(overlapping_cells_calc)
return(overlapping_cells_calc)

}


#control (+) images downsampled 1:16, adjust radius if images are downsampled 
calculate_controls_overlap_avp <- function(compiled_control_data) {
	
	split_control_data_samples <- split_compiled_data(compiled_control_data, "sample")
	overlapping_cells_calc <- list()
	
	for (sample in names(split_control_data_samples)) {
	
	data <- split_control_data_samples[[sample]]
	total_mcherry_cells <- nrow(data[grepl("mcherry", data$image),])
	mcherry_images <- unique(data[grepl("mcherry", data$image),14])
	antibody_images <- unique(data[!grepl("mcherry", data$image),14])
		
	overlap_data <- data.frame()
	for (mcherry_image in mcherry_images) {
		
		#get corresponding antibody cells
		base_name <- sub("_mcherry","",mcherry_image)
		antibody_image <- antibody_images[grepl(base_name, antibody_images)]
		
		
		if (length(antibody_image) > 0) {
			mcherry_cells <- data[data$image == mcherry_image, ]
			antibody_cells <- data[data$image == antibody_image, ]
			mcherry_cells$radius <- sqrt(mcherry_cells$area / pi)*16
    		antibody_cells$radius <- sqrt(antibody_cells$area / pi)*16
			
			mcherry_coords <- as.matrix(mcherry_cells[, c("x", "y")])
    			antibody_coords <- as.matrix(antibody_cells[, c("x", "y")])
    			combined_coords <- rbind(mcherry_coords, antibody_coords)
    			distance_matrix <- as.matrix(dist(combined_coords))
    			
    			distances <- distance_matrix[1:nrow(mcherry_coords), (nrow(mcherry_coords) + 1):ncol(distance_matrix)]
    			radius_sums <- outer(mcherry_cells$radius, antibody_cells$radius, "+")
    			overlap_idx <- which(distances < radius_sums, arr.ind = TRUE)
    			overlap_data <- rbind(overlap_data, mcherry_cells[overlap_idx[, 1], ])
    			overlap_data <- overlap_data[!duplicated(overlap_data[, c("x", "y")]), ]
	}
	}
	overlapping_cells_calc[[sample]] <- nrow(overlap_data)/total_mcherry_cells*100
		
}
print(overlapping_cells_calc)
return(overlapping_cells_calc)

}






##################### PLOT HEAT MAP OF "SIMILARITY" ########################
#### will assign colors based on p values for each region, first will calculate p value based on t-test 
##then will color in each region with the assigned heat map color
#more "similar" regions will be colored with bright colors

#input A/P coordinate to plot a single heatmap

load("/Users/snfreda/Documents/R/EPSatlas.RData")
load("/Users/snfreda/Documents/R/atlasIndex.RData")






plot_heatmap_similarity_ttest <- function(df1, df2, APcoord) {
	#APcoord1 <- unique(df1$AP)
	#APcoord2 <- unique(df2$AP)
	
	#all_APcoord <- unique(c(APcoord1,APcoord2))
	
	#criterion <- match.arg(APcoord, choices=all_APcoord)
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline, and in fiber tract regions, as well as cerebellar regions 
	#also remove PVN pixels, will pseudocolor PVN as highly similar
	proj_data1 <-df1[!df1$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	proj_data2 <- df2[!df2$id %in% c(0,129,1009,8,98,140,81,145, 512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	
	proj_data1<- proj_data1[!grepl("paraventricular hypothalamic", proj_data1$name,ignore.case=TRUE), ] #delete pixels found in PVN
	proj_data2 <- proj_data2[!grepl("paraventricular hypothalamic", proj_data2$name,ignore.case=TRUE), ] #delete pixels found in PVN
	
	#normalize pixel counts to total axonal pixels per sample
	proj_data1_samples <- split_compiled_data(proj_data1,"sample")
	proj_data2_samples <- split_compiled_data(proj_data2,"sample")
	
	
	
	
	proj_data1_normalized <- list()       #pixel counts normalized as % total output
	proj_data2_normalized <- list()			# i.e. axon pixels per parent region are divided by total # of axonal pixels in each sample
	
	for (i in seq_along(proj_data1_samples)) {
		data <- proj_data1_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id = numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), sqrt_norm_pixel_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_structure_id)) {
			rg_idx1 <- which(data$id == unq_structure_id[j])
			count <- sum(data$total.pixels[rg_idx1])
			norm_count <- count/total_pixel_count * 100
			sqrt_count <- sqrt(count/total_pixel_count)
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],pixel_count = count, normalized_pixel_count = norm_count,sqrt_norm_pixel_count=sqrt_count))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data1_normalized[[sample_name]] <- sorted_pixel_counts
		}
	
	
	for (i in seq_along(proj_data2_samples)) {
		data <- proj_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id= numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), sqrt_norm_pixel_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			rg_idx2 <- which(data$id == unq_structure_id[j])
			count <- sum(data$total.pixels[rg_idx2])
			norm_count <- count/total_pixel_count * 100
			sqrt_count <- sqrt(count/total_pixel_count)
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], pixel_count = count, normalized_pixel_count = norm_count,sqrt_norm_pixel_count=sqrt_count))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data2_normalized[[sample_name]] <- sorted_pixel_counts
	}
	
	
	proj_data1 <- proj_data1_normalized
	proj_data2 <- proj_data2_normalized
	
	#######calculate p-values (comparing normalized pixel counts in df1 vs df2) as the basis for the "similarity" heatmap values
	
	#1. organize data (combine lists into single data frame) and calculate p-values
	
	oxytocin_df <- bind_rows(proj_data1, .id = "sample") %>% mutate(group = "oxytocin")
	vasopressin_df <- bind_rows(proj_data2, .id = "sample") %>% mutate(group = "vasopressin")
	
	all_data <- bind_rows(oxytocin_df, vasopressin_df)
	
	
	#2. run stats (Wilcox t-test, and FDR correction for multiple comparisons)
	
	all_brain_regions <- unique(all_data$brain_region)
	all_structure_ids <- unique(all_data$structure_id)
	stat_results_list <- list()
	
	for (r in 1:length(all_brain_regions)) {
		sub_data <- all_data[all_data$brain_region == all_brain_regions[r], ]
		 ox_vals <- sub_data$normalized_pixel_count[sub_data$group == "oxytocin"]
  		vas_vals <- sub_data$normalized_pixel_count[sub_data$group == "vasopressin"]
  		
  		# calculate mean values
  		mean_ox <- mean(ox_vals, na.rm = TRUE)
  		mean_vas <- mean(vas_vals, na.rm = TRUE)
  		
  		p_val <- NA
  		if (length(ox_vals) > 0 && length(vas_vals) > 0) {
    		p_val <- wilcox.test(ox_vals, vas_vals)$p.value
  			}
  		
  		#set p-value to 0 if mean normalized pixels is < 0.05 (will help with graph visualization)
  		# Set p-value to 0 if both means are < 0.05
  		if (!is.na(mean_ox) && !is.na(mean_vas) && mean_ox < 0.04 && mean_vas < 0.04) {
   		 p_val <- NA
 			 }
		
	stat_results_list[[r]] <- data.frame(
    brain_region = all_brain_regions[r],
    mean_oxytocin = mean(ox_vals, na.rm = TRUE),
    median_oxytocin = median(ox_vals, na.rm = TRUE),
    mean_vasopressin = mean(vas_vals, na.rm = TRUE),
    median_vasopressin = median(vas_vals, na.rm = TRUE),
    p_value = p_val
  		)	
	}
	
	#combine stat results into one data frame
	results <- do.call(rbind, stat_results_list)
	
	# Adjust p-values for multiple comparisons (Benjamini-Hochberg FDR)
	results$p_adj <- p.adjust(results$p_value, method = "BH")
	
	#Replace 'NA' p-values with negative value, needed for plotting
	results$p_value[is.na(results$p_value)] <- -0.08
	results$p_adj[is.na(results$p_adj)] <- -0.08
	
	

	#organize data for plotting onto coronal alice

	plate_data_pvalues <- data.frame(p_value=results$p_adj,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_pvalues$p_value = results$p_adj
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	
	plate_data_pvalues$animal <- 'NA'
	
	
	#re-scale p-values to adjust heatmap spectrum
	#plate_data_pvalues$neglog_p <- -log10(plate_data_pvalues$p_value)
	#plate_data_pvalues$neglog_p[is.infinite(plate_data_pvalues$neglog_p)] <- NA
	
	######plate information to plot coronal outlines from wholebrain package (need to load two different .R files, see above)
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on p-value comparison
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	acronymlist<-acronym.from.id(structure_id)
	plateData<-data.frame(structure_id=structure_id, color=character(length(acronymlist)), acronym=acronymlist)
	plateData$color=rep('#000000',length(plateData$color))
	
	#generate "heat map" based on p-values
	g1<-ggplot(plate_data_pvalues, aes(x=animal, y=brain_region, fill=p_value)) + geom_tile() + scale_fill_viridis(option='B',name="p-value",na.value="black")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
	#if you want to rescale p-values using log scale
	#g1 <- ggplot(plate_data_pvalues, aes(x = animal, y = brain_region, fill = neglog_p)) + geom_tile() + scale_fill_viridis(option = 'B', name = expression(-log[10]~"(p-value)"),  
    #na.value = "black", limits=c(0,10), oob=squish) +
  	#theme(legend.position = "right")
  	#heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_pvalues$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
	
	
		
}



##### PLOT HEAT MAP OF CONVERGENT AXONAL SITES USING AXON DENSITY (NOT NORMALIZED OUTPUT)

#will normalize to reference region (LHA) and then calculate axonal density using region area
#uses only parent regions i.e. Prelimbic area and not (prelimbia area, layer 2)

#read in csv file with calculated areas for each brain region in Allen CCFv3
csv_dr = '/Users/snfreda/Documents/R/AtlasOntology.csv'
atlas_area_data <- read.csv(csv_dr, header = TRUE, stringsAsFactors = FALSE)

#correct for unfortunate bug in data
re_idx = which(atlas_area_data$id==181)
atlas_area_data$name[re_idx] = 'Nucleus of reunions'

plot_heatmap_convergence <- function(df1, df2, APcoord) {
	#APcoord1 <- unique(df1$AP)
	#APcoord2 <- unique(df2$AP)
	
	#all_APcoord <- unique(c(APcoord1,APcoord2))
	
	#criterion <- match.arg(APcoord, choices=all_APcoord)
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline, and in fiber tract regions, as well as cerebellar regions 
	#also remove PVN pixels, will pseudocolor PVN as highly similar
	proj_data1 <-df1[!df1$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	proj_data2 <- df2[!df2$id %in% c(0,129,1009,8,98,140,81,145, 512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	
	proj_data1<- proj_data1[!grepl("paraventricular hypothalamic", proj_data1$name,ignore.case=TRUE), ] #delete pixels found in PVN
	proj_data2 <- proj_data2[!grepl("paraventricular hypothalamic", proj_data2$name,ignore.case=TRUE), ] #delete pixels found in PVN
	
	#normalize pixel counts to total axonal pixels per sample
	proj_data1_samples <- split_compiled_data(proj_data1,"sample")
	proj_data2_samples <- split_compiled_data(proj_data2,"sample")
	
	
	
	
	proj_data1_normalized <- list()       #pixel counts normalized as % total output
	proj_data2_normalized <- list()			# i.e. axon pixels per parent region are divided by total # of axonal pixels in each sample
	
	ref_region = 194 #structure id for Lateral hypothalamic area (LHA)
	
	for (i in seq_along(proj_data1_samples)) {
		data <- proj_data1_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count_ref <- sum(data$total.pixels[data$id==194])
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id = numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), norm_density = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_structure_id)) {
			pr_parts <- strsplit(unq_brain_regions[j],',')[[1]]
			parent_region <- trimws(pr_parts[1])
			rg_idx1 <- which(data$name == parent_region)
			if (length(rg_idx1) == 0 ) {
				rg_idx1 = grep(parent_region, data$name)
				
			}
			count <- sum(data$total.pixels[rg_idx1])
			area <- sum(atlas_area_data$area[grepl(parent_region,atlas_area_data$name)])
			norm_count <- count/total_pixel_count_ref * 100
			norm_density <- norm_count/area
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],pixel_count = count, normalized_pixel_count = norm_count,norm_density=norm_density))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data1_normalized[[sample_name]] <- sorted_pixel_counts
		}
	
	
	for (i in seq_along(proj_data2_samples)) {
		data <- proj_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count_ref <- sum(data$total.pixels[data$id==194])
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id= numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), norm_density = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			pr_parts <- strsplit(unq_brain_regions[j],',')[[1]]
			parent_region <- trimws(pr_parts[1])
			rg_idx2 <- which(data$name == parent_region)
			if (length(rg_idx2) == 0 ) {
				rg_idx2 = grep(parent_region, data$name)
				
			}
			count <- sum(data$total.pixels[rg_idx2])
			area <- sum(atlas_area_data$area[grepl(parent_region,atlas_area_data$name)])
			norm_count <- count/total_pixel_count_ref * 100
			norm_density <- norm_count/area
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], pixel_count = count, normalized_pixel_count = norm_count,norm_density=norm_density))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data2_normalized[[sample_name]] <- sorted_pixel_counts
	}
	
	
	proj_data1 <- proj_data1_normalized
	proj_data2 <- proj_data2_normalized
	
	#######calculate p-values (comparing normalized pixel counts in df1 vs df2) as the basis for the "similarity" heatmap values
	
	#1. organize data (combine lists into single data frame) and calculate p-values
	
	all_regions <- bind_rows(
  map_dfr(proj_data1, ~ .x[, c("brain_region", "structure_id")]),
  map_dfr(proj_data2, ~ .x[, c("brain_region", "structure_id")])
) %>%
  distinct()
	
	#add in zeros for each sample for all identified brain regions (needed for proper calc of mean later) 
	expand_sample <- function(df, all_regions) {
  df %>%
    right_join(all_regions, by = c("brain_region", "structure_id")) %>%
    mutate(
      pixel_count = replace_na(pixel_count, 0),
      normalized_pixel_count = replace_na(normalized_pixel_count, 0),
      normalized_density = replace_na(norm_density, 0)
    )
}
	
	proj_data1_expanded <- lapply(proj_data1, expand_sample, all_regions = all_regions)
	proj_data2_expanded <- lapply(proj_data2, expand_sample, all_regions = all_regions)

	#combine data into single dataframe
	oxytocin_df <- bind_rows(proj_data1_expanded, .id = "sample") %>% mutate(group = "oxytocin")
	vasopressin_df <- bind_rows(proj_data2_expanded, .id = "sample") %>% mutate(group = "vasopressin")
	
	oxytocin_df$norm_density[is.infinite(oxytocin_df$norm_density)] <- NA
	vasopressin_df$norm_density[is.infinite(vasopressin_df$norm_density)] <- NA
	
	all_data <- bind_rows(oxytocin_df, vasopressin_df)
	
	
	#2. run stats (Wilcox t-test, and FDR correction for multiple comparisons)
	
	all_brain_regions <- unique(all_data$brain_region)
	all_structure_ids <- unique(all_data$structure_id)
	stat_results_list <- list()
	
	for (r in 1:length(all_brain_regions)) {
		sub_data <- all_data[all_data$brain_region == all_brain_regions[r], ]
		sub_data$norm_density[is.na(sub_data$norm_density)] <- 0
		 ox_vals <- sub_data$norm_density[sub_data$group == "oxytocin"]
  		vas_vals <- sub_data$norm_density[sub_data$group == "vasopressin"]
  		
  		mean_ox_norm_pixel <- mean(sub_data$normalized_pixel_count[sub_data$group == "oxytocin"], na.rm=TRUE)
  		mean_vas_norm_pixel <- mean(sub_data$normalized_pixel_count[sub_data$group == "vasopressin"], na.rm=TRUE)
  		
  		# calculate mean values
  		mean_ox <- mean(ox_vals, na.rm = TRUE)
  		mean_vas <- mean(vas_vals, na.rm = TRUE)
  		
  		p_val <- NA
  		 if (length(ox_vals[!is.na(ox_vals) & ox_vals != 0]) > 0 &&
      	length(vas_vals[!is.na(vas_vals) & vas_vals != 0]) > 0) {
    
    		p_val <- wilcox.test(ox_vals, vas_vals)$p.value
  		}
  		
  		#set p-value to NA if mean normalized pixels is < 0.3 (will help with graph visualization)
  		# Set p-value to NA if both means are < 0.3
  		if (!is.na(mean_ox_norm_pixel) && !is.na(mean_vas_norm_pixel)) {
 			if (mean_ox_norm_pixel < 0.3 && mean_vas_norm_pixel < 0.3) {
   			 p_val <- NA
 			 }
			}
		
	stat_results_list[[r]] <- data.frame(
    brain_region = all_brain_regions[r],
    mean_oxytocin = mean(ox_vals, na.rm = TRUE),
    median_oxytocin = median(ox_vals, na.rm = TRUE),
    mean_vasopressin = mean(vas_vals, na.rm = TRUE),
    median_vasopressin = median(vas_vals, na.rm = TRUE),
    p_value = p_val
  		)	
	}
	
	#combine stat results into one data frame
	results <- do.call(rbind, stat_results_list)
	
	# Adjust p-values for multiple comparisons (Benjamini-Hochberg FDR or Holm method)
	results$p_adj <- p.adjust(results$p_value, method = "BH")
	
	#Replace 'NA' p-values with negative value, needed for plotting
	results$p_value[is.na(results$p_value)] <- -0.05
	results$p_adj[is.na(results$p_adj)] <- -0.05
	
	results$p_adj[is.nan(results$p_adj)] <- 1

	#organize data for plotting onto coronal alice

	plate_data_pvalues <- data.frame(p_value=results$p_adj,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_pvalues$p_value = results$p_adj
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	
	plate_data_pvalues$animal <- 'NA'
	
	
	#re-scale p-values to adjust heatmap spectrum
	#plate_data_pvalues$neglog_p <- -log10(plate_data_pvalues$p_value)
	#plate_data_pvalues$neglog_p[is.infinite(plate_data_pvalues$neglog_p)] <- NA
	
	######plate information to plot coronal outlines from wholebrain package (need to load two different .R files, see above)
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on p-value comparison
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	acronymlist<-acronym.from.id(structure_id)
	plateData<-data.frame(structure_id=structure_id, color=character(length(acronymlist)), acronym=acronymlist)
	plateData$color=rep('#000000',length(plateData$color))
	
	#generate "heat map" based on p-values
	g1<-ggplot(plate_data_pvalues, aes(x=animal, y=brain_region, fill=p_value)) + geom_tile() + scale_fill_viridis(option='B',name="p-value",na.value="black")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
	#if you want to rescale p-values using log scale
	#g1 <- ggplot(plate_data_pvalues, aes(x = animal, y = brain_region, fill = neglog_p)) + geom_tile() + scale_fill_viridis(option = 'B', name = expression(-log[10]~"(p-value)"),  
    #na.value = "black", limits=c(0,10), oob=squish) +
  	#theme(legend.position = "right")
  	#heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_pvalues$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
	
	
		
}

####UPDATED HEATMAP SIMILARITY

#helper function to prevent missing zeros in the comparison of two brain areas:
#add in zeros for each sample for all identified brain regions (needed for proper calc of mean later) 
	expand_sample <- function(df, all_regions) {
  df %>%
    right_join(all_regions, by = c("brain_region", "structure_id")) %>%
    mutate(
      pixel_count = replace_na(pixel_count, 0),
      normalized_pixel_count = replace_na(normalized_pixel_count, 0),
      sqrt_norm_pixel_count = replace_na(sqrt_norm_pixel_count, 0)
    )
}


plot_heatmap_similarity_ttest <- function(df1, df2, APcoord) {
	#APcoord1 <- unique(df1$AP)
	#APcoord2 <- unique(df2$AP)
	
	#all_APcoord <- unique(c(APcoord1,APcoord2))
	
	#criterion <- match.arg(APcoord, choices=all_APcoord)
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline, and in fiber tract regions, as well as cerebellar regions 
	#also remove PVN pixels, will pseudocolor PVN as highly similar
	proj_data1 <-df1[!df1$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	proj_data2 <- df2[!df2$id %in% c(0,129,1009,8,98,140,81,145, 512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	
	proj_data1<- proj_data1[!grepl("paraventricular hypothalamic", proj_data1$name,ignore.case=TRUE), ] #delete pixels found in PVN
	proj_data2 <- proj_data2[!grepl("paraventricular hypothalamic", proj_data2$name,ignore.case=TRUE), ] #delete pixels found in PVN
	
	#normalize pixel counts to total axonal pixels per sample
	proj_data1_samples <- split_compiled_data(proj_data1,"sample")
	proj_data2_samples <- split_compiled_data(proj_data2,"sample")
	
	
	
	
	proj_data1_normalized <- list()       #pixel counts normalized as % total output
	proj_data2_normalized <- list()			# i.e. axon pixels per parent region are divided by total # of axonal pixels in each sample
	
	for (i in seq_along(proj_data1_samples)) {
		data <- proj_data1_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id = numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), sqrt_norm_pixel_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_structure_id)) {
			rg_idx1 <- which(data$id == unq_structure_id[j])
			count <- sum(data$total.pixels[rg_idx1])
			norm_count <- count/total_pixel_count * 100
			sqrt_count <- sqrt(count/total_pixel_count)
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],pixel_count = count, normalized_pixel_count = norm_count,sqrt_norm_pixel_count=sqrt_count))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data1_normalized[[sample_name]] <- sorted_pixel_counts
		}
	
	
	for (i in seq_along(proj_data2_samples)) {
		data <- proj_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id= numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), sqrt_norm_pixel_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			rg_idx2 <- which(data$id == unq_structure_id[j])
			count <- sum(data$total.pixels[rg_idx2])
			norm_count <- count/total_pixel_count * 100
			sqrt_count <- sqrt(count/total_pixel_count)
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], pixel_count = count, normalized_pixel_count = norm_count,sqrt_norm_pixel_count=sqrt_count))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data2_normalized[[sample_name]] <- sorted_pixel_counts
	}
	
	
	proj_data1 <- proj_data1_normalized
	proj_data2 <- proj_data2_normalized
	
	#######calculate p-values (comparing normalized pixel counts in df1 vs df2) as the basis for the "similarity" heatmap values
	
	
	
	#1. organize data (combine lists into single data frame) and calculate p-values
	
	all_regions <- bind_rows(
  map_dfr(proj_data1, ~ .x[, c("brain_region", "structure_id")]),
  map_dfr(proj_data2, ~ .x[, c("brain_region", "structure_id")])
) %>%
  distinct()
	
	
		
	proj_data1_expanded <- lapply(proj_data1, expand_sample, all_regions = all_regions)
	proj_data2_expanded <- lapply(proj_data2, expand_sample, all_regions = all_regions)

	#combine data into single dataframe
	oxytocin_df <- bind_rows(proj_data1_expanded, .id = "sample") %>% mutate(group = "oxytocin")
	vasopressin_df <- bind_rows(proj_data2_expanded, .id = "sample") %>% mutate(group = "vasopressin")
	
	
	all_data <- bind_rows(oxytocin_df, vasopressin_df)
		
	#2. run stats (Wilcox t-test, and FDR correction for multiple comparisons)
	
	all_brain_regions <- unique(all_data$brain_region)
	all_structure_ids <- unique(all_data$structure_id)
	stat_results_list <- list()
	
	for (r in 1:length(all_brain_regions)) {
		sub_data <- all_data[all_data$brain_region == all_brain_regions[r], ]
		 ox_vals <- sub_data$normalized_pixel_count[sub_data$group == "oxytocin"]
  		vas_vals <- sub_data$normalized_pixel_count[sub_data$group == "vasopressin"]
  		
  		# calculate mean values
  		mean_ox <- mean(ox_vals, na.rm = TRUE)
  		mean_vas <- mean(vas_vals, na.rm = TRUE)
  		
  		p_val <- NA
  		if (length(ox_vals) > 0 && length(vas_vals) > 0) {
    		p_val <- wilcox.test(ox_vals, vas_vals)$p.value
  			}
  		
  		#set p-value to 0 if mean normalized pixels is < 0.05 (will help with graph visualization)
  		# Set p-value to 0 if both means are < 0.05
  		if (!is.na(mean_ox) && !is.na(mean_vas) && mean_ox < 0.05 && mean_vas < 0.05) {
   		 p_val <- NA
 			 }
		
	stat_results_list[[r]] <- data.frame(
    brain_region = all_brain_regions[r],
    mean_oxytocin = mean(ox_vals, na.rm = TRUE),
    median_oxytocin = median(ox_vals, na.rm = TRUE),
    mean_vasopressin = mean(vas_vals, na.rm = TRUE),
    median_vasopressin = median(vas_vals, na.rm = TRUE),
    p_value = p_val
  		)	
	}
	
	#combine stat results into one data frame
	results <- do.call(rbind, stat_results_list)
	
	# Adjust p-values for multiple comparisons (Benjamini-Hochberg FDR)
	results$p_adj <- p.adjust(results$p_value, method = "BH")
	
	#Replace 'NA' p-values with negative value, needed for plotting
	results$p_value[is.na(results$p_value)] <- -0.08
	results$p_adj[is.na(results$p_adj)] <- -0.08
	
		

	#organize data for plotting onto coronal alice

	plate_data_pvalues <- data.frame(p_value=results$p_adj,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_pvalues$p_value = results$p_adj
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	
	plate_data_pvalues$animal <- 'NA'
	
	
	#re-scale p-values to adjust heatmap spectrum
	#plate_data_pvalues$neglog_p <- -log10(plate_data_pvalues$p_value)
	#plate_data_pvalues$neglog_p[is.infinite(plate_data_pvalues$neglog_p)] <- NA
	
	######plate information to plot coronal outlines from wholebrain package (need to load two different .R files, see above)
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on p-value comparison
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	acronymlist<-acronym.from.id(structure_id)
	plateData<-data.frame(structure_id=structure_id, color=character(length(acronymlist)), acronym=acronymlist)
	plateData$color=rep('#000000',length(plateData$color))
	
	#generate "heat map" based on p-values
	g1<-ggplot(plate_data_pvalues, aes(x=animal, y=brain_region, fill=p_value)) + geom_tile() + scale_fill_viridis(option='B',name="p-value",na.value="black")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
	#if you want to rescale p-values using log scale
	#g1 <- ggplot(plate_data_pvalues, aes(x = animal, y = brain_region, fill = neglog_p)) + geom_tile() + scale_fill_viridis(option = 'B', name = expression(-log[10]~"(p-value)"),  
    #na.value = "black", limits=c(0,10), oob=squish) +
  	#theme(legend.position = "right")
  	#heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_pvalues$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
	
	
		
}


####### NORMALIZE DATA BASED ON PARENT REGION ############
##### NORMALIZE AXON COUNTS FOR EACH PARENT REGION (E.G. cortex,hypothalamus) instead of normalizing to total pixel counts

normalize_axon_output_per_parent <- function(axon_df) {
	proj_data <-axon_df[!axon_df$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675,688), ]
	#proj_data <- axon_df[!axon_df$id %in% c(688),] #delete pixels in cerebral cortex for cortical analyses
	proj_data <- proj_data[!grepl("paraventricular hypothalamic", proj_data$name,ignore.case=TRUE), ] #delete pixels found in PVN
	
	# #bin brain region names into larger subregion catergories (i.e. Superior colliculus, motor related, intermediate gray layer --> Superior colliculus)
	proj_data$name <- split_relevant_brain_regions(proj_data$name)
	
	proj_data_samples <- split_compiled_data(proj_data,"sample")
	
	proj_data_normalized <- list()
	
	for (i in seq_along(proj_data_samples)) {
    data <- proj_data_samples[[i]]
    sample_name <- strsplit(data$image[1], "_")[[1]][2]

    # assign parent regions using existing function
    data$parent_region <- assign_parent_region(data$name)
    
     
     # sum total pixels per parent region for normalization
    parent_totals <- data %>%
      group_by(parent_region) %>%
      summarise(parent_total = sum(total.pixels), .groups = "drop")
      
      # merge totals back into data
    data <- left_join(data, parent_totals, by = "parent_region")
    
     # combine repeated regions (sum across all occurrences of the same brain region)
    region_summary <- data %>%
      group_by(brain_region = name, parent_region) %>%
      summarise(
        pixel_count = sum(total.pixels, na.rm = TRUE),
        parent_total = unique(parent_total),
        .groups = "drop"
      )
     
     # normalize within each parent region
    region_summary <- region_summary %>%
      mutate(
        normalized_pixel_count = (pixel_count / parent_total) * 100,
        sqrt_norm_pixel_count = sqrt(pixel_count / parent_total)
      )
      
       # sort for convenience
    sorted_pixel_counts <- region_summary %>%
    arrange(parent_region, desc(pixel_count))

    proj_data_normalized[[sample_name]] <- sorted_pixel_counts
  }       
      
     # combine all samples
  all_outputs <- bind_rows(proj_data_normalized, .id = "sample") %>%
    mutate(group = "outputs")

 # pivot to wide format: each sample is a column
  all_outputs_avg <- all_outputs %>%
    select(sample, brain_region, parent_region, normalized_pixel_count) %>%
    pivot_wider(
      names_from = sample,
      values_from = normalized_pixel_count,
      values_fill = 0
    ) %>%
    rowwise() %>%
    mutate(avg_norm_output = mean(c_across(-c(brain_region, parent_region)))) %>%
    ungroup()
    
  return(as.data.frame(all_outputs_avg))
    

	}	


########### normalize cortical axons for mPFC ########
normalize_mPFC_axons <- function(axon_df) {
	proj_data <-axon_df[!axon_df$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675,688), ]
	proj_data <- proj_data[!grepl("paraventricular hypothalamic", proj_data$name,ignore.case=TRUE), ] #delete pixels found in PVN
	
	#select only for pixels in the medial prefrontal cortex
	mpfc_regions <- c("Infralimbic","Prelimbic","Anterior cingulate area")
	mpfc_idx <- grepl(paste(mpfc_regions, collapse = "|"),
  		proj_data$name, ignore.case = TRUE)
  		
  	df_mpfc <- proj_data[mpfc_idx, ]
  	
  	mpfc_data_samples <- split_compiled_data(df_mpfc,"sample")
	mpfc_data_normalized <- list()
	
	for (i in seq_along(mpfc_data_samples)) {
		data <- mpfc_data_samples[[i]]
    	sample_name <- strsplit(data$image[1], "_")[[1]][2]
    	
    	total_mpfc_pixels = sum(data$total.pixels)
    	
    	mpfc_unique_areas = unique(data$name)
    	norm_vals <- numeric(length(mpfc_unique_areas))
    	for (j in 1:seq_along(mpfc_unique_areas)) {
    		area_idx <- data$name == mpfc_unique_areas[j]
    		pix_count_norm <- sum(data$total.pixels[area_idx])/total_mpfc_pixels*100
    		norm_vals[j] <- pix_count_norm
    		
    	}
    	
    	df_norm <- data.frame(mpfc_area=mpfc_unique_areas, norm_percent=norm_vals, stringsAsFactors = FALSE)
    	
    	mpfc_data_normalized[[sample_name]] <- df_norm
    	
		}
		
		mpfc_long <- bind_rows(
  mpfc_data_normalized,
  .id = "sample"
)

mpfc_wide <- mpfc_long %>%
  pivot_wider(
    names_from  = sample,
    values_from = norm_percent
  )

stopifnot(
  nrow(mpfc_wide) == length(unique(mpfc_long$mpfc_area))
)

mpfc_wide$mean_norm <- rowMeans(mpfc_wide[ , -1], na.rm = TRUE)

mpfc_wide_reshape <- mpfc_wide[
  order(mpfc_wide$mean_norm, decreasing = TRUE),
]
		
		#print(mpfc_wide)
		return(as.data.frame(mpfc_wide_reshape))
		

}
	
	
	

	






###########PLOTTING RECURRENT VS. TOP-DOWN VS. BOTTOM-UP MODULATION ###############################



### will be plotted as a heatmap with each brain region colored based on its type of "circuit"
## heatmap colors will be calculated as an average ratio  normalized_output+e/normalized_inputs+e (e= small number to account for areas that may have zero inputs, but axonal output)


##### PREPARE OUTPUT/AXON DATA FOR PLOTS #########

normalize_axon_output <- function(axon_df) {
	proj_data <-axon_df[!axon_df$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
	proj_data_samples <- split_compiled_data(proj_data,"sample")
	####UNLIKE BAR PLOTS, DOES NOT DELETE PVN AXONS 
	
	
	
	proj_data_normalized <- list()       #pixel counts normalized as % total output # i.e. axon pixels per parent region are divided by total # of axonal pixels in each sample
	
	
	for (i in seq_along(proj_data_samples)) {
		data <- proj_data_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_pixel_count <- sum(data$total.pixels)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		pixel_counts = data.frame(brain_region=character(),structure_id = numeric(), pixel_count = numeric(), normalized_pixel_count = numeric(), sqrt_norm_pixel_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_structure_id)) {
			rg_idx1 <- which(data$id == unq_structure_id[j])
			count <- sum(data$total.pixels[rg_idx1])
			norm_count <- count/total_pixel_count * 100
			sqrt_count <- sqrt(count/total_pixel_count)
			pixel_counts <- rbind(pixel_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],pixel_count = count, normalized_pixel_count = norm_count,sqrt_norm_pixel_count=sqrt_count))
			
		}
		
		sorted_pixel_counts <- pixel_counts %>% arrange(desc(pixel_count))
		
		proj_data_normalized[[sample_name]] <- sorted_pixel_counts
		}
		
	axon_data <- proj_data_normalized
	all_outputs <- bind_rows(axon_data, .id = "sample") %>% mutate(group = "outputs")
	
	#combine data into single dataframe, wide format:
	#each row = unique brain region with structure_id
	#each column = normalized_cell_count for one sample
	#add column for average normalized_cell_count for each region
	
	all_outputs_avg <- all_outputs %>%
  	select(sample, brain_region, structure_id, normalized_pixel_count) %>%	
  	pivot_wider(
    	names_from = sample,
    	values_from = normalized_pixel_count,
    	values_fill = 0  # fill missing with 0
 	 ) %>%
  
  	# Compute average across all sample columns
  	rowwise() %>%
  	mutate(avg_norm_output = mean(c_across(-c(brain_region, structure_id)))) %>%
  	ungroup()
  	all_outputs_avg <- as.data.frame(all_outputs_avg)
  	
  	return(all_outputs_avg)
	
	
}







### PREPARE INPUT DATA FOR PLOTS #########
calc_distance <- function(x1,y1,x2,y2) { sqrt((x1-x2)^2 + (y1-y2)^2) }

identify_starter_cells <- function(df) {
	GFP_images <- unique(df[grepl("GFP", df$image),14])
	rabies_images <- unique(df[!grepl("GFP", df$image),14])
	starter_cells <- data.frame()
	
	for (GFP_image in GFP_images) {
		
		#get corresponding non-GFP(rabies mcherry) cells
		base_name <- sub("_GFP","",GFP_image)
		rabies_image <- rabies_images[grepl(base_name, rabies_images)]
		
		if (length(rabies_image) > 0) {
			#filter data frame for specific GFP slice
			
			gfp_cells <- df[df$image == GFP_image,]
			rabies_cells <- df[df$image %in% rabies_image,]
			
			#calculate radius for each cell (assuming circularity)
			gfp_cells$radius <- sqrt(gfp_cells$area / pi) * 4  
			rabies_cells$radius <- sqrt(rabies_cells$area / pi) * 4			
			
			#check if each cell is overlapping and then save to data frame
			for (i in 1:nrow(gfp_cells)) {
				for (j in 1:nrow(rabies_cells)) {
				distance <- calc_distance(gfp_cells$x[i],gfp_cells$y[i],rabies_cells$x[j],rabies_cells$y[j])
				if (distance < (gfp_cells$radius[i] + rabies_cells$radius[j])) {
					starter_cells <- rbind(starter_cells, rabies_cells[j, ])
					#break  #move to next rabies cell to ensure each is only counted once
					}
				}
			}
		}
	}
	starter_cells <- starter_cells %>% distinct()   #remove duplicate rows in case any cells are added multiple times 
	starter_cells <- starter_cells[!starter_cells$id %in% c(523,181,668,797,756),] #remove any outliers
	row.names(starter_cells)<-c(1:nrow(starter_cells))
	return(starter_cells)
}

#function to remove starter cells and GFP+ cells from compiled data
remove_starter_cells <- function(df) {
	starter_cells <- identify_starter_cells(df)
	cleaned_compiled_data <- anti_join(df, starter_cells, by = c("animal", "x","y","area","ML","DV","acronym","name","image"))
	cleaned_rabies_data <- cleaned_compiled_data[!grepl("GFP", cleaned_compiled_data$image), ]
	return(cleaned_rabies_data)
}

normalize_rabies_input <- function(input_df) {
	#remove potential outliers i.e. located in ventricle or outside registered brain outline and remove starter cells/GFP cells
	cell_data <-input_df[!input_df$id %in% c(0,129,303,131,23,1015,276,284,411,1009,8), ] ###UNLIKE RABIES FIGURES, WILL REMOVE CELLS FROM FIBER TRACTS & BASIC CELL GROUPS
	cleaned_cell_data <- remove_starter_cells(cell_data)
	#split df into each sample, to normalize each sample separately
	cell_data_samples <- split_compiled_data(cleaned_cell_data,"sample")
	
	cell_data_normalized_counts <- list()       #cell counts normalized as % total input 
	
	for (i in seq_along(cell_data_samples)) {
		data <- cell_data_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		cell_counts = data.frame(brain_region=character(),structure_id = numeric(), cell_count = numeric(), normalized_cell_count = numeric(), sqrt_norm_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(grepl(unq_brain_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],cell_count = count, normalized_cell_count = norm_count,sqrt_norm_cell_count=sqrt_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		cell_data_normalized_counts[[sample_name]] <- sorted_cell_counts
		}
		
		cell_data <- cell_data_normalized_counts
		
		
	#combine data into single dataframe, wide format:
	#each row = unique brain region with structure_id
	#each column = normalized_cell_count for one sample
	#add column for average normalized_cell_count for each region
	
	all_inputs <- bind_rows(cell_data, .id = "sample") %>% mutate(group = "inputs")
	
	all_inputs_avg <- all_inputs %>%
  	select(sample, brain_region, structure_id, normalized_cell_count) %>%	
  	pivot_wider(
    	names_from = sample,
    	values_from = normalized_cell_count,
    	values_fill = 0  # fill missing with 0
 	 ) %>%
  
  	# Compute average across all sample columns
  	rowwise() %>%
  	mutate(avg_norm_count = mean(c_across(-c(brain_region, structure_id)))) %>%
  	ungroup()
  	all_inputs_avg <- as.data.frame(all_inputs_avg)
  	
  	return(all_inputs_avg)
		
}


########## PLOT HEATMAP OF RECURRENCE #####
### need to load wholebrain files before running function
load("/Users/snfreda/Documents/R/EPSatlas.RData")
load("/Users/snfreda/Documents/R/atlasIndex.RData")
library(wholebrain) ###warning need to load this package first, then dplry or else code will break
library(dplyr)
library(tidyverse)
library(scales)
library(RColorBrewer)

plot_heatmap_recurrence <- function(df_inputs, df_outputs, APcoord) {
	norm_inputs <- normalize_rabies_input(df_inputs)
	norm_outputs <- normalize_axon_output(df_outputs)
	
	
	#calculate ratio of outputs/inputs as a measure of recurrent network
	# R = 1 recurrent network
	# R > 1 "top-down" network
	# R < 1 bottom-up network
	# will include epsilon value to prevent divison by 0
	
	epsilon <- 0.5
	norm_outputs_clean <- norm_outputs %>%
  dplyr::select(
    brain_region,
    structure_id,
    avg_norm_output = avg_norm_output   # explicit name
  )

norm_inputs_clean <- norm_inputs %>%
  dplyr::select(
    brain_region,
    structure_id,
    avg_norm_count = avg_norm_count     # explicit name
  )
	
	
merged_df <- full_join(norm_outputs_clean, norm_inputs_clean,
                       by = c("brain_region", "structure_id"))
	
	
	
		#Replace missing with 0 (in case a region is only present in one df)	
	merged_df <- merged_df %>%
  		mutate(
    		avg_norm_output  = ifelse(is.na(avg_norm_output), 0, avg_norm_output),
    		avg_norm_count = ifelse(is.na(avg_norm_count), 0, avg_norm_count)
  		)
  		
  	#calculate ratio:
  	merged_df <- merged_df %>%
  		mutate(
    			ratio = (avg_norm_output + epsilon) / (avg_norm_count + epsilon),
    			log_ratio = log(ratio),
    
    			# apply cutoff: if either average < 0.1, replace with NA
    			ratio = ifelse(avg_norm_count <= 0.1 & avg_norm_output <= 0.1, NA, ratio),
    		log_ratio = ifelse(avg_norm_count <= 0.1 & avg_norm_output <= 0.1, NA, log_ratio)
  		)
  
  #all_structure_ids <- merged_df$structure_id
  #all_brain_regions <- merged_df$brain_region
	
	
	# coronal plot information
	coordinate <- APcoord
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    plate.info <- EPSatlas$plate.info[[k]]              
    scale.factor = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	
	#create "plateData" data frame, this will store all the structures found in your coronal slice outline and ultimately assign a color based on ratio comparison
	
	structure_id <- as.numeric(as.vector(plate.info$structure_id))
	plateData<-data.frame(structure_id=structure_id)
	plateData$color <- rep("#000000", nrow(plateData))
	

	#assign colors based on ratio/log_ratio values
	n_colors <- 50
	palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(n_colors)

	# Reverse palette so that high values are red
	palette <- rev(palette)	
	scaled_values <- as.integer(cut(merged_df$log_ratio, breaks = n_colors))
	merged_df$ratio_color <- palette[scaled_values]
	
	plateData <- plateData %>%
  left_join(merged_df %>% select(structure_id, ratio_color), by = "structure_id") %>%
  mutate(color = ifelse(is.na(ratio_color), "#000000", ratio_color)) %>%
  select(-ratio_color)  # remove the extra column
	
	
		
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#808080'
	
	
	#plot coronal heatmap
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
	
   
    lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "white", col=plateData$color[N], lwd = 1)
       })	
	
	
	
}

plot_recurrent_categories <- function(df_inputs, df_outputs) {
	norm_inputs <- normalize_rabies_input(df_inputs)
	norm_outputs <- normalize_axon_output(df_outputs)
	
	
	#calculate ratio of outputs/inputs as a measure of recurrent network
	# 0.67 < R < 1.5 recurrent network
	# R > 1.5 top-down network
	# R < 0.67 bottom-up network
	# will include epsilon value to prevent divison by 0
	
	epsilon <- 0.5
	norm_outputs_clean <- norm_outputs %>%
  dplyr::select(
    brain_region,
    structure_id,
    avg_norm_output = avg_norm_output   # explicit name
  )

norm_inputs_clean <- norm_inputs %>%
  dplyr::select(
    brain_region,
    structure_id,
    avg_norm_count = avg_norm_count     # explicit name
  )
	
	
merged_df <- full_join(norm_outputs_clean, norm_inputs_clean,
                       by = c("brain_region", "structure_id"))
	
	
	
		#Replace missing with 0 (in case a region is only present in one df)	
	merged_df <- merged_df %>%
  		mutate(
    		avg_norm_output  = ifelse(is.na(avg_norm_output), 0, avg_norm_output),
    		avg_norm_count = ifelse(is.na(avg_norm_count), 0, avg_norm_count)
  		)
  		
  	#calculate ratio:
  	merged_df <- merged_df %>%
  		mutate(
    			ratio = (avg_norm_output + epsilon) / (avg_norm_count + epsilon),
    			log_ratio = log(ratio),
    
    			# apply cutoff: if either average < 0.1, replace with NA
    			ratio = ifelse(avg_norm_count <= 0.05 & avg_norm_output <= 0.05, NA, ratio),
    		log_ratio = ifelse(avg_norm_count <= 0.05 & avg_norm_output <= 0.05, NA, log_ratio)
  		)
  		
  	#adjust paraventricular nucleus ratios to be 1:
  	pvn_idx <- grep("Paraventricular hypothalamic", merged_df$brain_region)
  	merged_df$ratio[pvn_idx] <- 1
  	merged_df$log_ratio[pvn_idx] <- 0
  	
  	#plot region proportions or categories (top-down, bottom-up, balanced/recurrent) as part of whole plot
  	#catergorize merged_df
  	#conservative approach to category cut-offs, set uses +/-1 or two-fold changes (log2) to represent category
  	merged_df <- merged_df %>%
  		filter(!is.na(ratio)) %>%
  		mutate(ratio_category = case_when(
    	ratio < 0.67 ~ "< 0.67",
    	ratio >= 0.67 & ratio <= 1.5 ~ "0.67 – 1.5",
    	ratio > 1.5 ~ "> 1.5"
  	))
	
	#count brain region in each of three categories:
	counts_df <- merged_df %>%
  	count(ratio_category) %>%
  	mutate(proportion = n / sum(n))
  	
  	quartz()
  	ggplot(counts_df, aes(x = "All Regions", y = proportion, fill = ratio_category)) +
  	geom_bar(stat = "identity", width = 0.5) +
  	scale_y_continuous(labels = scales::percent_format()) +
  	labs(x = "", y = "Proportion of Brain Regions", fill = "Ratio Category") +
  	theme_minimal()
	
}



###misc code to plot smoothed ratio distribution oxt vs avp:
#df <- merged_df %>% filter(!is.na(ratio))
# Define a common x-range that fully covers the data
#xrange <- range(df$ratio) + c(-0.5, 0.5)
## Compute density for each group manually
#densities <- df %>%
  #group_by(group) %>%
  #summarise(
    #dens = list(density(ratio, from = xrange[1], to = xrange[2], n = 512))
  #) %>%
  #mutate(
    #x = map(dens, "x"),
    #y = map(dens, "y")
  #) %>%
  #select(-dens) %>%
  #unnest(cols = c(x, y))
#to plot smoothed density curves:
#ggplot(densities, aes(x = x, y = y, color = group, fill = group)) +
  #geom_line(linewidth = 1.1) +
  #geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2, linetype = 0) +
  #geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  #labs(
   # x = "Output/Input Ratio",
   # y = "Density",
   # title = "Smoothed Distribution of Ratios",
   # color = "Group",
   # fill = "Group"
 # ) +
 # theme_classic(base_size = 14)
 
 
 
 
 
 ##### Generate supplemental tables, organize axon data for oxt/avp systems 
 #Column 1: all detected brain regions, organized alphabetically 
 #Column 2: Raw pixel counts
 #Column 3: Mean normalized %output
 #Column 4: SEM
 #Column 5: Norm mean % (male - female)
 ###repeat for Oxt & Avp system
 library(dplyr)
 library(tidyr)
 
 generate_anterograde_table <- function(oxt_outputs_df, avp_outputs_df) {
 	#delete outliers, ventricles, non-descript regions, and PVN
 	oxt_data <- oxt_outputs_df[!oxt_outputs_df$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
 	avp_data <- avp_outputs_df[!avp_outputs_df$id %in% c(0,129,1009,8,98,140,81,145,512, 10722, 10720, 10674, 10713, 10710, 10708 ,10672, 10711,10692, 10690,10687,10689,10677,10675), ]
 	
 	oxt_data<- oxt_data[!grepl("paraventricular hypothalamic", oxt_data$name,ignore.case=TRUE), ] #delete pixels found in PVN
 	avp_data<- avp_data[!grepl("paraventricular hypothalamic", avp_data$name,ignore.case=TRUE), ] #delete pixels found in PVN
 	
 	#organize pixel data, first find all unique brain regions (in this case use larger subregion category, not smaller subregion classification)
 	# (i.e. Inferior colliculus, external nucleus --> inferior colliclus )
 	
 	oxt_data$name <- split_relevant_brain_regions(oxt_data$name)
 	avp_data$name <- split_relevant_brain_regions(avp_data$name)
 	
 	#add column for sample identifier 
 	oxt_data <- add_sample_name(oxt_data)
 	avp_data <- add_sample_name(avp_data)
 	
 	#add column for sex identifier
 	oxt_data <- add_sex_identifier(oxt_data)
 	avp_data <- add_sex_identifier(avp_data)
 	
 	#remove columns that are not needed (AP, id, acronym) and add system label
 	#oxt <- oxt_data %>%
  	#select(animal, total.pixels, name, sex) %>%
  	#mutate(system = "Oxt")
 	
 	#avp <- avp_data %>%
 	#select(animal,total.pixels,name,sex) %>%
 	#mutate(system = "Avp")
 	
 	#combine into one data frame
 	data_all <- bind_rows(oxt_data, avp_data)
 	
 	#rename column
 	names(data_all)[names(data_all) == "name"] <- "brain_region"
 	
 	#organize based on unique name (brain region) x animal x system
 	all_regions <- unique(data_all$brain_region)
 	
 	
 	oxt_samples <- unique(oxt_data$animal)
	oxt_data_samples <- split_compiled_data(oxt_data, "sample")
	
	normalized_oxt_output <- list()
	##### organize data, based on each sample #####
	
	for (i in seq_along(oxt_data_samples)) {
		data <- oxt_data_samples[[i]]
		
		sample_name <- oxt_samples[i]
		sum_total_pixels <- sum(data$total.pixels)
		
		#data$parent_region <- assign_parent_region(data$name)
		#data$name <- split_relevant_brain_regions(data$name)
		unq_brain_regions <- unique(data$name)
		
		
		all_norm_pixels = data.frame(brain_region=character(), pix_count=numeric(), norm_pixel = numeric(), stringsAsFactors=FALSE)
		
		for (j in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[j], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				#parent_region <-  brain_region_data$parent_region[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[j],pix_count = pix_data, norm_pixel=norm_pix))
				}
			}
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(norm_pixel))
		
		normalized_oxt_output[[sample_name]] <- sorted_pixels
			
		}
 	
 	combined_data_oxt <- data.frame()  
  
  for (i in seq_along(normalized_oxt_output)) {
    df <- normalized_oxt_output[[i]]
    df$sample <- names(normalized_oxt_output)[i]  # Add sample name
    combined_data_oxt <- rbind(combined_data_oxt, df)  # Combine with the existing data
  }
 	
 	# reshape the data frame so that each sample has a column
	df_wide_oxt <- reshape(combined_data_oxt, idvar = c("brain_region"), timevar = "sample", direction = "wide")
	df_wide_oxt[is.na(df_wide_oxt)] <- 0
	
	# calculate the average norm_pixels, total pixel counts and SEM for each brain_region across samples
	pixel_avg_oxt <- rowMeans(df_wide_oxt[, grepl("norm_pixel", names(df_wide_oxt))], na.rm = TRUE)  
	total_pixel_count <- rowSums(df_wide_oxt[, grepl("pix_count", names(df_wide_oxt))], na.rm = TRUE)
	
	pixel_avg_oxt_female <- rowMeans(df_wide_oxt[, grepl("norm_pixel.F", names(df_wide_oxt))], na.rm = TRUE)
	pixel_avg_oxt_male <- rowMeans(df_wide_oxt[, grepl("norm_pixel.M", names(df_wide_oxt))], na.rm = TRUE)
	pixel_avg_diff <- pixel_avg_oxt_female - pixel_avg_oxt_male
	
	n_samples_oxt <- length(oxt_samples)
	sd_per_row <- apply(df_wide_oxt[, grepl("norm_pixel", names(df_wide_oxt))], 1, sd, na.rm = TRUE)
	sem_pixel_oxt <- sd_per_row / sqrt(n_samples_oxt)
	
	oxt_summary_data <- data.frame(brain_region= df_wide_oxt$brain_region, total_pix_count=total_pixel_count, norm_pixel = pixel_avg_oxt, sem_norm_pixel = sem_pixel_oxt, pixel_avg_diff_sex = pixel_avg_diff)
	
	
	
	######################## do the same for vasopressin axon data ###############################
	
	avp_samples <- unique(avp_data$animal)
	avp_data_samples <- split_compiled_data(avp_data, "sample")
	
	normalized_avp_output <- list()
	##### organize data, based on each sample #####
	
	for (j in seq_along(avp_data_samples)) {
		data <- avp_data_samples[[j]]
		
		sample_name <- avp_samples[j]
		sum_total_pixels <- sum(data$total.pixels)
		
		#data$parent_region <- assign_parent_region(data$name)
		#data$name <- split_relevant_brain_regions(data$name)
		unq_brain_regions <- unique(data$name)
		
		
		all_norm_pixels = data.frame(brain_region=character(), pix_count=numeric(), norm_pixel = numeric(), stringsAsFactors=FALSE)
		
		for (k in 1:length(unq_brain_regions)) {
			brain_region_data <- data[data$name == unq_brain_regions[k], ]
			
			if (nrow(brain_region_data) >0) {
				pix_data <- sum(brain_region_data$total.pixels)
				norm_pix <- pix_data / sum_total_pixels*100
				#parent_region <-  brain_region_data$parent_region[1]
				
				all_norm_pixels <- rbind(all_norm_pixels, data.frame(brain_region=unq_brain_regions[k],pix_count = pix_data, norm_pixel=norm_pix))
				}
			}
		
		sorted_pixels <- all_norm_pixels %>% arrange(desc(norm_pixel))
		
		normalized_avp_output[[sample_name]] <- sorted_pixels
			
		}
 	
 	combined_data_avp <- data.frame()  
  
  for (i in seq_along(normalized_avp_output)) {
    df <- normalized_avp_output[[i]]
    df$sample <- names(normalized_avp_output)[i]  # Add sample name
    combined_data_avp <- rbind(combined_data_avp, df)  # Combine with the existing data
  }
 	
 	# reshape the data frame so that each sample has a column
	df_wide_avp <- reshape(combined_data_avp, idvar = c("brain_region"), timevar = "sample", direction = "wide")
	df_wide_avp[is.na(df_wide_avp)] <- 0
	
	# calculate the average norm_pixels, total pixel counts and SEM for each brain_region across samples
	pixel_avg_avp <- rowMeans(df_wide_avp[, grepl("norm_pixel", names(df_wide_avp))], na.rm = TRUE)  
	total_pixel_count <- rowSums(df_wide_avp[, grepl("pix_count", names(df_wide_avp))], na.rm = TRUE)
	
	pixel_avg_avp_female <- rowMeans(df_wide_avp[, grepl("norm_pixel.F", names(df_wide_avp))], na.rm = TRUE)
	pixel_avg_avp_male <- rowMeans(df_wide_avp[, grepl("norm_pixel.M", names(df_wide_avp))], na.rm = TRUE)
	pixel_avg_diff <- pixel_avg_avp_female - pixel_avg_avp_male
	
	n_samples_avp <- length(avp_samples)
	sd_per_row <- apply(df_wide_avp[, grepl("norm_pixel", names(df_wide_avp))], 1, sd, na.rm = TRUE)
	sem_pixel_avp <- sd_per_row / sqrt(n_samples_avp)
	
	avp_summary_data <- data.frame(brain_region= df_wide_avp$brain_region, total_pix_count=total_pixel_count, norm_pixel = pixel_avg_avp, sem_norm_pixel = sem_pixel_avp, pixel_avg_diff_sex = pixel_avg_diff)

	
### combine oxt and avp summary data into one table and order by brain region in alphabetical order
	 combined_summary_data <- full_join(oxt_summary_data, avp_summary_data,
                              by = "brain_region",
                              suffix = c(".oxt", ".avp"))
     
     #rearrange to be alphabetical
     combined_summary_data <- combined_summary_data %>%
  arrange(brain_region)                        
	
	
	
	return(combined_summary_data)
 		
 	}
 
 
 


