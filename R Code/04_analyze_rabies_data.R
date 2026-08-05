setwd('') #set working directory to folder containing analyzed data from rabies workflow
dir_path= getwd()


library(dplyr)

subfolders <- list.dirs(path=dir_path, full.names=TRUE, recursive=FALSE)


# compile all .R data files containing "dataset" variable with identified cells
# do this for both rabies mcherry and GFP channels

folders <- filtered_subfolders

#initialize empty data frames to store compiled datasets
compiled_data <- data.frame()

#function to load R data files from each subfolder
load_and_compile_r_data_files <- function(subfolder_path,compiled_data) {
	
	r_data_path <- file.path(subfolder_path,"R data")
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
change_dataset_image_name <- function(subfolder_path) {
	r_data_path <- file.path(subfolder_path, "R data")
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
resave_csv_file <- function(subfolder_path) {
	all_data_folders <- list.dirs(path = subfolder_path, full.names=TRUE, recursive=TRUE)
	csv_data_path <- grep("csv", all_data_folders, ignore.case=TRUE, value=TRUE)
	r_data_path <- file.path(subfolder_path, "R data")
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



#iterate through each main folder, then into each GFP/rabies folder to modify the dataset$image name
for (folder in folders) {
	#list GFP and rabies mcherry folders
	gfp_rabies_folders <- list.dirs(path=folder, full.names = TRUE, recursive= FALSE)
	for (subfolder in gfp_rabies_folders) {
		change_dataset_image_name(subfolder)
	}	
}


#iterate through each main folder, then into each GFP/rabies folder to compile data
compiled_data <- data.frame()
for (folder in folders) {
	#list GFP and rabies mcherry folders
	gfp_rabies_folders <- list.dirs(path=folder, full.names = TRUE, recursive= FALSE)
	for (subfolder in gfp_rabies_folders) {
		compiled_data <-load_and_compile_r_data_files(subfolder, compiled_data)
	}	
}

#iterate through each main folder, then into each GFP/rabies folder to resave csv files with current R data
for (folder in folders) {
	#list GFP and rabies mcherry folders
	gfp_rabies_folders <- list.dirs(path=folder, full.names = TRUE, recursive= FALSE)
	for (subfolder in gfp_rabies_folders) {
		resave_csv_file(subfolder)
	}	
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




#function to identify/remove starter cell population from rabies data set
#uses centroid x/y coordinates of image and area of each recorded cell, then calculates the Euclidean distance between mCherry+ and GFP+ cells in the same slice
#will output based on each unique sample identified
#multiply radi by scale factor (images are scaled down 1/4 to calculate pixel area)


#function to calculate euclidean distance
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
#and another function to display starter cell counts for each sample, as a table
remove_starter_cells <- function(df) {
	starter_cells <- identify_starter_cells(df)
	cleaned_compiled_data <- anti_join(df, starter_cells, by = c("animal", "x","y","area","ML","DV","acronym","name","image"))
	cleaned_rabies_data <- cleaned_compiled_data[!grepl("GFP", cleaned_compiled_data$image), ]
	return(cleaned_rabies_data)
}


count_starter_cells <- function(df) {
	starter_cells <- identify_starter_cells(df)
	sample_names <- get_sample_names(starter_cells)
	
	sample_name_idx <- paste0("_", sample_names, "_")
	
	starter_cell_counts <- data.frame(sample = character(), total_sc = integer(), stringsAsFactors = FALSE)
	for (i in 1:length(sample_names)) {
	
	total_sc <- sum(grepl(sample_name_idx[i], starter_cells$image))
	
	starter_cell_counts <- rbind(starter_cell_counts, data.frame(sample= sample_names[i], total_sc = total_sc, stringsAsFactors=FALSE))
	
	}	
	
	return(starter_cell_counts)
	print(starter_cell_counts)
}



##function to count starter cells for each sample within a given AP range
#run identify_starter_cells to generate df of all starter cells
count_starter_cells_APrange <- function(startercell_df,AP1,AP2) {
	sample_names <- get_sample_names(startercell_df)
	sample_name_idx <- paste0("_", sample_names, "_")
	
	min_AP <- min(c(AP1,AP2))
	max_AP <- max(c(AP1,AP2))
	
	starter_cell_counts <- data.frame(sample=character(), sc_count = integer(), stringsAsFactors = FALSE)
	
	for (i in 1:length(sample_names)) {
		sample_data <- startercell_df[grepl(sample_name_idx[i],startercell_df$image),]
		cell_count <- length(which(sample_data$AP<= max_AP & sample_data$AP>= min_AP))
		starter_cell_counts <- rbind(starter_cell_counts, data.frame(sample= sample_names[i], sc_count = cell_count, stringsAsFactors=FALSE))
		
	}
	return(starter_cell_counts)
	print(starter_cell_counts)
}







########################################### PLOTTING FUNCTIONS FOR RABIES DATA ########################################


setwd('/Users/snfreda/Desktop/Avp Rabies Data/Figures/Revisions')

library(ggplot2)

total_input_cells <- remove_starter_cells(compiled_data)
row.names(total_input_cells)<-c(1:nrow(total_input_cells))

#function to combine lists of cell counts and add column specifying sex
combine_data <- function(list, sex) { 
	do.call(rbind, lapply(list, function(df) {
		df$sex <- sex
		df
	}))	
}


############   bar plots of % total input of all samples, compared by sex #######################

plot_percent_total_input_by_sex <- function(total_input_cells) {
	
	split_data <- split_compiled_data(total_input_cells, "sex")
	female_data <- split_data$female
	male_data <- split_data$male
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline
	female_data <-female_data[!female_data$id %in% c(0,129,303,131,23,1015,276,284,411), ] #for avp
	male_data <- male_data[!male_data$id %in% c(0,129,303,131,23,1015,276,284,411), ] #for avp
	
	total_input_cells_by_sample_female <- split_compiled_data(female_data,"sample")
	total_input_cells_by_sample_male <- split_compiled_data(male_data,"sample")
	
	female_normalized_cell_counts <- list()       #cell counts normalized as % total input 
	male_normalized_cell_counts <- list()			# i.e. cells per parent region are divided by total # of input cells
	
	for (i in seq_along(total_input_cells_by_sample_female)) {
		data <- total_input_cells_by_sample_female[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		data <- data[!data$id %in% c(0,129,23,411),]  #remove any cells located in ventricle or outside registered brain outline
		total_cell_count <- nrow(data)
		regions <- data$name
		all_parent_regions <- unlist(lapply(regions, function(x)strsplit(x,",")[[1]][1]))
		unq_parent_regions <- unique(all_parent_regions)
		
		cell_counts = data.frame(parent_region=character(),cell_count = numeric(), normalized_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_parent_regions)) {
			count <- sum(grepl(unq_parent_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			cell_counts = rbind(cell_counts, data.frame(parent_region = unq_parent_regions[j],cell_count = count, normalized_cell_count = norm_count))
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		female_normalized_cell_counts[[sample_name]] <- sorted_cell_counts
		}
	
	for (k in seq_along(total_input_cells_by_sample_male)) {
		data <- total_input_cells_by_sample_male[[k]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		data <- data[!data$id %in% c(0,129,23,411),]  #remove any cells located in ventricle or outside registered brain outline
		total_cell_count <- nrow(data)
		regions <- data$name
		all_parent_regions <- unlist(lapply(regions, function(x)strsplit(x,",")[[1]][1]))
		unq_parent_regions <- unique(all_parent_regions)
		
		cell_counts = data.frame(parent_region=character(),cell_count = numeric(), normalized_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_parent_regions)) {
			count <- sum(grepl(unq_parent_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			cell_counts = rbind(cell_counts, data.frame(parent_region = unq_parent_regions[j],cell_count = count, normalized_cell_count = norm_count))}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		male_normalized_cell_counts[[sample_name]] <- sorted_cell_counts
		}
	
	#save(male_normalized_cell_counts, female_normalized_cell_counts, file='normalized_cell_counts.RData')
	
	############# plotting ####################
	
	#first combine data into single data frame with an additional column specifying sex
	
	female_data <- combine_data(female_normalized_cell_counts, "female")
	male_data <- combine_data(male_normalized_cell_counts, "male")
	
	combined_data <- rbind(female_data,male_data)
	combined_data_save <- combined_data
	#calculate SEM for each brain region within each sex
	combined_data <- combined_data %>% 
	group_by(parent_region,sex) %>%
	summarise(mean_norm_count = mean(normalized_cell_count), sem_norm_count = sd(normalized_cell_count)/sqrt(n()))
	
	#order brain regions based on total cell counts in females
	brain_region_order <- combined_data %>% filter(sex == 'female') %>% arrange(desc(mean_norm_count)) %>%
	pull(parent_region)
	
	combined_data$parent_region <- factor(combined_data$parent_region, levels = brain_region_order)
	
	#create box plots with SEM error bars
	p<-ggplot(combined_data, aes(x=parent_region, y = mean_norm_count, fill=sex)) + geom_boxplot(position = position_dodge(width=0.8)) + 
	geom_errorbar(aes(ymin= mean_norm_count - sem_norm_count, ymax= mean_norm_count + sem_norm_count), position= position_dodge(width=0.8), width=0.25) +
	labs(x = "Brain Region", y="Normalized Cell Count (% total input)", fill="Sex") + theme(axis.text.x = element_text(angle=90, hjust=1))	
	
	save(combined_data_save, file='combined_input_data.RData')
	print(p)
	
}


########### density contour plots ###################

library(wholebrain)
library(tidyverse)
#library(RColorBrewer)
library(viridis)

#first, you need to load two .R data files saved under "Documents/R" folder 
	#1. EPSatlas.RData
	#2. AtlasIndex.RData
	
load("/Users/snfreda/Documents/R/EPSatlas.RData")
load("/Users/snfreda/Documents/R/atlasIndex.RData")	
	

plot_density_contour_by_sex <- function(total_input_cells, sex) {
	
	opt_sex <- c("female","male")
	criterion <- match.arg(sex, choices=opt_sex)
	
	
	split_data <- split_compiled_data(total_input_cells, "sex")
	female_data <- split_data$female
	male_data <- split_data$male
	
	if (criterion == "female") {
		cell_data <- female_data
	}
	
	else if (criterion == "male") {
		cell_data <- male_data
		
	}
	
	#plot density contour for each unique AP coordinate/brain slice
	#plot will automatically adjust cell locations to be in the right hemisphere
	
	unq.APcoord<-unique(cell_data[,2])
	
	for (i in 1:length(unq.APcoord)) {
	
	idx<-which(cell_data$AP == unq.APcoord[i])
	dataset <- cell_data[idx,]
	dataset$right.hemisphere <- TRUE
	dataset$ML <- abs(dataset$ML)
	
	
	#information needed for wholebrain plotting of atlas outlines
	coordinate <- unq.APcoord[i]
	k <- which(abs(coordinate - atlasIndex$mm.from.bregma[1:132]) == min(abs(coordinate - atlasIndex$mm.from.bregma[1:132])))
    plate <- atlasIndex$plate.id[which(abs(coordinate - 
                  atlasIndex$mm.from.bregma) == min(abs(coordinate - 
                  atlasIndex$mm.from.bregma)))]
    scale = 0.9579832
	numPaths <- EPSatlas$plates[[k]]@summary@numPaths
	bregmaY = 200
	bregmaX = 5640
	xmin <- min(EPSatlas$plates[[k]][[1]]@paths$path@x) - 97440/2
	
	
	#initialize data for density contour plot, ML/DV values need to be converted back to "outline" coordinates
	newX = (dataset$ML *1000-bregmaY) * 97440/456/25 + 97440/2
	newY = (8210+ dataset$DV * 1000-bregmaY) *97440/456/25
	
	plot.data = data.frame(newX,newY)
	data.plot = ggplot(plot.data, aes(x=newX, y=newY))
	

	
	
	#initialize plot, then use ggplot to create density contours 
	#once ggplot is generated, overlay outlines ***warning: in newer MACs graph might not display immediately
	
	quartz()
	plot(c(0, 97400), c(0, 68234.56), ylim = c(0, 68234.56), 
        xlim = c(0, 97440), asp = 1, axes = F, xlab = "", ylab = "", 
        col = 0, main = paste("Bregma:", coordinate, 
            "mm"), font.main = 1)
		
	#option 1:	
	data.plot + geom_density_2d_filled(contour_var="ndensity") + theme(legend.position="none") + 
	coord_fixed(ratio=1,xlim=c(0,97440),ylim=c(0,68234.56))
	
	#option 2: if geom_density cuts off contours:
	
	#pad <- 2000
	#padded_data <- rbind(
  	#plot.data,
  	#data.frame(newX = c(-pad, max(plot.data$newX) + pad), 
             #newY = c(-pad, max(plot.data$newY) + pad))
	#)
	#data.plot = ggplot(padded_data, aes(x=newX, y=newY))
	
	#data.plot + geom_density_2d_filled(contour_var="ndensity") + theme(legend.position="none") + 
	#coord_fixed(ratio=1,xlim=c(0,97440),ylim=c(0,68234.56))
	 
	
	lapply(1:numPaths, function(N) {
             polygon(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                xmin, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "black", lwd = 1)
            polygon(-(EPSatlas$plates[[k]][[N]]@paths$path@x - 
                 xmin - 97440/2) + 97440/2, EPSatlas$plates[[k]][[N]]@paths$path@y, 
                border = "black", lwd = 1)
       })
	

	}
	

	
}

################################# SCHEMATIC PLOTS ############################

plot_schematic_sample <- function(total_input_cells, sample) {
	split_data <- split_compiled_data(total_input_cells, "sample")
	sample_opt = names(split_data)
	
	criterion <- match.arg(sample, choices=sample_opt)
	cell_data <- split_data[[criterion]]
	
	
	unq.APcoord <- unique(cell_data$AP)
	
	#### assign color for all somatic signals
	### all somas will be automatically plotted on the right hemisphere
	
	#cell_data$color <- '#8DC63F' #green color for male samples, oxt
	#cell_data$color <- '#005A88' #teal/blue color for female samples, oxt
	
	#cell_data$color <- '#330066' #dark purple for female samples, avp
	#cell_data$color <- '#F08000' #orange for male samples, avp
	
	#cell_data$color <- '#B72467' #magenta/purple for avp samples, all sexes
	cell_data$color <- '#00ADDC' #light blue color for oxt samples, all sexes
	
	
	
	for (j in 1:length(unq.APcoord)) {
		idx<-which(cell_data$AP == unq.APcoord[j])
		dataset <- cell_data[idx,]
		dataset$right.hemisphere <- TRUE
		dataset$ML <- abs(dataset$ML)
		schematic.plot(dataset, title=TRUE, scale.bar=TRUE, mm.grid=TRUE, pch=21, col=gray(0.1))
		
	}

}

##################### PLOTS MALE VS FEMALE, ALL SAMPLES ########################
###will plot all somatic signals for either male or female and assign a specific color to each sample
### all somas will be automatically plotted on the right hemisphere

plot_schematic_by_sex <- function(total_input_cells, sex) {
	opt_sex <- c("female","male")
	criterion <- match.arg(sex, choices=opt_sex)
	
	
	split_data <- split_compiled_data(total_input_cells, "sex")
	female_data <- split_data$female
	male_data <- split_data$male
	
	if (criterion == "female") {
		cell_data <- female_data
	}
	
	else if (criterion == "male") {
		cell_data <- male_data
		
	}


	unq.APcoord<-unique(cell_data[,2])
	
	
	
	unq.samples <- unique(cell_data$image)
	all.sample.names <- sapply(strsplit(unq.samples,"_"), "[[",2)
	unq.samples <- unique(all.sample.names)
	
	colors <- viridis(length(unq.samples),option="viridis")
	
	
	#assign colors of each sample in the cell_data variable
	
	
	for (i in seq_along(unq.samples)) {
		idx <- grepl(paste0("_",unq.samples[i],"_"),cell_data$image)
		cell_data$color[idx] <- colors[i]
	}
		
	
	for (k in 1:length(unq.APcoord)) {
	
	idx<-which(cell_data$AP == unq.APcoord[k])
	dataset <- cell_data[idx,]
	dataset$right.hemisphere <- TRUE
	dataset$ML <- abs(dataset$ML)
	schematic.plot(dataset, title=FALSE, scale.bar=TRUE, mm.grid=FALSE, pch=21, col=gray(0.1))
	}	
	
}

##################### PLOT HEAT MAP OF "SIMILARITY" ########################
#### will assign colors based on p values for each region, first will calculate p value based on t-test 
##then will color in each region with the assigned heat map color
#more "similar" regions will be colored with bright colors

#input A/P coordinate to plot a single heatmap






plot_heatmap_similarity_ttest <- function(df1, df2, APcoord) {
	#APcoord1 <- unique(df1$AP)
	#APcoord2 <- unique(df2$AP)
	
	#all_APcoord <- unique(c(APcoord1,APcoord2))
	
	#criterion <- match.arg(APcoord, choices=all_APcoord)
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline
	cell_data1 <-df1[!df1$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	cell_data2 <- df2[!df2$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	
	
	#normalize cell counts to total cell counted
	cell_data1_samples <- split_compiled_data(cell_data1,"sample")
	cell_data2_samples <- split_compiled_data(cell_data2,"sample")
	
	
	
	
	cell_data1_normalized_counts <- list()       #cell counts normalized as % total input 
	cell_data2_normalized_counts <- list()			# i.e. cells per parent region are divided by total # of input cells
	
	for (i in seq_along(cell_data1_samples)) {
		data <- cell_data1_samples[[i]]
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
		
		cell_data1_normalized_counts[[sample_name]] <- sorted_cell_counts
		}
	
	
	for (i in seq_along(cell_data2_samples)) {
		data <- cell_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		cell_counts = data.frame(brain_region=character(),structure_id= numeric(), cell_count = numeric(), normalized_cell_count = numeric(), sqrt_norm_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(grepl(unq_brain_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], cell_count = count, normalized_cell_count = norm_count,sqrt_norm_cell_count=sqrt_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		cell_data2_normalized_counts[[sample_name]] <- sorted_cell_counts
	}
	
	
	
	
	
	cell_data1 <- cell_data1_normalized_counts
	cell_data2 <- cell_data2_normalized_counts
	
	#######calculate p-values (comparing normalized cell counts in df1 vs df2) as the basis for the "similarity" heatmap values
	
	#1. organize data and calculate p-values
	
	p_values <- c()
	data1_num_samples <- length(cell_data1_samples)
	data2_num_samples <- length(cell_data2_samples)
	
	
	all_brain_regions <- unique(c(
	unlist(lapply(cell_data1, function(df) df$brain_region)),
	unlist(lapply(cell_data2, function(df) df$brain_region))
	))
	
	all_structure_ids <- unique(c(
	unlist(lapply(cell_data1, function(df) df$structure_id)),
	unlist(lapply(cell_data2, function(df) df$structure_id))
	))

	
	for (k in 1:length(all_brain_regions)) {
		data1_region <- data.frame(brain_region=rep(all_brain_regions[k], data1_num_samples),structure_id=rep(all_structure_ids[k],data1_num_samples), cell_count=rep(0, data1_num_samples), normalized_cell_count=rep(0,data1_num_samples))
		row.names(data1_region) <- names(cell_data1)
		
		data2_region <- data.frame(brain_region=rep(all_brain_regions[k], data2_num_samples),structure_id=rep(all_structure_ids[k],data2_num_samples), cell_count=rep(0, data2_num_samples), normalized_cell_count=rep(0,data2_num_samples))
		row.names(data2_region) <- names(cell_data2)
		
		for (s in 1:length(cell_data1)) {
			sample_data <- cell_data1[[s]]
			
			row_idx <- which(sample_data$brain_region == data1_region$brain_region[1])
			
			if (length(row_idx) >0) {
				data1_region$cell_count[s] <- sample_data$cell_count[row_idx]
				data1_region$normalized_cell_count[s] <- sample_data$normalized_cell_count[row_idx]
				}
		}
		
		for (r in 1:length(cell_data2)) {
			sample_data <- cell_data2[[r]]
			
			row_idx2 <- which(sample_data$brain_region == data2_region$brain_region[1])
			
			if (length(row_idx2) >0) {
				data2_region$cell_count[r] <- sample_data$cell_count[row_idx2]
				data2_region$normalized_cell_count[r] <- sample_data$normalized_cell_count[row_idx2]
				}
		}
	
	##### determine p-values for comparing each region ****IF TOTAL CELL COUNT IS LOW, THE AREA WILL BE DESREGARDED FROM T-TEST AND GIVEN ZERO VALUE
	if (mean(data1_region$normalized_cell_count)<1 & mean(data2_region$normalized_cell_count) < 0.2) {
		p_values[k] <- 0
		}
	else {
	test_result <- t.test(data1_region$normalized_cell_count, data2_region$normalized_cell_count)
	p_values[k] <- test_result$p.value
		}
					
	}

	plate_data_pvalues <- data.frame(p_value=p_values,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_pvalues$p_value = p_values
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	
	plate_data_pvalues$p_value[is.nan(plate_data_pvalues$p_value)] <- 0
	plate_data_pvalues$animal <- 'NA'
	
	
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
	g1<-ggplot(plate_data_pvalues, aes(x=animal, y=brain_region, fill=p_value)) + geom_tile() + scale_fill_viridis(option='B',name="p-value")+theme(legend.position="right")
	heatmap_data<-layer_data(g1)
	
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




plot_heatmap_similarity_adjustpval <- function(df1, df2, APcoord) {
	#APcoord1 <- unique(df1$AP)
	#APcoord2 <- unique(df2$AP)
	
	#all_APcoord <- unique(c(APcoord1,APcoord2))
	
	#criterion <- match.arg(APcoord, choices=all_APcoord)
	
	#remove potential outliers i.e. located in ventricle or outside registered brain outline
	cell_data1 <-df1[!df1$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	cell_data2 <- df2[!df2$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	
	
	#normalize cell counts to total cell counted
	cell_data1_samples <- split_compiled_data(cell_data1,"sample")
	cell_data2_samples <- split_compiled_data(cell_data2,"sample")
	
	
	
	
	cell_data1_normalized_counts <- list()       #cell counts normalized as % total input 
	cell_data2_normalized_counts <- list()			# i.e. cells per parent region are divided by total # of input cells
	
	for (i in seq_along(cell_data1_samples)) {
		data <- cell_data1_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		cell_counts = data.frame(brain_region=character(),structure_id = numeric(), cell_count = numeric(), normalized_cell_count = numeric(), sqrt_norm_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(data$id == unq_structure_id[j])
			norm_count <- count/total_cell_count * 100
			sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j],structure_id=unq_structure_id[j],cell_count = count, normalized_cell_count = norm_count,sqrt_norm_cell_count=sqrt_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		cell_data1_normalized_counts[[sample_name]] <- sorted_cell_counts
		}
	
	
	for (i in seq_along(cell_data2_samples)) {
		data <- cell_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		cell_counts = data.frame(brain_region=character(),structure_id= numeric(), cell_count = numeric(), normalized_cell_count = numeric(), sqrt_norm_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(data$id == unq_structure_id[j])
			norm_count <- count/total_cell_count * 100
			sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], cell_count = count, normalized_cell_count = norm_count,sqrt_norm_cell_count=sqrt_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		cell_data2_normalized_counts[[sample_name]] <- sorted_cell_counts
	}
	
	cell_data1 <- cell_data1_normalized_counts
	cell_data2 <- cell_data2_normalized_counts
	
	#######calculate p-value using Mann-Whitney U (Wilcox rank sum)
	
	#1. organize data and calculate p-values
	
	p_values <- c()
	data1_num_samples <- length(cell_data1_samples)
	data2_num_samples <- length(cell_data2_samples)
	
	
	all_brain_regions <- unique(c(
	unlist(lapply(cell_data1, function(df) df$brain_region)),
	unlist(lapply(cell_data2, function(df) df$brain_region))
	))
	
	all_structure_ids <- unique(c(
	unlist(lapply(cell_data1, function(df) df$structure_id)),
	unlist(lapply(cell_data2, function(df) df$structure_id))
	))
	
	all_parent_regions <- unique(sapply(strsplit(all_brain_regions, ","),'[',1))
	
	for (k in 1:length(all_brain_regions)) {
		data1_region <- data.frame(brain_region=rep(all_brain_regions[k], data1_num_samples),structure_id=rep(all_structure_ids[k],data1_num_samples), cell_count=rep(0, data1_num_samples), normalized_cell_count=rep(0,data1_num_samples))
		row.names(data1_region) <- names(cell_data1)
		
		data2_region <- data.frame(brain_region=rep(all_brain_regions[k], data2_num_samples),structure_id=rep(all_structure_ids[k],data2_num_samples), cell_count=rep(0, data2_num_samples), normalized_cell_count=rep(0,data2_num_samples))
		row.names(data2_region) <- names(cell_data2)
		
		for (s in 1:length(cell_data1)) {
			sample_data <- cell_data1[[s]]
			
			row_idx <- which(sample_data$brain_region == data1_region$brain_region[1])
			
			if (length(row_idx) >0) {
				data1_region$cell_count[s] <- sample_data$cell_count[row_idx]
				data1_region$normalized_cell_count[s] <- sample_data$normalized_cell_count[row_idx]
				}
		}
		
		for (r in 1:length(cell_data2)) {
			sample_data <- cell_data2[[r]]
			
			row_idx2 <- which(sample_data$brain_region == data2_region$brain_region[1])
			
			if (length(row_idx2) >0) {
				data2_region$cell_count[r] <- sample_data$cell_count[row_idx2]
				data2_region$normalized_cell_count[r] <- sample_data$normalized_cell_count[row_idx2]
				}
		}
	
	##### determine p-values for comparing each region ****IF TOTAL CELL COUNT IS LOW, THE AREA WILL BE DESREGARDED FROM T-TEST AND GIVEN ZERO VALUE
	if (mean(data1_region$normalized_cell_count)<1 & mean(data2_region$normalized_cell_count) < 0.2) {
		p_values[k] <- 0
		}
	else {
	test_result <- wilcox.test(data1_region$normalized_cell_count, data2_region$normalized_cell_count)
	p_values[k] <- test_result$p.value
		}
					
	}
	
	##apply FDR to all brain regions
	
	adj_p_values = mutate(p_adj= p.adjust(p_values, method="fdr")) ## Benjamini–Hochberg FDR

	plate_data_pvalues <- data.frame(p_value=p_values,structure_id=all_structure_ids,brain_region=all_brain_regions)
	plate_data_pvalues$p_value = p_values
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	
	plate_data_pvalues$p_value[is.nan(plate_data_pvalues$p_value)] <- 0
	plate_data_pvalues$animal <- 'NA'	
	
	
	
	
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
	g1<-ggplot(plate_data_ksvalues, aes(x=animal, y=brain_region, fill=ks_value)) + geom_tile() + scale_fill_viridis(option='B')
	heatmap_data<-layer_data(g1)
	
	for (t in 1:nrow(plateData)) {
	index=which(plateData$structure_id[t]==plate_data_pvalues$structure_id)
	if (length(index) >0){
	plateData$color[t]<-heatmap_data$fill[index]
		}
	}
	
	
	idx.ft=which(plateData$structure_id==1009)
	plateData$color[idx.ft]<-'#3b3b3b'
	
	
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



######################################
########### EXTRA MISC. CODE ##########################

plot_heatmap_similiarity_parentregion(df1,df2,APcoord) {

	cell_data1 <-df1[!df1$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	cell_data2 <- df2[!df2$id %in% c(0,129,303,131,23,1015,276,284,411), ]
	
	
#normalize cell counts to total cell counted
	cell_data1_samples <- split_compiled_data(cell_data1,"sample")
	cell_data2_samples <- split_compiled_data(cell_data2,"sample")
	
	
	
	
	cell_data1_normalized_counts <- list()       #cell counts normalized as % total input 
	cell_data2_normalized_counts <- list()			# i.e. cells per parent region are divided by total # of input cells
	
	for (i in seq_along(cell_data1_samples)) {
		data <- cell_data1_samples[[i]]
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
		
		cell_data1_normalized_counts[[sample_name]] <- sorted_cell_counts
		}
	
	
	for (i in seq_along(cell_data2_samples)) {
		data <- cell_data2_samples[[i]]
		sample_name <- strsplit(data$image[1],"_")[[1]][2]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		unq_structure_id <- unique(data$id)
		
		cell_counts = data.frame(brain_region=character(),structure_id= numeric(), cell_count = numeric(), normalized_cell_count = numeric(), sqrt_norm_cell_count = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(grepl(unq_brain_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j], structure_id =unq_structure_id[j], cell_count = count, normalized_cell_count = norm_count,sqrt_norm_cell_count=sqrt_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		cell_data2_normalized_counts[[sample_name]] <- sorted_cell_counts
	}
	
	cell_data1 <- cell_data1_normalized_counts
	cell_data2 <- cell_data2_normalized_counts
	

	
	#1. organize data and calculate p-values
	
	p_values <- c()
	data1_num_samples <- length(cell_data1_samples)
	data2_num_samples <- length(cell_data2_samples)
	
	
	all_brain_regions <- unique(c(
	unlist(lapply(cell_data1, function(df) df$brain_region)),
	unlist(lapply(cell_data2, function(df) df$brain_region))
	))
	
	all_structure_ids <- unique(c(
	unlist(lapply(cell_data1, function(df) df$structure_id)),
	unlist(lapply(cell_data2, function(df) df$structure_id))
	))




all_parent_regions <- unique(sapply(strsplit(all_brain_regions, ","),'[',1))
	
	for (k in 1:length(all_parent_regions)) {
		data1_region <- data.frame(brain_region=rep(all_parent_regions[k], data1_num_samples), cell_count=rep(0, data1_num_samples), normalized_cell_count=rep(0,data1_num_samples))
		row.names(data1_region) <- names(cell_data1)
		
		data2_region <- data.frame(brain_region=rep(all_parent_regions[k], data2_num_samples), cell_count=rep(0, data2_num_samples), normalized_cell_count=rep(0,data2_num_samples))
		row.names(data2_region) <- names(cell_data2)
		
		for (s in 1:length(cell_data1)) {
			sample_data <- cell_data1[[s]]
			
			row_idx <- grep(all_parent_regions[k],sample_data$brain_region)
			
			if (length(row_idx) >0) {
				data1_region$cell_count[s] <- sum(sample_data$cell_count[row_idx])
				data1_region$normalized_cell_count[s] <- sum(sample_data$cell_count[row_idx])/sum(sample_data$cell_count)*100
				}
		}
		
		for (r in 1:length(cell_data2)) {
			sample_data <- cell_data2[[r]]
			
			row_idx2 <- grep(all_parent_regions[k],sample_data$brain_region)
			
			if (length(row_idx2) >0) {
				data2_region$cell_count[r] <- sum(sample_data$cell_count[row_idx2])
				data2_region$normalized_cell_count[r] <- sum(sample_data$cell_count[row_idx2])/sum(sample_data$cell_count)*100
				}
				
		}		
		if (mean(data1_region$normalized_cell_count)<0.3 & mean(data2_region$normalized_cell_count) < 0.3) {
		p_values[k] <- -0.1
			}
		else {
		test_result <- t.test(data1_region$normalized_cell_count, data2_region$normalized_cell_count)
		p_values[k] <- test_result$p.value
			}
							

	}
	
	plate_data_pvalues <- data.frame(structure_id=all_structure_ids,brain_region=all_brain_regions)
	#plate_data_pvalues$p_value = p_values
	plate_data_pvalues$structure_id = all_structure_ids
	plate_data_pvalues$brain_region = all_brain_regions
	plate_data_pvalues$p_value=-0.1
	
	plate_data_pvalues$p_value[is.nan(plate_data_pvalues$p_value)] <- -0.1
	plate_data_pvalues$animal <- 'NA'
	
	#will assign p-value from "parent region" to each subregion within
	#i.e. paraventricular parvicellular, mangocellular, etc. subregions will be assigned p-value calculated for all of PVN
	for (r in 1:length(all_parent_regions)) {
		region_idx = grep(all_parent_regions[r],all_brain_regions)
		plate_data_pvalues$p_value[region_idx] = p_values[r]
	}
	
	
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
	g1<-ggplot(plate_data_pvalues, aes(x=animal, y=structure_id, fill=p_value)) + geom_tile() + scale_fill_viridis(option='B')
	heatmap_data<-layer_data(g1)
	
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




#### generate table of summary data oxt vs. avp #############

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

#function to add sample identifier to each row of dataframe
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

oxt_inputs <- total_input_cells[!total_input_cells$id %in% c(0,129,131,411), ] #delete outliers
avp_inputs <- total_input_cells[!total_input_cells$id %in% c(0,129,303,131,23,1015,276,284,411), ]


generate_retrograde_table <- function(oxt_inputs, avp_inputs) {
	
	#re-name brain regions to be larger category (i.e.  Paraventricular hypothalamic nucleus, magnocellular division --> Paraventricular hypothalamic nucleus)
	oxt_inputs$name <- split_relevant_brain_regions(oxt_inputs$name)
 	avp_inputs$name <- split_relevant_brain_regions(avp_inputs$name)
 	
 	#add column for sample identifier 
 	oxt_inputs <- add_sample_name(oxt_inputs)
 	avp_inputs <- add_sample_name(avp_inputs)
 	
 	#add column for sex identifier
 	oxt_inputs <- add_sex_identifier(oxt_inputs)
 	avp_inputs <- add_sex_identifier(avp_inputs)
	
	
	data_all <- bind_rows(oxt_inputs, avp_inputs)
	all_regions <- unique(data_all$brain_region)
	
	oxt_samples <- unique(oxt_inputs$animal)
	oxt_data_samples <- split_compiled_data(oxt_inputs, "sample")
	
	normalized_oxt_inputs <- list()
	
	#organize data based on samples
	for (i in seq_along(oxt_data_samples)) {
		data <- oxt_data_samples[[i]]
		sample_name <- data$animal[1]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		
		cell_counts = data.frame(brain_region=character(), cell_count = numeric(), normalized_inputs = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(grepl(unq_brain_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			#sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j],cell_count = count, normalized_inputs = norm_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		normalized_oxt_inputs[[sample_name]] <- sorted_cell_counts
		}
	
	combined_data_oxt <- data.frame()  
  
  for (i in seq_along(normalized_oxt_inputs)) {
    df <- normalized_oxt_inputs[[i]]
    df$sample <- names(normalized_oxt_inputs)[i]  # Add sample name
    combined_data_oxt <- rbind(combined_data_oxt, df)  # Combine with the existing data
  }
 	
 	# reshape the data frame so that each sample has a column
	df_wide_oxt <- reshape(combined_data_oxt, idvar = c("brain_region"), timevar = "sample", direction = "wide")
	df_wide_oxt[is.na(df_wide_oxt)] <- 0
	
	
	# calculate the average norm_pixels, total pixel counts and SEM for each brain_region across samples
	norm_counts_avg_oxt <- rowMeans(df_wide_oxt[, grepl("normalized", names(df_wide_oxt))], na.rm = TRUE)  
	total_cell_count <- rowSums(df_wide_oxt[, grepl("cell_count", names(df_wide_oxt))], na.rm = TRUE)
	
	norm_avg_oxt_female <- rowMeans(df_wide_oxt[, grepl("normalized_inputs.F", names(df_wide_oxt))], na.rm = TRUE)
	norm_avg_oxt_male <- rowMeans(df_wide_oxt[, grepl("normalized_inputs.M", names(df_wide_oxt))], na.rm = TRUE)
	norm_avg_diff <- norm_avg_oxt_female - norm_avg_oxt_male
	
	n_samples_oxt <- length(oxt_samples)
	sd_per_row <- apply(df_wide_oxt[, grepl("normalized", names(df_wide_oxt))], 1, sd, na.rm = TRUE)
	sem_input_oxt <- sd_per_row / sqrt(n_samples_oxt)
	
	oxt_summary_data <- data.frame(brain_region= df_wide_oxt$brain_region, total_cell_count=total_cell_count, norm_input = norm_counts_avg_oxt, sem_norm_input = sem_input_oxt, norm_avg_diff_sex = norm_avg_diff)
	
	
	
	######### organize vasopressin data in the same way ###########################
	
	avp_samples <- unique(avp_inputs$animal)
	avp_data_samples <- split_compiled_data(avp_inputs, "sample")
	
	normalized_avp_inputs <- list()
	
	#organize data based on samples
	for (i in seq_along(avp_data_samples)) {
		data <- avp_data_samples[[i]]
		sample_name <- data$animal[1]
		total_cell_count <- nrow(data)
		unq_brain_regions <- unique(data$name)
		
		cell_counts = data.frame(brain_region=character(), cell_count = numeric(), normalized_inputs = numeric(), stringsAsFactors=FALSE)
		for (j in 1:length(unq_brain_regions)) {
			count <- sum(grepl(unq_brain_regions[j], data$name))
			norm_count <- count/total_cell_count * 100
			#sqrt_count <- sqrt(count/total_cell_count)
			cell_counts <- rbind(cell_counts, data.frame(brain_region = unq_brain_regions[j],cell_count = count, normalized_inputs = norm_count))
			
		}
		
		sorted_cell_counts <- cell_counts %>% arrange(desc(cell_count))
		
		normalized_avp_inputs[[sample_name]] <- sorted_cell_counts
		}
	
	combined_data_avp <- data.frame()  
	
	 for (i in seq_along(normalized_avp_inputs)) {
    df <- normalized_avp_inputs[[i]]
    df$sample <- names(normalized_avp_inputs)[i]  # Add sample name
    combined_data_avp <- rbind(combined_data_avp, df)  # Combine with the existing data
  }
 	
 	# reshape the data frame so that each sample has a column
	df_wide_avp <- reshape(combined_data_avp, idvar = c("brain_region"), timevar = "sample", direction = "wide")
	df_wide_avp[is.na(df_wide_avp)] <- 0
	
	
	# calculate the average norm_pixels, total pixel counts and SEM for each brain_region across samples
	norm_counts_avg_avp <- rowMeans(df_wide_avp[, grepl("normalized", names(df_wide_avp))], na.rm = TRUE)  
	total_cell_count <- rowSums(df_wide_avp[, grepl("cell_count", names(df_wide_avp))], na.rm = TRUE)
	
	norm_avg_avp_female <- rowMeans(df_wide_avp[, grepl("normalized_inputs.F", names(df_wide_avp))], na.rm = TRUE)
	norm_avg_avp_male <- rowMeans(df_wide_avp[, grepl("normalized_inputs.M", names(df_wide_avp))], na.rm = TRUE)
	norm_avg_diff <- norm_avg_avp_female - norm_avg_avp_male
	
	n_samples_avp <- length(avp_samples)
	sd_per_row <- apply(df_wide_avp[, grepl("normalized", names(df_wide_avp))], 1, sd, na.rm = TRUE)
	sem_input_avp <- sd_per_row / sqrt(n_samples_avp)
	
	avp_summary_data <- data.frame(brain_region= df_wide_avp$brain_region, total_cell_count=total_cell_count, norm_input = norm_counts_avg_avp, sem_norm_input = sem_input_avp, norm_avg_diff_sex = norm_avg_diff)

	
	
################# combine oxt and avp data into single data frame, organize by brain region in alphabetical order ######################

	combined_summary_data <- full_join(oxt_summary_data, avp_summary_data,
                              by = "brain_region",
                              suffix = c(".oxt", ".avp"))
     
     #rearrange to be alphabetical
     combined_summary_data <- combined_summary_data %>%
  arrange(brain_region)                        
	
	
	
	return(combined_summary_data)
	


}
