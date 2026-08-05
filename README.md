# Oxt_Avp_mapping
# Overview

Code accompanying the manuscript: “Comparative atlas of PVN oxytocin and vasopressin systems reveals sex-conserved input-output wiring”

Authors: Sara N. Freda, Yasmeen F. Lowe, Michael F. Priest, Deanna Badong, Yuejun Liu, Lei Xiao, Yevgenia Kozorovitskiy	 		 		
	 	 		
This repository utilizes a modified version of the open-source R package wholebrain, which performs registration of mouse brain sections to the Allen Common Coordinate Framework. The original package is available at: https://github.com/tractatus/wholebrain.

# Repository Structure

## R Code files:

01_segment_register_rabiesimages.R
Segments somatic signals from rabies tracing and registers coronal brain sections to the Allen Common Coordinate Framework using the wholebrain package. Contains a few plotting functions to visualize detected somatic signals.

02_segment_register_axonimages.R
Segments axonal projections and registers coronal brain sections to the Allen Common Coordinate Framework using a modified version of the wholebrain package. Contains a few plotting functions to visualize detected axons.

03_quantify_axon_pixels.R
Script to modify wholebrain package to register detected axon pixels. Needed to run ‘02_segment_register_axonimages.R’ code.

04_analyze_rabies_data.R
Script to compile, analyze, and plot data from registered brain slices of axonal projection tracing.

05_analyze_projection_data.R
Script to compile, analyze, and plot data from registered brain slices of rabies tracing.

## Quantification files:

.RData files of segmented cells and axons for each subject analyzed
The repository contains subfolders for each experiment type (e.g., retrograde or anterograde tracing) for each system (oxytocin or vasopressin). Each subject has a single folder of all analyzed files (.RData). Folders are named after the specimen number and contains either F or M string to specify whether the specimen is from "female" or "male" cohort. 
