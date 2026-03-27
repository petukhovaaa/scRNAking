#!/usr/bin/Rscript --vanilla

library(Matrix)
library(devtools)
library(Seurat)
library(dplyr)
library(grid)
library(gridExtra)
library(tibble)
library(ggplot2)
library(hdf5r)
library(Signac)
library(tidyverse)
library(ape)
library(dendextend)
library(dendsort)

#library(SeuratDisk)

library(SeuratObject)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
	          print("The arguments are needed!")
  stop("Requires command line argument.")
}

setwd(args[1])
getwd()

pc_num=20
min_genes=500
res=0.6
project=args[2]
dir_seurat_objects=args[3]

print(project)
print(dir_seurat_objects)

dir.create(dir_seurat_objects)

# reading input

folder_list <- list.dirs(path = getwd(), full.names = TRUE, recursive = FALSE)
folder_list <- folder_list[grepl("Samineni_SR", basename(folder_list))]
print(folder_list)

for (main_folder in folder_list) {
	sample_folders <- list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
	sample_folders <- sample_folders[grepl("SR", basename(sample_folders))]
	print(sample_folders)
	for (sample_path in sample_folders) {
		h5_path <- file.path(sample_path, "outs", "filtered_feature_bc_matrix.h5")
		if (!file.exists(h5_path)) {
	            message("File not found: ", h5_path)
		next
	        }
		sample_name <- basename(sample_path)
		print(c(main_folder,h5_path,sample_name))
		
		counts <- Read10X_h5(h5_path)
		seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
		seurat_obj_name = paste0('seurat_obj_',sample_name)
		assign(seurat_obj_name, seurat_obj, envir = .GlobalEnv)
		
	        }
}

seurat_list <- mget(ls(pattern = "seurat_obj_SR")) 

seurat_mat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = project)

saveRDS(seurat_mat, paste0(dir_seurat_objects,'/',args[4]))

q() 
