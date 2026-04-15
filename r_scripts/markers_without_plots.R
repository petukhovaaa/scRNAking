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

#library(qqconf)
#library(metap)
#library(reticulate)

#library(SeuratDisk)
#library(BPCells)

library(SeuratObject)
#library(Azimuth)
#library(MuDataSeurat)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
	  print("The arguments are needed!")
  stop("Requires command line argument.")
}

library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize = 350 * 1024^3)

setwd(args[1])

project=args[2]

seurat_mat <- readRDS(args[3])
seurat_mat_idents <- args[4]

### finding markers in defined groups

Idents(seurat_mat) <- seurat_mat_idents

seurat_mat_subset <- subset(seurat_mat, downsample = 250)

seurat_mat_subset[["RNA"]] <- JoinLayers(seurat_mat_subset[["RNA"]])

seurat_mat_markers <- FindAllMarkers(seurat_mat_subset,
				     assay = "RNA",
				     only.pos = TRUE,
				     verbose = TRUE)

write.table(seurat_mat_markers, paste('all_markers_',project,'.csv',sep=''), sep = ',')

seurat_mat_markers_025 <- subset(seurat_mat_markers, subset = pct.1 > 0.25) # defining reasonabe cutoff for percentage of the cell expressing a gene to consider it to be a marker

write.csv(seurat_mat_markers_025,paste('all_markers_025_',project,'.csv', sep = ''))

q() 
