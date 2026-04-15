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

library(future)
plan("multisession", workers = 3)
options(future.globals.maxSize = 350 * 1024^3)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
	          print("The arguments are needed!")
  stop("Requires command line argument.")
}

setwd(args[1])
getwd()

project=args[2]
print(args[3])
print(args[4])
seurat_mat_1=readRDS(args[3])
seurat_mat_2=readRDS(args[4])

#seurat_mat_1[["RNA"]] <- JoinLayers(seurat_mat_1[["RNA"]])
#seurat_mat_2[["RNA"]] <- JoinLayers(seurat_mat_2[["RNA"]])

seurat_mat_merged <- merge(x = seurat_mat_1,
			   y = seurat_mat_2,
			   merge.data = TRUE,
			   merge.dr = FALSE,
			   project = project
				       )

seurat_mat_merged[['RNA']] <- JoinLayers(seurat_mat_merged[['RNA']])

seurat_mat_merged[["RNA"]] <- split(seurat_mat_merged[["RNA"]], f = seurat_mat_merged$orig.ident)

saveRDS(seurat_mat_merged,args[5])
