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

setwd(args[1])

pc_num=20
min_genes=500
res=0.6
project=args[2]


seurat_mat <- readRDS(args[3])
#seurat_mat[["RNA"]] <- JoinLayers(seurat_mat[["RNA"]])

DefaultAssay(object = seurat_mat) <- "RNA"

#seurat_mat.split <- SplitObject(seurat_mat, split.by = "Sample_name") 
seurat_mat.split <- SplitObject(seurat_mat, split.by = "orig.ident") 

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)

run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
	# Print sample number
	#print(paste0("Sample ", unique(seu_sample_subset[['Sample_name']]), '...........')) 
	print(paste0("Sample ", unique(seu_sample_subset[['orig.ident']]), '...........')) 

	if(is.null(multiplet_rate)){
		print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
		# 10X multiplet rates table
      		#https://rpubs.com/kenneditodd/doublet_finder_example
      		multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.002, 0.004, 0.008, 0.012, 0.0160, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04, 0.044, 0.048, 0.052, 0.056, 0.06, 0.064, 0.068, 0.072, 0.076, 0.08),
					'Loaded_cells' = c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050, 14500, 15950, 17400, 18850, 20300, 21750, 23200, 24650, 26100, 27550, 29000),
					'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000))
          
          	print(multiplet_rates_10x)

          	multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
		  	dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
		  	dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
	      	print(paste('Setting multiplet rate to', multiplet_rate))
	        }
    
	# Pre-process seurat object with standard seurat workflow --- 
	sample <- NormalizeData(seu_sample_subset)
	sample <- FindVariableFeatures(sample)
	sample <- ScaleData(sample)
	sample <- RunPCA(sample, nfeatures.print = 10)
  
	# Finish pre-processing with pc_num
	sample <- RunUMAP(sample, dims = 1:pc_num)
	sample <- FindNeighbors(object = sample, dims = 1:pc_num)
	sample <- FindClusters(object = sample, resolution = res)

	# pK identification (no ground-truth) 
	#introduces artificial doublets in varying props, merges with real data set and 
	# preprocesses the data + calculates the prop of artficial neighrest neighbours, 
	# provides a list of the proportion of artificial nearest neighbours for varying
	# combinations of the pN and pK
	sweep_list <- paramSweep(sample, PCs = 1:pc_num, sct = FALSE)   
	sweep_stats <- summarizeSweep(sweep_list)
	bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
	# Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
	optimal.pk <- bcmvn %>% 
		dplyr::filter(BCmetric == max(BCmetric)) %>%
		dplyr::select(pK)
	optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
		        
        ## Homotypic doublet proportion estimate
        annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
        homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
	  
	nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
	nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
    
    	# run DoubletFinder
	sample <- doubletFinder(seu = sample, 
				PCs = 1:pc_num, 
	                        pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
	                        nExp = nExp.poi.adj) # number of expected real doublets
	# change name of metadata column with Singlet/Doublet information
	colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
      
	# Subset and save
	# head(sample@meta.data['doublet_finder'])
      	# singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
      	# singlets$ident
      	double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
      	double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
        return(double_finder_res)
}

seurat_mat.split <- lapply(seurat_mat.split, run_doubletfinder_custom)

doublet_finder_res <- data.frame(bind_rows(seurat_mat.split))
rownames(doublet_finder_res) <- doublet_finder_res$row_names
doublet_finder_res$row_names <- NULL
seurat_mat <- AddMetaData(seurat_mat, doublet_finder_res, col.name = 'doublet_finder')

seurat_mat <- subset(seurat_mat, subset = doublet_finder == "Singlet") 

saveRDS(seurat_mat, file = args[4])

q() 
