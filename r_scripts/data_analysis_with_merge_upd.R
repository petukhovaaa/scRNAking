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
dir_seurat_objects_raw=args[3]
dir_seurat_objects=args[4]

print(project)
print(dir_seurat_objects_raw)
print(dir_seurat_objects)

dir.create(dir_seurat_objects)
dir.create(dir_seurat_objects_raw)

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
		#seurat_obj_name = paste0('seurat_obj_',sample_name)
		#assign(seurat_obj_name, seurat_obj, envir = .GlobalEnv)
		# Save RDS
		saveRDS(seurat_obj, file = file.path(dir_seurat_objects_raw, paste0("seurat_mat_", sample_name, "_before_qc.RDS")))
	        }
}

#seurat_list <- mget(ls(pattern = "seurat_obj_SR")) 

files_list <- list.files(path = dir_seurat_objects_raw, full.names = TRUE, recursive = FALSE)
seurat_list <- lapply(files_list, readRDS)
names(seurat_list) <- sapply(files_list, function(x) tools::file_path_sans_ext(basename(x)))

seurat_mat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = project)

saveRDS(seurat_mat, paste0(dir_seurat_objects,'/',args[5]))

DefaultAssay(object = seurat_mat) <- "RNA"

seurat_mat[["percent.mt"]] <- PercentageFeatureSet(seurat_mat, pattern = "^mt-")
seurat_mat <- NormalizeData(seurat_mat)

seurat_mat <- subset(seurat_mat, subset = (nFeature_RNA < 15000) & (nFeature_RNA > min_genes))
seurat_mat <- subset(seurat_mat, subset = percent.mt < 5)

seurat_mat <- NormalizeData(seurat_mat)

seurat_mat <- FindVariableFeatures(seurat_mat, assay = 'RNA')

# scale data
seurat_mat <- ScaleData(object = seurat_mat)

# PCA
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(object = seurat_mat),verbose=T)
Embeddings(seurat_mat, reduction="pca") %>% is.na() %>% sum

library(future)
plan("multisession", workers = 3)
options(future.globals.maxSize = 350 * 1024^3)

# dimension reduction and cluster
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:pc_num)
seurat_mat <- FindClusters(seurat_mat, resolution = res)
print(levels(Idents(seurat_mat)))

# visualization
seurat_mat <- RunTSNE(seurat_mat, dims = 1:pc_num)
seurat_mat <- RunUMAP(seurat_mat, dims = 1:pc_num)

seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["tsne"]]),colnames(Embeddings(seurat_mat[["tsne"]])))
seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["umap"]]),colnames(Embeddings(seurat_mat[["umap"]])))

pdf(paste("Umap_clusters+libs_",project,".pdf",sep=''),heigh=15,width=30)
DimPlot(seurat_mat, reduction = "umap",label = T,raster = F) + DimPlot(seurat_mat, reduction = "umap",label = F,raster = F, group.by = 'orig.ident')
dev.off()

pdf(paste("VlnPlot_nFeature_nCounts_mt_libs_",project,".pdf",sep=''),heigh=15,width=45)
VlnPlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), pt.size = 0, group.by = 'orig.ident')
dev.off()

pdf(paste("Featureplot_nFeature_nCounts_mt_",project,".pdf",sep=''),heigh=15,width=15)
FeaturePlot(seurat_mat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),reduction='umap',coord.fixed=T,raster = F)
dev.off()

saveRDS(seurat_mat, file = paste0(dir_seurat_objects,'/',args[6]))

q() 
