#!/usr/bin/Rscript --vanilla

options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages('harmony')

library(harmony)
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
library(SeuratObject)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
	  print("The arguments are needed!")
  stop("Requires command line argument.")
}

setwd(args[1])

pc_num=50
min_genes=500
res=1
project=args[2]

library(future)
plan("multisession", workers = 3)
options(future.globals.maxSize = 350 * 1024^3)

seurat_mat <- readRDS(args[3])

#seurat_mat[["RNA"]] <- split(seurat_mat[["RNA"]], f = seurat_mat$orig.ident)
seurat_mat[["RNA"]] <- JoinLayers(seurat_mat[["RNA"]])
#seurat_mat[["RNA"]] <- split(seurat_mat[["RNA"]], f = seurat_mat$Sample_name)
#seurat_mat[["RNA"]] <- split(seurat_mat[["RNA"]], f = seurat_mat$batch)
seurat_mat[["RNA"]] <- split(seurat_mat[["RNA"]], f = seurat_mat$condition)

#seurat_mat  <- IntegrateLayers(object = seurat_mat, assay = 'RNA', method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = T, k.weight = 1, dims = 1:2)

seurat_mat  <- IntegrateLayers(object = seurat_mat, assay = 'RNA', method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = T, k.weight = 30)

# re-join layers after integration
#seurat_mat[["RNA"]] <- JoinLayers(seurat_mat[["RNA"]])

seurat_mat  <- FindNeighbors(seurat_mat , reduction = "harmony", dims = 1:pc_num)
seurat_mat  <- FindClusters(seurat_mat , resolution = res, cluster.name = "harmony_clusters")

seurat_mat  <- RunUMAP(seurat_mat , reduction = "harmony", dims = 1:pc_num,reduction.name = 'umap.integrated')

seurat_mat = AddMetaData(seurat_mat,Embeddings(seurat_mat[["umap.integrated"]]),colnames(Embeddings(seurat_mat[["umap.integrated"]])))

#DimPlot(seurat_mat , reduction = "umap.integrated",label = T) + DimPlot(seurat_mat , reduction = "umap.integrated",label = F, group.by = 'orig.ident')

pdf(paste("Umap_clusters+libs_",project,".pdf"),heigh=15,width=30)
DimPlot(seurat_mat, reduction = "umap",label = T,raster = F, group.by = 'harmony_clusters') + DimPlot(seurat_mat, reduction = "umap",label = F,raster = F, group.by = 'orig.ident')
dev.off()

pdf(paste("Umap_clusters+libs_integrated_",project,".pdf"),heigh=15,width=30)
DimPlot(seurat_mat, reduction = "umap.integrated",label = T,raster = F, group.by = 'harmony_clusters') + DimPlot(seurat_mat, reduction = "umap.integrated",label = F,raster = F, group.by = 'orig.ident')
dev.off()

saveRDS(seurat_mat, file = args[4])

q() 
