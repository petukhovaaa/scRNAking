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
seurat_mat=readRDS(args[3])

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

saveRDS(seurat_mat, file = args[4])

q() 
