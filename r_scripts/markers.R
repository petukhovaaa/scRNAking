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

seurat_mat_markers_035 <- subset(seurat_mat_markers, subset = pct.1 > 0.35) # defining reasonabe cutoff for percentage of the cell expressing a gene to consider it to be a marker

top3_markers = as.data.frame(seurat_mat_markers_035 %>%
				       group_by(cluster) %>% 
				       top_n(n = 3, wt = c(avg_log2FC)))

top5_markers = as.data.frame(seurat_mat_markers_035 %>% 
				       group_by(cluster) %>% 
				       top_n(n = 5, wt = avg_log2FC))

write.csv(top3_markers,paste('top3_markers_',project,'.csv', sep = ''))

write.csv(top5_markers,paste('top5_markers_',project,'.csv', sep = ''))

##### markers dotplots

pdf(file = paste('DotPlot_',project,'_top3_markers.pdf',sep=''), width = 30, height = 15) 
Seurat::DotPlot(seurat_mat_subset, assay = 'RNA', features = unique(top3_markers$gene)) +
	  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
							     vjust = 1,
							     size = 15, 
							     hjust = 1)) 
dev.off()

pdf(file = paste('DotPlot_',project,'_top5_markers.pdf',sep = ''), width = 40, height = 15)                                
Seurat::DotPlot(seurat_mat_subset, assay = 'RNA', features = unique(top5_markers$gene)) +
	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
							   vjust = 1,
							   size = 15,                                                   
							   hjust = 1)) 
dev.off()

#### markers heatmaps

seurat_mat_heat <- ScaleData(object = seurat_mat_subset, assay = 'RNA', features = unique(top3_markers$gene))

pdf(file = paste('heatmap_top3_markers_',project,'.pdf',sep = ''), width = 30, height = 20)
DoHeatmap(seurat_mat_heat, assay = 'RNA', features = unique(top3_markers$gene), size = 7, lines.width = 7, raster= F) + NoLegend() + theme(axis.text.y = element_text(size = 11)) + theme(plot.margin = margin(r = 100,  # Right margin
                                                 b = 0,  # Bottom margin
						 l = 10))
dev.off()


seurat_mat_heat <- ScaleData(object = seurat_mat_subset, assay = 'RNA', features = unique(top5_markers$gene))

pdf(file = paste('heatmap_top5_markers_',project,'.pdf',sep = ''), width = 35, height = 25)
DoHeatmap(seurat_mat_heat, assay = 'RNA', features = unique(top5_markers$gene), size = 7, lines.width = 7, raster= F) + NoLegend() + theme(axis.text.y = element_text(size = 11)) + theme(plot.margin = margin(r = 100,  # Right margin
						 b = 0,  # Bottom margin
						 l = 10))
dev.off()

q() 
