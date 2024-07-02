#Scenic_likai_Lung_cd3
#Likai 
##------ Wed May  1 06:46:45 2024 ------##

library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

setwd("/home/big/tanlikai/Lung/")


#Load Seurat Obj of lung GTD

GDTlung.rds <- 'GDTlung2023july_7p.rds'
GDTlung_s <- readRDS(GDTlung.rds)

exprMat <- GDTlung_s@assays$RNA@data
cellInfo <- GDTlung_s@meta.data

loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

loom <- build_loom("GDTlung_s_loom_2024MAY.loom", dgem=exprMat_filter)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loomPath <- system.file(package="SCENIC", "data/Demo_GDTlung_sSeurat_SCT_Preprocess_FilterLQCells.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)

exprMat <- get_dgem(loom)

loom <- open_loom('/home/big/zheng/scenic/scenic/pySCENIC_Lung/lung_output_tracks.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')



