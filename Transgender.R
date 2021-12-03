
# Transgender  ------------------------------------------------------------


library(Seurat)
library(ggplot2)
library(purrr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(cowplot)
library(ggrastr)
library(ggrepel)
library(xlsx)
library(stringr)
library(magrittr)
library(openxlsx)
library(harmony)
source('/home/big/tanlikai/script/rscripts/funcs.r')
file.edit('/home/big/tanlikai/script/rscripts/funcs.r')


# set working dr ----------------------------------------------------------

setwd('/home/big/tanlikai/Transgender/')

TG_Casar <- readRDS('all_samples_cca_ds_subclustering_tcr.rds')

TG_Casar@meta.data

Feature_rast(TG_Casar, 'cloneType')
 


# creat project -----------------------------------------------------------


dir_all <- list.files('cellranger_out/')
# TCR library
dir_TCR <- str_subset(dir_all, 'tcr')
# gdTCR library
dir_gdTCR<- str_subset(dir_all, 'gdTCR')
# gene expression library
dir_GEX <- str_subset(dir_all, 'tcr|gdTCR', negate = T)

samplesheet <- data.frame(sample = dir_GEX, donor  = str_extract(dir_GEX, 'sample0\\d\\d'), tp = str_extract(dir_GEX, '\\w\\w$'))




# gene expression ---------------------------------------------------------

# read raw gene expression data  
RAW_GEX <- map(dir_GEX, ~ Read10X_h5(paste0('cellranger_out/', .x,'/outs/filtered_feature_bc_matrix.h5')) ) %>%  setNames(dir_GEX)


names(RAW_GEX)

names(RAW_GEX$sample001_BL)
# [1] "Gene Expression"  "Antibody Capture"

RAW_GEX$sample001_BL$`Gene Expression` %>%  str()

map(RAW_GEX, ~ dim(.x$`Gene Expression`))


RAW_GEX$sample001_BL$`Antibody Capture` %>%  str()


rownames(RAW_GEX$sample001_BL$`Antibody Capture`)



# generate Seurat  project ------------------------------------------------



TG_pos_inte <- map2(RAW_GEX, dir_GEX, ~ CreateSeuratObject(counts = .x$`Gene Expression`, project = .y, min.features =30) %>% 
                      AddMetaData(col.name = c('donor', 'tp'), metadata = c(str_extract(dir_GEX, 'sample0\\d\\d'), str_extract(dir_GEX, '\\w\\w$')) ))


str(TG_pos_inte$sample001_BL)
summary(TG_pos_inte$sample001_BL)



for (x in dir_GEX) {
  TG_pos_inte[[x]][['CITE']] <- CreateAssayObject(RAW_GEX[[x]]$`Antibody Capture`)
  
}


# QC ----------------------------------------------------------------------

# calculate mitochonrial, ribosomal, and hispa  gene content

TG_pos_inte <- map(TG_pos_inte, ~ PercentageFeatureSet(.x, '^MT', col.name =  'percent.mito') %>% 
                   PercentageFeatureSet('^RP', col.name = 'percent.ribo')    %>% 
                   PercentageFeatureSet('^HSPA', col.name =  'percent.hspa') )

TG_pos_inte$sample001_BL@meta.data


dir.create('figs')

QCvio <- map(TG_pos_inte, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito','percent.ribo'),
                                    colors ='blue'  ,box = T, jitter = T , ncol = 4  )) %>% 
  PG(labels = dir_GEX, ncol = 2) %T>%
  figsave('beforeQC_violin.pdf', 220, 270)

QCvio


QC_scatter <- map(TG_pos_inte, ~ Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', 
                                            noaxis = F, axis.number = T)+grd+
                    geom_smooth(method = "lm")+
                    
                    scale_x_continuous(breaks = seq(0, 15000, 1000), limits = c(0,15000))+
                    scale_y_continuous(breaks = seq(0, 3000, 500), limits = c(0,3000))+
                    geom_hline(yintercept = c(400,2000))+geom_vline(xintercept = c(1200,7000)))  %>% 
  PG(labels = dir_GEX)%T>%
  figsave('beforeQC_scatter.pdf',600,400)
QC_scatter



