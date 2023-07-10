
# lung project 6 patients -------------------------------------------------





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
# library(rstatix)
library(stringr)
library(magrittr)
library(Nebulosa)

# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')


file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

# future::plan(strategy = "sequential")
# ?parallelly::supportsMulticore
#single core
setwd("/home/big/tanlikai/Lung")

GDTlung.rds <- 'GDTlung_08032022.rds'

GDTlung_s <- readRDS(GDTlung.rds)




# creat 10x ---------------------------------------------------------------


raw_data_1 <- Read10X_h5('raw/surface_Falk1_gd/outs/filtered_feature_bc_matrix.h5')
raw_data_1$`Antibody Capture`

plot(x = raw_data_1$`Antibody Capture`[3,], y = raw_data_1$`Antibody Capture`[4,])

plot(x = raw_data_1$`Antibody Capture`[1,], y = raw_data_1$`Antibody Capture`[2,])

GDTlung_1 <- CreateSeuratObject(counts = raw_data_1$`Gene Expression`,project = 'p25',
                                meta.data = t(raw_data_1$`Antibody Capture`), min.features = 20)
GDTlung_1@meta.data %>% dim()
GDTlung_1@assays$RNA %>% dim()
GDTlung_1$bc_backup <- rownames(GDTlung_1@meta.data)
colnames(GDTlung_1@meta.data)

rwa_data_2_3 <- Read10X_h5('raw/surf_322_325_453_456_gd/outs/filtered_feature_bc_matrix.h5')

plot(x = rwa_data_2_3$`Antibody Capture`[4,], y = rwa_data_2_3$`Antibody Capture`[5,])
rwa_data_2_3$`Antibody Capture`

GDTlung_2_3 <- CreateSeuratObject(counts = rwa_data_2_3$`Gene Expression`,project = 'p32&p45',
                                  meta.data = t(rwa_data_2_3$`Antibody Capture`), min.features = 20)

GDTlung_2_3[['HTO']] <- CreateAssayObject(rwa_data_2_3$`Antibody Capture`[3:6,])

GDTlung_2_3$bc_backup <- rownames(GDTlung_2_3@meta.data)


rawdata_4 <- Read10X_h5('raw/surface_falk3_gd/outs/filtered_feature_bc_matrix.h5')
Read10X_h5('raw/surface_falk3_gd/outs/filtered_feature_bc_matrix.h5')$`Antibody Capture` %>% rownames()
rownames(rawdata_4$`Antibody Capture`)


rawdata_4$`Antibody Capture`[1:6, ] 

GDTlung_4_6 <- CreateSeuratObject(rawdata_4$`Gene Expression`, project = 'p71&p31&p27',
                                  meta.data = t(rawdata_4$`Antibody Capture`))

GDTlung_4_6@meta.data

GDTlung_4_6[['HTO']] <- CreateAssayObject(rawdata_4$`Antibody Capture`[1:6,])
GDTlung_4_6[['CITE']] <- CreateAssayObject(rawdata_4$`Antibody Capture`[7:23,])



GDTlung_4_6$bc_backup <- rownames(GDTlung_4_6@meta.data)


rawdata_7 <- Read10X('raw/surface_falk77_73_gd//outs/filtered_feature_bc_matrix/')

Read10X('raw/surface_falk77_73_gd//outs/filtered_feature_bc_matrix/')$`Antibody Capture` %>% rownames()

GDTlung_7_8 <- CreateSeuratObject(rawdata_7$`Gene Expression`, project = 'p73p77',
                                  meta.data = t(rawdata_7$`Antibody Capture`))

GDTlung_7_8[['HTO']] <- CreateAssayObject(rawdata_7$`Antibody Capture`[1:4,])
GDTlung_7_8[['CITE']] <- CreateAssayObject(rawdata_7$`Antibody Capture`[5:21,])



# assign tissue
# demultiplexing ----------------------------------------------------------
Feature_rast(GDTlung_1, g = 'nCount_RNA', d1 = "X1_TotalSeqC_LuPa" , d2 = "X9_TotalSeqC_LN",
             noaxis = F, axis.number = T)+xlim(0,5000)+ylim(0, 10000)+
  geom_vline(xintercept = 100)+geom_hline(yintercept = 180)+grd
GDTlung_1$X9_TotalSeqC_LN
GDTlung_1@meta.data <- GDTlung_1@meta.data %>% 
  mutate(tissue = case_when(X1_TotalSeqC_LuPa >= 100 & X9_TotalSeqC_LN < 180 ~ 'lung',
                            X1_TotalSeqC_LuPa < 100 & X9_TotalSeqC_LN >= 180 ~ 'luLN',
                            X1_TotalSeqC_LuPa >= 100 & X9_TotalSeqC_LN >= 180 ~ 'Doublet')) %>%
  mutate(tissue = factor(tissue, levels = c('luLN', 'lung', 'Doublet'))) %>% 
  `rownames<-`(GDTlung_1$bc_backup)
GDTlung_1$patient <- 'p25'
dim(GDTlung_1)
GDTlung_1 <-GDTlung_1 %>%  subset( subset = tissue %in% c('lung', 'luLN'))

table(GDTlung_1$tissue,useNA = 'always')

table(GDTlung_1$tissue,useNA = 'always')


tissuehashtag <- Feature_rast(GDTlung_1, g = 'tissue', d1 = "X1_TotalSeqC_LuPa" , d2 = "X9_TotalSeqC_LN",
                              noaxis = F, axis.number = T)+xlim(0,5000)+ylim(0, 10000)+
  geom_vline(xintercept = 100)+geom_hline(yintercept = 180)
figsave(tissuehashtag,'tissuehashtag_beforetrimming.pdf', 150,150)

GDTlung_2_3 <- NormalizeData(GDTlung_2_3, assay = "HTO", normalization.method = "CLR") %>% 
  ScaleData(assay='HTO', features = rownames(GDTlung_2_3[['HTO']]))
GDTlung_2_3 <- HTODemux(GDTlung_2_3, assay = "HTO", positive.quantile = 0.99) 

GDTlung_2_3$HTO_classification.global %>% table
GDTlung_2_3$HTO_maxID %>% table
GDTlung_2_3$HTO_classification.global

table(GDTlung_2_3$hash.ID ,useNA = 'always')

Hash_hist_2_3 <- RidgePlot(GDTlung_2_3, assay = "HTO", slot = 'scale.data', features = rownames(GDTlung_2_3[["HTO"]]), ncol = 2)
Hash_hist_2_3

figsave(Hash_hist_2_3, 'Hash_hist_2_3.pdf', 200, 150)

VlnPlot(GDTlung_2_3, 'nFeature_RNA',g = 'hash.ID')

GDTlung_2_3 %<>% RunUMAP(assay= 'HTO', features = rownames(GDTlung_2_3[['HTO']]))
DimPlot(GDTlung_2_3)
# add tissue and patien ID

GDTlung_2_3@meta.data %<>% 
  mutate(tissue = case_when(grepl('Lu', hash.ID) ~ 'lung', 
                            grepl('LN', hash.ID)~ 'luLN'),
         patient = case_when(grepl('32', hash.ID) ~ 'p32', 
                             grepl('45', hash.ID)~ 'p45'))%>%
  `rownames<-`(GDTlung_2_3$bc_backup)
GDTlung_2_3 <-GDTlung_2_3 %>%  subset( subset = tissue %in% c('lung', 'luLN'))
GDTlung_s$hash.ID %>% unique()


GDTlung_4_6 %<>% NormalizeData( assay = "HTO", normalization.method = "CLR") %>% 
  ScaleData(assay='HTO', features = rownames(GDTlung_4_6[['HTO']])) %>% 
  RunUMAP(assay= 'HTO', features = rownames(GDTlung_4_6[['HTO']])) %>% 
  RunPCA(assay= 'HTO', features = rownames(GDTlung_4_6[['HTO']])) %>% 
  HTODemux( assay = "HTO", positive.quantile = 0.99) 
DimPlot(GDTlung_4_6,reduction = 'pca')  

table(GDTlung_4_6$hash.ID, useNA = 'always')
Hash_hist_4_6<-  RidgePlot(GDTlung_4_6, assay = "HTO", features = rownames(GDTlung_4_6[["HTO"]]), ncol = 2)
figsave(Hash_hist_4_6, 'Hash_hist_4_6.pdf', 200, 250)

GDTlung_4_6@meta.data %<>% 
  mutate(tissue = case_when(grepl('P', hash.ID) ~ 'lung', 
                            grepl('L', hash.ID)~ 'luLN'),
         patient = case_when(grepl('71', hash.ID) ~ 'p71', 
                             grepl('31', hash.ID)~ 'p31',
                             grepl('27', hash.ID)~ 'p27'))%>%
  `rownames<-`(GDTlung_4_6$bc_backup)

GDTlung_4_6 %<>% subset( subset = tissue %in% c('lung', 'luLN'))

GDTlung_7_8 %<>% NormalizeData( assay = "HTO", normalization.method = "CLR") %>% 
  ScaleData(assay='HTO', features = rownames(GDTlung_7_8[['HTO']])) %>% 
  RunUMAP(assay= 'HTO', features = rownames(GDTlung_7_8[['HTO']])) %>% 
  RunPCA(assay= 'HTO', features = rownames(GDTlung_7_8[['HTO']])) 
GDTlung_7_8 %<>%  HTODemux( assay = "HTO", positive.quantile = 0.85,init = 4) 

GDTlung_7_8 %<>% RunUMAP(assay= 'HTO', features = rownames(GDTlung_7_8[['HTO']]))


DimPlot(GDTlung_7_8,reduction = 'pca',group.by = 'HTO_classification',cols = umap.colors)  
DimPlot(GDTlung_7_8,reduction = 'umap',group.by = 'HTO_classification',cols = umap.colors)  


GDTlung_7_8$HTO_classification %>% table
GDTlung_7_8@active.ident %>% table()

rownames(GDTlung_7_8[["HTO"]])
# "1-TotalSeqC-737L" ,"2-TotalSeqC-736P" ,"3-TotalSeqC-772L", "4-TotalSeqC-771P"

RidgePlot(GDTlung_7_8, assay = "HTO", features = rownames(GDTlung_7_8[["HTO"]]), ncol = 2)
VlnPlot(GDTlung_7_8,assay = 'HTO', rownames(GDTlung_7_8[['HTO']]))
GDTlung_7_8[['HTO']] <- CreateAssayObject(rawdata_7$`Antibody Capture`[1:4,])
GDTlung_7_8[['CITE']] <- CreateAssayObject(rawdata_7$`Antibody Capture`[5:21,])


FeatureScatter(GDTlung_7_8,feature1 = "1-TotalSeqC-737L",feature2 = "2-TotalSeqC-736P"  )


FeatureScatter(GDTlung_7_8,"1-TotalSeqC-737L", "2-TotalSeqC-736P" )

# manually annotation 

ViolinPlot(GDTlung_7_8, 'nFeature_RNA', group.by = 'HTO_classification')

Feature_rast(GDTlung_7_8, d1="X1_TotalSeqC_737L", d2="X2_TotalSeqC_736P", do.label = F,noaxis = F, axis.number = T)+
  xlim(0,2500)+ylim(0, 2500)+geom_hline(yintercept = 80)+geom_vline(xintercept = 150)


# library(scGate)
# 
# m1 <- gating_model(name = 'lung_73', signature = c("1-TotalSeqC-737L" ,"2-TotalSeqC-736P-" ,"3-TotalSeqC-772L-", "4-TotalSeqC-771P-"))
# m2 <- gating_model(name = 'lung_77', signature = c("1-TotalSeqC-737L-" ,"2-TotalSeqC-736P-" ,"3-TotalSeqC-772L", "4-TotalSeqC-771P-"))
# m3 <- gating_model(name = 'luLN_73', signature = c("1-TotalSeqC-737L-" ,"2-TotalSeqC-736P" ,"3-TotalSeqC-772L-", "4-TotalSeqC-771P-"))
# m4 <- gating_model(name = 'luLN_77', signature = c("1-TotalSeqC-737L-" ,"2-TotalSeqC-736P-" ,"3-TotalSeqC-772L-", "4-TotalSeqC-771P+"))
# 
# GDTlung_7_8$lung_73_UCell
# 
# GDTlung_7_8 <- scGate(GDTlung_7_8, m1, assay = 'HTO' )
# GDTlung_7_8 <- scGate(GDTlung_7_8, m2, assay = 'HTO')
# GDTlung_7_8 <- scGate(GDTlung_7_8, m3, assay = 'HTO')
# GDTlung_7_8 <- scGate(GDTlung_7_8, m4, assay = 'HTO')
# 
# GDTlung_7_8$luLN_77_UCell
# 
# 
# DimPlot(GDTlung_7_8, group.by = 'lung_77_UCell', reduction = 'pca')


GDTlung_7_8@meta.data[ ,c("X1_TotalSeqC_737L","X2_TotalSeqC_736P", "X3_TotalSeqC_772L","X4_TotalSeqC_771P")] <-  t(GDTlung_7_8[["HTO"]]@scale.data)

GDTlung_7_8$bc_backup <- GDTlung_7_8@meta.data %>% rownames()


ViolinPlot(GDTlung_7_8,'lung_77_UCell' , group.by = 'HTO_classification')
Feature_rast(GDTlung_7_8, d1="X1_TotalSeqC_737L", d2="X2_TotalSeqC_736P", noaxis = F, axis.number = T)+
  geom_vline(xintercept = 1)+geom_hline(yintercept = 2)

GDTlung_7_8@meta.data$manul
GDTlung_7_8@meta.data  %<>%  mutate(manul= case_when(X1_TotalSeqC_737L > 1 & X2_TotalSeqC_736P >2 ~'DP' ,
                                                     X1_TotalSeqC_737L > 1 ~ 'luLN_73',
                                                     X2_TotalSeqC_736P >2 ~ 'lung_73'  )) %>%  `rownames<-`(GDTlung_7_8$bc_backup)
Feature_rast(GDTlung_7_8, 'manul')

Feature_rast(GDTlung_7_8, d1="X1_TotalSeqC_737L", d2="X3_TotalSeqC_772L", noaxis = F, axis.number = T)+
  geom_vline(xintercept = 1)+geom_hline(yintercept = 0.5)
GDTlung_7_8@meta.data  %<>% 
  mutate(
    manul = replace(manul,  X3_TotalSeqC_772L > 0.5, 'luLN_77'),
    manul = replace(manul, X1_TotalSeqC_737L >1 & X3_TotalSeqC_772L > 0.5, 'DP')
  )

DimPlot(GDTlung_7_8, group.by = 'manul', reduction = 'pca')


Feature_rast(GDTlung_7_8, d1="X1_TotalSeqC_737L", d2="X4_TotalSeqC_771P", noaxis = F, axis.number = T)+
  geom_vline(xintercept = 1)+geom_hline(yintercept = -0.8)
GDTlung_7_8@meta.data  %<>% 
  mutate(
    manul = replace(manul,  X4_TotalSeqC_771P > -0.8, 'lung_77'),
    manul = replace(manul, X1_TotalSeqC_737L >1 & X4_TotalSeqC_771P > -0.8, 'DP')
  )




Feature_rast(GDTlung_7_8, d1="X3_TotalSeqC_772L", d2="X2_TotalSeqC_736P", noaxis = F, axis.number = T)+
  geom_vline(xintercept = 0.5)+geom_hline(yintercept = 2)

GDTlung_7_8@meta.data  %<>% 
  mutate(
    manul = replace(manul, X3_TotalSeqC_772L >0.5 & X2_TotalSeqC_736P > 2, 'DP')
  )

Feature_rast(GDTlung_7_8, d1="X4_TotalSeqC_771P", d2="X2_TotalSeqC_736P", noaxis = F, axis.number = T)+
  geom_vline(xintercept = -0.8)+geom_hline(yintercept = 2)

GDTlung_7_8@meta.data  %<>% 
  mutate(
    manul = replace(manul, X4_TotalSeqC_771P >-0.8 & X2_TotalSeqC_736P > 2, 'DP')
  )


Feature_rast(GDTlung_7_8, d1="X3_TotalSeqC_772L", d2="X4_TotalSeqC_771P", noaxis = F, axis.number = T)+
  geom_vline(xintercept = 0.5)+geom_hline(yintercept = -0.8)
GDTlung_7_8@meta.data  %<>% 
  mutate(
    manul = replace(manul, X3_TotalSeqC_772L >0.5 & X4_TotalSeqC_771P > -0.8, 'DP')
  )

GDTlung_7_8$manul %>% table

GDTlung_7_8@meta.data %<>% 
  mutate(tissue = case_when(grepl('lung', manul) ~ 'Pulmonary', 
                            grepl('luLN', manul)~ 'LN'),
         patient = case_when(grepl('77', manul) ~ 'p77', 
                             grepl('73', manul)~ 'p73'))%>%
  `rownames<-`(GDTlung_7_8$bc_backup)
GDTlung_7_8 %<>% subset( subset = tissue %in% c('Pulm', 'LN'))



DimPlot(GDTlung_7_8, group.by = 'manul', reduction = 'pca')

DefaultAssay(GDTlung_7_8) <- 'RNA'


GDTlung <- list(GDTlung_1, GDTlung_2_3, GDTlung_4_6,GDTlung_7_8) %>% setNames(c('p1', 'p23', 'p456', 'p78'))