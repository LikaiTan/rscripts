
# Pulm CD4 and CD8 project
library(Seurat)
library(ggplot2)
library(Nebulosa)
library(ggrastr)
library(cowplot)
library(purrr)
library(ggrepel)
library(dplyr)
library(magrittr)
library(tibble)
library(stringr)
library(Cairo)
library(Matrix)
library(ggalluvial)
library(gridExtra)
library(openxlsx)
library(rstatix)
# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')


file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

multicores()
 

#paralell computation

# work dir ----------------------------------------------------------thx

setwd('/home/big/tanlikai/Lung/')
CD4CD8RDS <- 'CD4CD8_integrated_2022_8p_withTCR_2023AUG.rds'
figpath_ni <-  "/home/big/googledrive/Lungfigures/"




CD4CD8 <- readRDS(CD4CD8RDS)
Feature_rast(CD4CD8)
# read_rawdata ---------------------------------------------------------------

# directories to raw gene expression and CITEseq data

dirs <- list.dirs(path = 'abt/raw') %>% str_subset('surface.+outs$' ) %>% grep(314, . ,value = T,invert = T)



proj <- c('p32', 'p45','p73', 'p77', 'p25_CD4', 'p25_CD8','p27','p71')
CD4 <- c('p32', 'p45','p73', 'p77','p25', 'p25','p27','p71')

# read raw GEX and CITEseq data 
rawdata <- dirs %>% map(~ Read10X_h5(paste0(.x, '/filtered_feature_bc_matrix.h5')) ) %>% setNames(proj)
# creat seurat project of gene expression datasepeately into a big list

PulmTcell <- map2(rawdata,proj, ~
                   CreateSeuratObject(counts = .x$`Gene Expression`,project = .y,
                                       min.features = 20) %>% 
                    AddMetaData(str_extract(.y, 'p\\d\\d'),col.name = 'patient')
                 )

PulmTcell$p32


map(PulmTcell, ~ .@meta.data %>% colnames )
map(PulmTcell, ~ .@meta.data %>% ncol )

map(PulmTcell, ~ dim(.))


# demultiplexing ----------------------------------------------------------


rawdata
#  citeseq and hashtaq data

hash <- unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )) %>%  str_subset('^\\d_')
hash

cite <- setdiff(unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )), hash)
cite

# integrate HTO and CITESEQdata into seurat projects 

for (x in proj) {
  PulmTcell[[x]][['HTO']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(hash, rownames(rawdata[[x]]$`Antibody Capture`)),])
  PulmTcell[[x]][['CITE']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(cite, rownames(rawdata[[x]]$`Antibody Capture`)),])

  
}



  proj


# change names of HTO to more recognizable names 

  PulmTcell[["p32"]][["HTO"]]@counts
  PulmTcell[["p32"]][["HTO"]]@data
  
  PulmTcell[["p32"]][['HTO']]["counts"]
  
  
  Feature_rast(PulmTcell$p32, d1 = "1-TotalSeqC-1502C", d2 ="9-TotalSeqC-1502E")
  
  
  rownames(PulmTcell[["p27"]][['HTO']]["counts"]) %>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+') %>% paste0('S',.)
  
for (x in proj) {
  rownames(PulmTcell[[x]][['HTO']]["counts"]) <- rownames(PulmTcell[[x]][['HTO']]["counts"])  %>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+') %>% paste0('S',.)
  rownames(PulmTcell[[x]][['HTO']]["data"])  <- rownames(PulmTcell[[x]][['HTO']]["data"])  %>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+')%>% paste0('S',.)
  
  
}

hash_sample <- PulmTcell %>% map(~ .[['HTO']] %>% rownames)
hash_sample
PulmTcell$p32@meta.data






Feature_rast(CD4CD8 %>% subset(patient == "p32"), "tissue", d1 = "S1502C", d2 = "S1502E", noaxis = F, colorset = c("blue", "red"), mythe = F)

PulmTcell$p32


# normalization
for (i in proj) {
  PulmTcell[[i]] %<>%  NormalizeData( assay = "HTO", normalization.method = "CLR") %>% 
    ScaleData(assay='HTO', features = rownames( .[['HTO']]), vars.to.regress= 'nCount_HTO') %>% 
    HTODemux( assay = "HTO", positive.quantile = 0.9, nstarts = 20)
}


for (i in proj) {
  PulmTcell[[i]] %<>%  NormalizeData( assay = "CITE", normalization.method = "CLR") %>% 
    ScaleData(assay='CITE', features = rownames( .[['CITE']]), vars.to.regress= 'nCount_CITE') 
    
}


map(PulmTcell, ~ .@meta.data %>% colnames )
map(PulmTcell, ~ .@meta.data %>% ncol )

for (x in proj) {
  DefaultAssay(PulmTcell[[x]] ) <- 'HTO'
}

saveRDS(PulmTcell,'abt/lungTcell_pre_inte.rds')

PulmTcell <- readRDS('abt/lungTcell_pre_inte.rds')


CD4CD8 %>% subset(patient == "p32")
Feature_rast(CD4CD8 %>% subset(patient == "p25"), "tissue", d1 = "S1502C", d2 = "S1502E")
Feature_rast(CD4CD8, "tissue", d1 = "S1502C", d2 = "S1502E")


PulmTcell$p32$hash.ID
map2(PulmTcell,hash_sample,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]], assay = 'HTO',
                                             g = 'hash.ID',noaxis =F, axis.number=T)+
       # xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.65)+
       geom_vline(xintercept = 0.7)) %>% PG(ncol =4) %T>%
  figsave('hash.id.automatic.pdf', 400, 100)

# the automiatic demultiplexing by seurat is not very accurate 

# add HTO to metadata

for (x in proj) {
 PulmTcell[[x]]  %<>% AddMetaData(t(PulmTcell[[x]]@assays$HTO@data), col.name = rownames(PulmTcell[[x]]@assays$HTO@data))
}

for (x in proj) {
  PulmTcell[[x]]  %<>% AddMetaData( rownames(PulmTcell[[x]]@meta.data), 
                                    col.name = 'bc_backup')
}

PulmTcell$p73@meta.data<- PulmTcell$p73@meta.data[,c(1:14, 16,15, 17)] 

PulmTcell$p77@meta.data<- PulmTcell$p77@meta.data[,c(1:14, 16,15, 17)] 




PulmTcell$p32@assays$HTO@data %>% rownames()

#  demultiplexing manually

aa<- map(PulmTcell,~ rownames(.@assays$HTO@data)
    )


bb<- map(PulmTcell,~     colnames(.@meta.data)[15:16]

)

map2(aa, bb, ~ .x == .y)

syms(hash_sample)
newmeta <- list()

# marker tissues based on HTO
for (x in proj) {
  DefaultAssay(PulmTcell[[x]] ) <- 'RNA'
}

map2(PulmTcell,bb,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]], assay = 'RNA',
                                    g = 'hash.ID',noaxis =F, axis.number=T)+
       # xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.65)+
       geom_vline(xintercept = 1)) %>% PG(ncol =4)%T>%
  figsave('hash.id.automatic.pdf', 400, 100)




newmeta <- map(PulmTcell,  ~
                .x@meta.data  %<>%     
                  mutate(tissue = case_when(.[[15]] >= 1 & .[[16]] < 0.55 ~ 'Pulm',
                                            .[[15]]  < 1 & .[[16]] >= 0.65 ~ 'LN',
                                            .[[15]]  >= 1 & .[[16]] >= 0.65 ~ 'DP'))%>%
  `rownames<-`(.x$bc_backup)
)


newmeta$p32  %<>%    
  mutate(tissue = case_when(.[[15]] >= 1 & .[[16]] < 0.55 ~ 'LN',
                            .[[15]]  < 1 & .[[16]] >= 0.65 ~ 'Pulm',
                            .[[15]]  >= 1 & .[[16]] >= 0.65 ~ 'DP'))%>%
  `rownames<-`(.$bc_backup)


map2(newmeta,bb,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]], 
                                    g = 'tissue',noaxis =F, axis.number=T)+
       # xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.7)+
       geom_vline(xintercept = 1)) %>% PG(ncol =4)


# newmeta$p77  %<>%    
#   mutate(tissue = case_when(.[[15]] >= 1 & .[[16]] < 0.7 ~ 'LN',
#                             .[[15]]  < 1 & .[[16]] >= 0.7 ~ 'Pulm',
#                             .[[15]]  >= 1 & .[[16]] >= 0.7 ~ 'DP'))%>%
#   `rownames<-`(.$bc_backup)
# 
# 
# 
# newmeta$p73  %<>%    
#   mutate(tissue = case_when(.[[15]] >= 0.6 & .[[16]] < 0.8 ~ 'LN',
#                             .[[15]]  < 0.6 & .[[16]] >= 1 ~ 'Pulm',
#                             .[[15]]  >= 0.6 & .[[16]] >= 1 ~ 'DP'))%>%
#   `rownames<-`(.$bc_backup)
# 
# newmeta$p25_CD8 %<>%    
#   mutate(tissue = case_when(.[[15]] >= 1 & .[[16]] < 0.5 ~ 'Pulm',
#                             .[[15]]  < 1 & .[[16]] >= 1 ~ 'LN',
#                             .[[15]]  >= 1 & .[[16]] >= 0.7 ~ 'DP'))%>%
#   `rownames<-`(.$bc_backup)
# 
# 
# newmeta$p71 %<>%    
#   mutate(tissue = case_when(.[[15]] >= 1 & .[[16]] < 0.3 ~ 'Pulm',
#                             .[[15]]  < 1.2 & .[[16]] >= 0.5 ~ 'LN',
#                             .[[15]]  >= 1 & .[[16]] >= 0.5 ~ 'DP'))%>%
#   `rownames<-`(.$bc_backup)
# 
# newmeta$p32  %<>%
#   mutate(tissue = case_when(.[[15]] >= 0.8 & .[[16]] < 0.5 ~ 'LN',
#                             .[[15]]  < 0.8 & .[[16]] >= 0.5 ~ 'Pulm',
#                             .[[15]]  >= 0.8 & .[[16]] >= 0.5 ~ 'DP'))%>%
#   `rownames<-`(.$bc_backup)
# 
# map2(newmeta,hash_sample,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]], assay = 'HTO',
#                                              g = 'tissue',noaxis =F, axis.number=T)+
#        # xlim(0,5)+ylim(0,5)+
#        geom_hline(yintercept = 0.7)+
#        geom_vline(xintercept = 1)) %>% PG(ncol =4)
#                 
  
for (i in proj) {
  PulmTcell[[i]]@meta.data <- newmeta[[i]]
}

map2(newmeta,bb,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]],
                                             g = 'tissue',noaxis =F, axis.number=T)+
       ggtitle(unique(.x$orig.ident))) %>% PG(ncol =4) %T>%
  figsave('hash.id.mannual.pdf', 400, 200)


# gene expression ---------------------------------------------------------


for (x in proj) {
  DefaultAssay(PulmTcell[[x]] ) <- 'RNA'
}




# QC ----------------------------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# calculate percent of mitochondrial  ribosam  and hspa genes

PulmTcell <- map(PulmTcell, ~ PercentageFeatureSet(.x, '^MT', col.name =  'percent.mito') %>% 
             PercentageFeatureSet('^RP', col.name = 'percent.ribo') )  
PulmTcell$p32$tissue
  
# QC plot

QCvio <- map(PulmTcell, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                colors =umap.colors  ,box = T, jitter = F , ncol = 3  )) %>% 
    PG(labels = proj, ncol = 2) %T>% 
  figsave('beforeQC_violin.pdf',300, 250)

  QCvio
  
  
QC_scatter <- map(PulmTcell, ~ Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', 
                                      noaxis = F, axis.number = T)+grd+
                    geom_smooth(method = "lm")+
                    
                    scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
                    scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
                    geom_hline(yintercept = c(600,2000))+geom_vline(xintercept = c(1000,6000)))  %>% 
  PG(labels = proj, ncol = 2)%T>%
  figsave('beforeQC_scatter.pdf',600,600)
QC_scatter


# data cleaning 
for (x in proj) {
  PulmTcell[[x]] <- subset(PulmTcell[[x]], subset =  nCount_RNA %in% 1000:6000 &
                       nFeature_RNA %in% 800:2000 &  percent.mito <15 &
                       tissue %in% c('Pulm','LN'))
}
QCvio_clean <- map(PulmTcell, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                              colors =umap.colors  ,box = T, jitter = T , ncol = 3  )) %>% 
  PG(labels = proj, ncol = 2) %T>%
  figsave('afterQC_violin_CD4CD8.pdf',200,270)
QCvio_clean

map(PulmTcell,~ dim(.x))



# $CD4
# [1] 33538  2398
# 
# $CD8
# [1] 33538  3626


# normalization and scaling -----------------------------------------------
# scale and find high var genes 
multicores()

PulmTcell %<>% map(~   NormalizeData(.x, normalization.method = 'LogNormalize',
                               scale.factor = 10000, assay = 'RNA')%>% 
                     CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F) %>%
                     FindVariableFeatures(assay = 'RNA',nfeatures = 3500, selection.method = 'vst') )

PulmTcell$p32@assays$RNA@scale.data


for (x in proj) {
  PulmTcell[[x]]@assays$RNA@var.features <- PulmTcell[[x]]@assays$RNA@var.features%>%  
    str_subset('^RP|^MT|^HIST|^TRA|^TRB|^HSP', negate = T)
}

saveRDS(PulmTcell, 'lungTcell_pre_inte.rds')

map(PulmTcell, ~ .x@assays$RNA@var.features %>% length() )

# integration -------------------------------------------------------------
anchors <- FindIntegrationAnchors(PulmTcell, dims = 1:60)

CD4CD8 <- IntegrateData(anchors, dims = 1:60)

# dimensional reduction by PCA and UMAP -----------------------------------


DefaultAssay(CD4CD8) <- "integrated"
# scale data and run PCA
CD4CD8 %<>% ScaleData( vars.to.regress = c('patient', 
                                           "percent.mito",
                                                   "S.Score",
                                                   'G2M.Score',
                                                   # "percent.ribo",
                                                   # 'nCount_RNA',
                                                   'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)




ElbowPlot(CD4CD8, ndims = 50)

# jackstraw is a way to choose significant PCs
# CD4CD8 %<>% JackStraw( num.replicate = 100, dims = 80)%>%
#   ScoreJackStraw(dims = 1:80) 
# CD4CD8 %>% JackStrawPlot( dims = 1:50 ) %T>%
#   figsave('CD4CD8_Pulm_5p.jackstraw.pdf' , w = 400, h = 400)

CD4CD8 <- RunUMAP(object = CD4CD8, dims = 1:43, 
                   reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:43))

  CD4CD8 <- FindClusters(CD4CD8, resolution = 0.8)

Feature_rast(CD4CD8,'tissue',facets = 'patient')
saveRDS(CD4CD8, CD4CD8RDS)
# CD4CD8$integrated_snn_res.0.6
Feature_rast(CD4CD8, c('CCR6', 'RORC', 'FOXP3', 'DPP4', 'IL23R', 'IL17A', 'IL17F', 'AREG', 'IFNG', 'IL22'), assay = 'RNA')


# CITEseq -----------------------------------------------------------------

CD4CD8 %<>% 
  ScaleData(assay='CITE', features = rownames( .[['CITE']]), vars.to.regress= c('nCount_CITE', 'patient')) 


CD4CD8@assays$CITE@data




CD4CD8 %<>% 
  ScaleData(assay='CITE', vars.to.regress= c('nCount_CITE', 'patient')) 

# RNA vs CITE
surfacemakers <- c('CD4','CD8A','CD8B', 'CXCR3','CCR6',  'ITGAE',   'CD69','CCR7','KLRB1','CD27','KLRG1', 'IL7R',
                   'DPP4', 'ITGA1', 'KLRD1', 'PTPRC', 'PDCD1')


surface_RNA <- Feature_rast(subset(CD4CD8, patient %in% c('p27', 'p71')), surfacemakers, sz = 0.3)

surface_CITE <- Feature_rast(subset(CD4CD8, patient %in% c('p27', 'p71')),sz = 0.3, CD4CD8@assays$CITE@data %>% rownames(), assay = 'CITE', color_grd = 'grd')

RNAvsCITE  <- PG(list(surface_RNA, surface_CITE), labels = c('RNA', 'CITE'), ncol = 1,vjust = -10) %T>%  figsave('RNAvsCITE.pdf',200, 290)

PG(list(surface_RNA, surface_CITE), labels = c('RNA', 'CITE'), ncol = 1,label_y = -2)

map(PulmTcell, ~ .x@assays$CITE@counts %>% rownames )

ViolinPlot(CD4CD8, c('CD49a.protein', 'CD103.protein', 'CD4.protein', 'CD45RA.protein'), assay = 'CITE', group.by = 'patient'  ,colors = umap.colors, box = T)


CITEpatients <- c('p27', 'p71', 'p77', 'p73')
CD4CD8_cite <- subset(CD4CD8, patient %in% CITEpatients) 
DefaultAssay(CD4CD8_cite) <- 'CITE'

CD4CD8_cite

Feature_rast(CD4CD8, "PTPRC")

Feature_rast(CD4CD8_cite, "CD103.protein", assay = "CITE")

Feature_rast(CD4CD8_cite, CD4CD8@assays$CITE %>%  rownames(), sz = 0.1, ncol = 5,assay = 'CITE', colorgrd = "grd2")  %T>%
  figsave('CD4CD8_Citeseq_allmarkers_rast.pdf', 210, 170)


Feature_rast(CD4CD8_cite,  c('CD103.protein', 'CD49a.protein'), sz = 0.1)


Feature_rast(GDTlung_cite, "IL7R.protein")

# mark CD4 and CD8 based on CITEseq ---------------------------------------

CD4CD8@assays$CITE@data
# change name of cite seq antibodies

CD4CD8 %<>% ScaleData( assay = 'CITE',vars.to.regress = c('patient', 'nCount_CITE'))
rownames(CD4CD8@assays$CITE@data) %<>% str_replace('-TotalSeqC', '.protein')
rownames(CD4CD8@assays$CITE@counts) %<>% str_replace('-TotalSeqC', '.protein')
rownames(CD4CD8@assays$CITE@scale.data) %<>% str_replace('-TotalSeqC', '.protein')



Feature_rast(CD4CD8, d1='CD4.protein', d2='CD8.protein', assay = 'CITE', 
             color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8, d1='CD4', d2='CD8B',  color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8,  g='patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T) +facet_grid(~ patient)+geom_hline(yintercept = 1)+geom_vline(xintercept = 1)

dim(CD4CD8)


Feature_rast(CD4CD8, c('CD4', 'CD8A', 'CD8B'))


CD4CD8 %<>% AddMetaData(FetchData(CD4CD8, c('CD4.protein', 'CD8.protein'), slot = 'data'))
# CD4 CD8 demultiplexing 


Feature_rast(CD4CD8, c('ident', 'CD4CD8', 'tissue'), ncol = 2 ) %T>% 
  figsave('CD4CD8_UMAP_rawcluster.pdf', 400, 400) 



CD4CD8@meta.data %<>% mutate(CD4CD8= case_when(orig.ident == 'p25_CD4'~ 'CD4', orig.ident == 'p25_CD8'~ 'CD8', 
     cite_CD4.protein >= 1.2 & cite_CD8.protein >= 1 ~ 'DP',
     cite_CD8.protein >= 1 ~ 'CD8',
     cite_CD4.protein >= 1.2 ~'CD4',
     cite_CD4.protein < 1.2 & cite_CD8.protein < 1 ~ 'DN'))

(Feature_rast(CD4CD8 %>% subset(patient != 'p25'),  g='CD4CD8',  d1='CD4.protein', d2='CD8.protein', assay = 'CITE',colorset =  alpha(c('purple','green' ), 0.2), noaxis = F,  axis.number = T) +
    # facet_wrap(~ patient,ncol = 2)+
    ggtitle('CITEseq on CD4 and CD8')+geom_hline(yintercept = 1)+geom_vline(xintercept = 1) )%T>%  figsave('CD4andCD8staining_cite_new_allinone.pdf', 200, 200)

table(CD4CD8$CD4CD8)

Feature_rast(CD4CD8, 'CD4CD8','ident',  colorset = 'gg')

Feature_rast(CD4CD8, c('CD4.protein', 'CD8.protein'), facets = 'CD4CD8', assay = 'CITE')

Feature_rast(CD4CD8, 'ident',facets = 'CD4CD8',  colorset = 'gg')

table(CD4CD8$CD4CD8)

ViolinPlot(CD4CD8, 'nCount_RNA', group.by = 'CD4CD8')



saveRDS(CD4CD8, CD4CD8RDS)

# DP cells are significantly larger than other cells, thus we remove it as doublets
# CD4CD8  %<>% subset(subset = CD4CD8 %in% c('CD4', 'CD8', 'DN'))


Feature_rast(CD4CD8,c('CD103.protein', 'CD49a.protein', 'CD69.protein'), assay = 'CITE', color_grd = 'grd')




CD4CD8$Cell_cluster<- Idents(CD4CD8)

Feature_rast(CD4CD8, c('CD4.protein', 'CD8.protein'), assay = 'CITE', color_grd = 'grd')
Feature_rast(CD4CD8, 'CD4CD8', sz = 0.5, do.label = F,
        colorset =      alpha(c('purple','green' ), 0.8)
             )


Feature_rast(CD4CD8)


Feature_density(CD4CD8,c('AREG', 'CTLA4', 'CSF1', 'IL2RA', 'TNFRSF18', 'GATA3', 'LGALS3'))


# CD45RARO ----------------------------------------------------------------
CD4CD8_cite@assays$CITE@scale.data
(Feature_rast(CD4CD8_cite,
              g = 'CD45RO.protein',
              #w facets = c('CD4CD8', 'tissue'),
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              # facets = 'CD4CD8', 
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.5, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
    ) %T>%  figsave('CD4CD8_CD45RA_CD27.pdf', 60, 60)


(Feature_rast(CD4CD8_cite,
              g = "T_pheno",
            facets = c('tissue'),  
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
                
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.5, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
)

bc_naive <- colnames(   
  subset(CD4CD8_cite,CD27.protein > 0.5 &  CD45RA.protein > 1.2 )
  )


# test <-  subset(CD4CD8_cite,CD27.protein > 0.5 &  CD45RA.protein > 1.2 )

bc_cm <- colnames(   
  subset(CD4CD8_cite, CD27.protein > 0.5 &  CD45RA.protein <= 1.2 )
)

bc_em <- colnames(   
  subset(CD4CD8_cite,CD27.protein <= 0.5 &  CD45RA.protein <= 1.2 )
)

bc_tmra <- colnames(   
  subset(CD4CD8_cite,CD27.protein <= 0.5 &  CD45RA.protein > 1.2 )
)



CD4CD8_cite@meta.data  %<>%  mutate(T_pheno= case_when( 
  bc_backup %in% bc_naive ~ 'naive',
  bc_backup %in% bc_cm ~ 'Tcm',
  bc_backup %in% bc_em ~ 'Tem',
  bc_backup %in% bc_tmra ~ 'Temra'
  
  )   )

Feature_rast(CD4CD8_cite, 'T_pheno', facets = c('tissue'), 
             colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev()
) 


(Feature_rast(CD4CD8_cite,
              g = "T_pheno",
              facets = c('tissue'),  
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
              
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.5, linewidth = 0.2)+
    geom_vline(xintercept = 1.2,linewidth = 0.2)
)



Feature_rast(CD4CD8_cite, 'T_pheno', facets = c('tissue', 'CD4CD8'), 
             colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev()
             )  %T>%
  figsave('CD4CD8_TMRA_UMAP.pdf', 200, 200)


Feature_rast(CD4CD8_cite, 'T_pheno', facets = c('tissue'), colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev()
             ) 


(Feature_rast(CD4CD8_cite,
              g = 'ID',
              # sz=0.2,
              mythe =F,
              # facets = c( 'tissue'),
              d1 ='CD49a.protein', d2 =  'CD103.protein',
              colorset = c(brewer.pal(9,'Blues')[6:9], brewer.pal(9,'Reds')[6:9]),
              # colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 6) %>% rev(),
              do.label = F,
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 1.2, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
    # theme(axis.text = element_text(size = 10))
)



(Feature_rast(CD4CD8_cite,
              g = 'T_pheno',
              #w facets = c('CD4CD8', 'tissue'),
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
              
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.5, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
)

# clean and 2nd clustering ------------------------------------------------
# remove cells that not T cells 
ViolinPlot(CD4CD8, 'nFeature_RNA')

Feature_rast(CD4CD8, noaxis = F, axis.number = T)
CD4CD8  %<>% subset(UMAP_1 >= -5 & CD4CD8 %in% c('CD4', 'CD8') )

CD4CD8  %<>% subset(nFeature_RNA >= 800 )


DefaultAssay(CD4CD8) <- "RNA"

# CD4CD8<- PercentageFeatureSet(CD4CD8, '^HSPA', col.name =  'percent.hspa') 
DefaultAssay(CD4CD8) <- "integrated"
CD4CD8 %<>% ScaleData( vars.to.regress = c(
                                           "percent.mito",
                                           "S.Score",
                                           'G2M.Score',
                                           "percent.ribo",
                                           'patient',
                                           'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)

ElbowPlot(CD4CD8, ndim = 100)

CD4CD8 %<>% JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
CD4CD8 %>% JackStrawPlot( dims = 1:50 ) %T>%
  figsave('CD4CD8_Pulm_8p.jackstraw.pdf' , w = 400, h = 400)
saveRDS(CD4CD8, CD4CD8RDS)

# set.seed(123)
CD4CD8 <- RunUMAP(object = CD4CD8, dims = 1:40, seed.use = 1007,
                  reduction = 'pca', min.dist = 0.1) %>%
  FindNeighbors(dims = c(1:40))
for (i in seq(0.7,1.8,0.1) %>% rev()) {
  CD4CD8 <- FindClusters(CD4CD8, resolution = i, random.seed = 123)
}
Feature_rast(CD4CD8, paste0('integrated_snn_res.',seq(0.7,1.5,0.1)), sz = 0.2, ncol = 4)  %T>%
  figsave('CD4CD8_clustering_allresolution.pdf', 300, 300)

Feature_rast(CD4CD8, 'CD4CD8')

CD4CD8@reductions$umap@cell.embeddings[, 'UMAP_1'] <-  -CD4CD8@reductions$umap@cell.embeddings[, 'UMAP_1']

table(CD4CD8@active.ident)




CD4CD8$integrated_snn_res.1.8

ClusterCompare(CD4CD8, '1', '11', group.by = 'integrated_snn_res.1.2', assay = 'RNA')



CD4CD8$Cell_cluster<- factor(paste0('c', (as.numeric(CD4CD8$integrated_snn_res.1.2 )+1)),
                         levels = paste0('c', 1:18)
)
# CD4CD8$Cell_cluster <-  Idents(CD4CD8)
Feature_rast(CD4CD8, 'Cell_cluster')
Feature_rast(CD4CD8, "integrated_snn_res.1.2")
table(Idents(CD4CD8))




ClusterCompare(CD4CD8, 'c6', 'c17')

# Cluster adjustment ------------------------------------------------------

Feature_rast(CD4CD8, facets = 'CD4CD8')

Idents(CD4CD8)


Feature_rast(CD4CD8, c('TOX', 'HAVCR2', 'PDCD1'))

# ClusterCompare(CD4CD8, '4', '12')

ViolinPlot(CD4CD8, 'nCount_RNA')

Feature_rast(CD4CD8, c("integrated_snn_res.1.2", "Cell_cluster", "Cell_pheno"))


Cell_pheno <-  c(
  'TMRA_1_LG',
  'TCM_2_M',
  'TRM_3_LG',
  'Th17/Tc17_M',
  'Tfh_LN',
  'Naive_2_LN',
  'TEMRA_2_LG',
  'TRM_2_LG',
  'unidentified_LG',
  'TCM_3_M',
  'TRM_1_LG',
  'TCM_1_M',
  'TCM_3_M',
  'Naive_2_LG',
  'TEM_M',
  'TEMRA_1_LG',
  'Naive_1_LN',
  'Treg_LN'
  
  
  
)



phenotable <- data.frame(Cell_cluster = paste0('c', 1:18) , Cell_pheno) %>% arrange(Cell_pheno)

phenotable$Cluster_pheno <- paste0(phenotable$Cell_cluster,': ', phenotable$Cell_pheno) 
phenotable %<>% mutate(Cluster_pheno = factor(Cluster_pheno, levels = unique(Cluster_pheno)))

phenotable
CD4CD8$Cell_pheno <- Idents(CD4CD8)

# CD4CD8@meta.data %<>%  left_join(phenotable, by = 'integrated_snn_res.1.4') %>% `rownames<-`(CD4CD8$bc_backup)
CD4CD8@meta.data %<>%  left_join(phenotable, by = 'Cell_cluster', suffix = c('',''))  %>% 
  mutate(Cluster_pheno = factor(Cluster_pheno, levels = unique(Cluster_pheno))) %>% 
  
  `rownames<-`(CD4CD8$bc_backup,"Cell_cluster")

Feature_rast(CD4CD8)

CD4CD8@meta.data %<>%  
  mutate(Cell_pheno = str_replace(Cell_pheno, "TMRA", "Temra")) 

CD4CD8@meta.data %<>%   mutate(
  Cell_pheno = factor(Cell_pheno, levels = sort(unique(Cell_pheno)))
  
  
)

CD4CD8$Cell_pheno


Feature_rast(CD4CD8, 'Cell_pheno'
)


CD4CD8  %<>%  StashIdent(save.name = 'old_cluster') 
Idents(CD4CD8) <- CD4CD8$Cell_pheno


Feature_rast(CD4CD8, c('ident', 'old_cluster')
)

CD4CD8$Cell_pheno %>%  as.numeric() %>%  unique()


CD4CD8 

CD4CD8$cluster_no <- paste0('c', CD4CD8$Cell_pheno       %>%  as.numeric())

CD4CD8@meta.data  %<>%  mutate(cluster_no = factor(cluster_no,
                                                   levels = paste0('c', 1:15)))

Feature_rast(CD4CD8, 'cluster_no', noaxis = F, mythe = F)


CD4CD8$Cell_cluster <- paste0(CD4CD8$cluster_no, ': ', CD4CD8$Cell_pheno)

Feature_rast(CD4CD8, 'Cell_pheno', noaxis = F)

# Feature_rast(CD4CD8,'Cluster_pheno', colorset = 'gg' )

# Idents(CD4CD8) <- CD4CD8$Cell_cluster


UMAP_CD4CD8 <-(Feature_rast(CD4CD8, noaxis = F,sz = 0.3) +
                 ggtitle('T cells from Pulm tissue and lymph nodes')+
                 # scale_color_manual(labels = levels(CD4CD8$Cluster_pheno), values = umap.colors)+
                 guides(color = guide_legend(ncol = 2, title = 'phenotypes', override.aes = list(size = 1.5)))) %T>% 
  figsave('CD4CD8umap_CD4CD8_with_phenotypes.pdf',  180, 100  )



# change L to LN and P to Lu, pulm to Lung , LN to LLN

CD4CD8@meta.data$Cell_pheno  %<>% str_replace("_P", "_Lu")
CD4CD8@meta.data$Cell_pheno  %<>% str_replace("_Lu", "_LG")

Feature_rast(CD4CD8, "Cell_pheno")
CD4CD8@meta.data$Cell_pheno  %<>% str_replace("_L$", "_LN")
CD4CD8@meta.data$tissue  %<>% str_replace("Pulm", "Lung")
CD4CD8@meta.data$tissue  %<>% str_replace("LN", "LLN")

CD4CD8@meta.data %<>%   mutate(
  Cell_pheno = factor(Cell_pheno, levels = sort(unique(Cell_pheno)))
  
  
)

Feature_rast(CD4CD8, "Cell_pheno")
Idents(CD4CD8) <- CD4CD8$Cell_pheno

saveRDS(CD4CD8, CD4CD8RDS)

# !!!!!!!!!!!since here we adjusted and renamed the clustering, some of scripts above from line 432 to 711 need to be run again






# # distribution between tissue and CD4CD9 --------------------------------


library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[3:9], brewer.pal(9,'Reds')[3:9])

# tissue

CD4CD8$ID <- paste0(CD4CD8$tissue,'_',CD4CD8$patient)
CD4CD8$ID
#calculate tissue composition table

comp_tissue_CD4CD8 <-CD4CD8@meta.data %>% group_by(tissue, ID, Cell_pheno, CD4CD8) %>%  
  summarise(n = n() ) %>% mutate(n = case_when(tissue == 'LLN'~ -n,
                                               tissue == 'Lung' ~ n)) %>% group_by(ID,CD4CD8) %>% 
  mutate(percent = case_when(tissue == 'LLN'~ -(n/sum(n)*100),
                             tissue == 'Lung' ~ n/sum(n)*100  )) %>% as.data.frame() 

comp_tissue_CD4CD8

# tissue composition_before cluster adjustment
cl_comp_flow <- ggplot(comp_tissue_CD4CD8, aes(y = n, x = Cell_pheno, fill = ID, color = ID,
                                     stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = ID_cl  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number of cluster per donor')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,linewidth = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+facet_wrap(~CD4CD8,ncol = 1)+
  NULL
cl_comp_flow 




# CD4CD8

CD4CD8$ID_CD <- paste0(CD4CD8$CD4CD8,'_',CD4CD8$patient)



comp_CD <-CD4CD8@meta.data %>% filter(CD4CD8 %in% c('CD4', 'CD8') ) %>% 
  group_by(CD4CD8, patient, Cell_cluster) %>%  
  summarise(n = n() ) %>%group_by(patient, Cell_cluster)  %>% 
  mutate(percent = case_when(CD4CD8 == 'CD8'~ -(n/sum(n)*100),
                             CD4CD8 == 'CD4' ~ n/sum(n)*100  ))  %>%
  mutate(n = case_when(CD4CD8 == 'CD8'~ -n,CD4CD8 == 'CD4' ~ n) ) %>% ungroup() %>%  tidyr::complete(Cell_cluster, patient,CD4CD8,fill = list(n = 0, percent = 0))%>%group_by(patient, Cell_cluster)  %>% mutate(ID_CD = paste0(CD4CD8, '_', patient), total = sum(abs(n))) %>% filter(total>0)
# comp_CD %>% filter(seurat_clusters == 17) %>% as.data.frame()


comp_CD
library('tidyr')
library('broom')

t_test()
comp_CD %>% group_by(seurat_clusters, patient) %>%  
  t_test(formula = abs(percent) ~ CD4CD8)
cl_comp_flow_CD <- ggplot(comp_CD, aes(y = n, x = Cell_cluster,
                                       fill = ID_CD, color = ID_CD,
                                        stratum = ID_CD  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = CD_cl  )+
  scale_color_manual(values = CD_cl)+
  theme_minimal() + 
  ylab('Cell_number')+
  xlab(NULL)+
  scale_y_continuous(labels = abs, )+
  # ylim(-4000,4000)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+
  NULL

cl_comp_flow_CD

PG(list(cl_comp_flow,cl_comp_flow_CD), ncol = 1, rh = c(2, 1.3)) %T>%
  figsave('CD4CD8_clustercomposition.pdf', 130, 270)


meant_comp_CD<- comp_CD  %>% group_by(Cell_cluster, CD4CD8) %>% summarise( SD = sd(percent),percent = mean(percent))   %>% as.data.frame()
meant_comp_CD


comp_CD%>%ungroup %>%  group_by(Cell_cluster, patient)   %>%
  summarise_each(pv=funs(t.test(.[CD4CD8 == "CD4"], .[CD4CD8 == "CD8"])$p.value), vars=abs(percent))

library(ggpubr)

ggplot(meant_comp_CD, aes(y = percent, x = Cell_cluster  ))+
  geom_bar(stat = 'identity', aes(fill = CD4CD8), color='black', size=0.3, width = 0.95)+
  # geom_errorbar(color = 'black', size = 0.2,aes(ymax = percent+SD, ymin=percent-SD))+
  theme_bw(base_line_size = 0)+
  geom_point(data=comp_CD, aes(color=ID_CD, y = percent, x = Cell_cluster, group = patient), position = position_dodge(width = 0.8)  )+  
  scale_color_manual(values = c(ID_cl)  )+
  scale_y_continuous(labels = abs)+
  scale_fill_manual(values = alpha(c('blue','red' ), 0.2))+
  stat_compare_means(paired = T,method = 't.test')
comp_CD

ggplot(comp_CD, aes(y = percent, x = Cell_cluster  , group = patient))+
  
geom_point(data=comp_CD, aes(color=ID_CD, y = percent, x = Cell_cluster, group = patient), position = position_dodge(width = 0.8)  )+  
  scale_color_manual(values = c(ID_cl)  )+
  scale_y_continuous(labels = abs)+
  scale_fill_manual(values = alpha(c('blue','red' ), 0.2))+
  stat_compare_means(paired = T,method = 't.test')

comp_CD


comp_CD %>% group_by(Cell_cluster) %>% t.test(percent ~ CD4CD8, data = .)

comp_CD2<- comp_CD %>% mutate(percent = abs(percent) )  %>% ungroup() %>% group_by(Cell_cluster)  %>% mutate(nrows = NROW(CD4CD8)) %>% filter(nrows >4 ) 

pvs <- comp_CD2 %>% ungroup() %>% group_by(Cell_cluster) %>% summarise(pv = t.test(abs(percent) ~ CD4CD8, paired = T)$p.value)


pvs$p.adjust <- p.adjust(pvs$pv, method = 'fdr', n = length(pvs$pv))

pvs

comp_CD2   %<>% ungroup %>%  group_split(Cell_cluster) 

map(comp_CD2, ~ nrow(.x))

comp_CD2
comp_CD %>% mutate(percent = abs(percent) ) %>% ungroup() %>% group_by(Cell_cluster)  %>% mutate(nrows = NROW(CD4CD8))


comp_CD %>% group_by(Cell_cluster)  %>% mutate(nrow = NROW(CD4CD8))

comp_CD2%>% map(~ t.test(percent ~ CD4CD8, data = .x, paired = TRUE)$p.value)



t.test(percent ~ CD4CD8, data = comp_CD)$p.value

pv$p.value


# DEG CD4CD8 --------------------------------------------------------------

CD4CD8$CD_pheno <-   paste0(CD4CD8$CD4CD8, '_',CD4CD8$Cell_pheno)
CD4CD8$CD_pheno


CD4CD8_TEMRA <- ClusterCompare(CD4CD8, "CD4_Temra_1_P", "CD8_Temra_1_P", group.by = "CD_pheno",log2fc = 0.5, rm = "TRA|TRB|MT|RP") 
                                  
                               
CD4CD8_TEMRA$plot

CD4CD8_TRM1 <- ClusterCompare(CD4CD8, "CD4_TRM_1_P", "CD8_TRM_1_P", group.by = "CD_pheno",log2fc = 0.5, rm = "TRA|TRB|MT|RP") 
CD4CD8_TRM1$plot
CD4CD8_TRM1$table

ClusterCompare(CD4CD8, "CD8_TRM_1_P", "CD8_TRM_2_P", group.by = "CD_pheno",log2fc = 0.5, rm = "TRA|TRB|MT|RP|HSP") 

# gene signatrue scors ----------------------------------------------------




colnames(CD4CD8@meta.data)

CD4CD8$Cell_cluster %>% unique() %>% mixedsort()

CD4CD8$Cell_pheno %>% unique() %>% mixedsort()

CD4CD8$Cell_pheno %<>% factor(CD4CD8$Cell_pheno %>% unique() %>% mixedsort())

Ident(CD4CD8)

Feature_rast(CD4CD8, "Cell_pheno")



CD4CD8@meta.data <- CD4CD8@meta.data[, 1:50]

sigtable <- openxlsx::read.xlsx('abt/abd5778_Table_S3.xlsx') %>% 
  `colnames<-`(str_remove(colnames(.), '.\\(.+\\)') ) %>%  as.list() %>% 
  map(~ na.exclude(.x) %>% as.vector) %T>% print()
 
names(sigtable)
sigtable$Tissue.resident

TRM <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "CD69", "GPR25")

#caculate scores
for (i in names(sigtable)) {
  CD4CD8 <- AddModuleScore(CD4CD8, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}


CD4CD8 <- AddModuleScore(CD4CD8, features = list(TRM), name = "TRM_score", assay = 'RNA')

CD4CD8$Tissue.resident <- CD4CD8$TRM1

colnames(CD4CD8@meta.data) %<>% str_replace("(?<=\\w)1$", '')
CD4CD8@meta.data
# CD4CD8$Th1 <- CD4CD8$Th

sigtable %>% names()

Feature_rast(CD4CD8,c(names(sigtable)[c(1:9, 14,15)], 'ident', 'tissue', 'CD4CD8'), color_grd = 'grd')

ViolinPlot(CD4CD8, names(sigtable)[c(1,2,4,5,7,8,9,10,12,13, 14,15)], 
           colors = umap.colors, 
           box = T,
           jitter = T,
           x.angle = 90) %T>% 
  figsave("CD4CD8_GM_vlnplot.pdf", 300, 350)

CD4CD8@assays$GM <- NULL
GMS  <- CD4CD8@meta.data[,c(names(sigtable))] %>% as.data.frame() %>% t()
scaleGM <- scale(t(CD4CD8@meta.data[,c(names(sigtable))]))
scaleGMassay <- CreateAssayObject(data = scaleGM)
CD4CD8@assays$GM <- scaleGMassay
CD4CD8@assays$GM@key <- "GM_"

DEM <- 
  FindAllMarkers(CD4CD8, test.use = 'wilcox', assay = "GM", 
                 
                 only.pos = T )

DEM

DoHeatmap(subset(CD4CD8, downsample = 700), 
          raster =T, draw.lines = T, angle = 45,
          lines.width = 10,group.colors = umap.colors,
          
          assay = 'GM', features = names(sigtable), slot = 'data', size = gs(8)) +hmp2 + mytheme+
  theme(legend.position = 'bottom',
        legend.key.height = unit(2,'mm'))+
  guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))


names(sigtable)
CD4CD8@assays$GM@data

FetchData(subset(CD4CD8, downsample = 700), c("CD4CD8", "Cell_pheno", names(sigtable))) %>% 
  write.csv( "CD4CD8_GMS_downsample.csv")


GMS_long <- FetchData(CD4CD8, c("CD4CD8", "Cell_pheno", names(sigtable))) %>% 
  reshape2::melt(value.name = "score" ) %>% 
  mutate(variable = if_else(variable =="CD8.Cytotoxictiy","Cytotoxicity",variable ) )

head(GMS_long)
GMS_long$variable %>% unique()

ggplot(GMS_long %>% filter(variable %in% c("Effectors", "Tissue.resident", "Exhaustion", "Th17", "Cytotoxicity") & Cell_pheno != "unidentified_P" ) %>% 
         mutate(variable = factor(variable, 
                                  level = c("Effectors", "Tissue.resident", "Exhaustion", "Th17", "Cytotoxicity"))), 
       aes(x = Cell_pheno, y = score, fill = Cell_pheno))+
  geom_violin_rast(size = 0.1) +facet_grid(variable~ CD4CD8,scales = "free_y"  )+
  ylab("module score")+xlab(NULL)+
  geom_boxplot( alpha = 0.5, size = 0.3,  width = 0.5, outlier.alpha = 0)+
  fill_m()+theme_classic()+mytheme+
  theme(axis.text.x = element_blank())
  




saveRDS(CD4CD8, CD4CD8RDS)
# tissue residency --------------------------------------------------------

ViolinPlot(CD4CD8, c('CD49a.protein', 'CD103.protein', 'CD4.protein', 'CD45RA.protein'), assay = 'CITE', group.by = 'patient'  ,colors = umap.colors, box = T)

CD4CD8_cite <- subset(CD4CD8, patient %in% c('p27', 'p71', 'p77', 'p73'))

CD4CD8@assays$CITE %>%  rownames()

Feature_rast(CD4CD8_cite, c('CD49a.protein', 'CD103.protein', 'KLRG1.protein', 'ident'), assay = 'CITE')

ViolinPlot(CD4CD8_cite, c('CD49a.protein', 'CD103.protein'), assay = 'CITE', colors = umap.colors, box = T)

Feature_rast(CD4CD8, c('ITGA1', 'ITGAE', 'ZNF683'))

CD4CD8_cite <- AddModuleScore(CD4CD8_cite,assays = 'CITE', features = c('CD49a.protein', 'CD103.protein'))
CD4CD8_cite@assays$CITE


CD4CD8_103_49a_CITE <-(Feature_rast(CD4CD8_cite, 
             d1 ='CD103.protein', d2 =  'CD49a.protein', facets = 'CD4CD8', noaxis = F,
             axis.number = T)+
  geom_hline(yintercept = 0.8)+
  geom_vline(xintercept = 1.2)) %T>% figsave('CD4CD8_103_49a_CITE.pdf', 200,100) 

# define Trm
TRM_BC <- CD4CD8_cite %>%  subset( CD103.protein >= 1.2 &  CD49a.protein >= 0.8) %>% colnames()

CD4CD8_cite@meta.data  %<>% 
  mutate(TRM = case_when(bc_backup %in% TRM_BC & CD4CD8 == 'CD4' ~ 'CD4 Trm' ,
                         bc_backup %in% TRM_BC & CD4CD8 == 'CD8' ~ 'CD8 Trm'                           
                                                    ))

Feature_rast(CD4CD8_cite, 'TRM', facets = 'tissue') 


# percentage of Trm in each cluster,#

CD4CD8_cite@meta.data %>% group_by(seurat_clusters) %>% 
  count(TRM) %>%  group_by(seurat_clusters) %>% 
  mutate(percent = n/sum(n)*100) %>% 
  ggplot(aes(x = seurat_clusters, y = percent, fill = TRM, group = TRM))+
  geom_col()+
  fill_m()



CD4CD8@meta.data %>% group_by(tissue, patient) %>%
  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)


Feature_rast(CD4CD8, c('RORC','CCR6', 'CD40LG', 'D'))







# # distribution between tissue and CD4CD9 --------------------------------


library(ggalluvial)
library(RColorBrewer)
RColorBrewer::display.brewer.all()
ID_cl <- c(brewer.pal(9,'Blues')[3:9], brewer.pal(9,'Reds')[3:9])
CD_cl <- c(brewer.pal(9,'RdPu')[3:9], brewer.pal(9,'Greens')[3:9])

# tissue

CD4CD8$ID <- paste0(CD4CD8$tissue,'_',CD4CD8$patient)
Feature_rast(CD4CD8, 'ID', colorset = CD_cl, do.label = F, noaxis = F)
Umap_donor_tissue <- Feature_rast(CD4CD8, 'ID', colorset = ID_cl, do.label = F, noaxis = F)+ ggtitle(NULL)+NoLegend()
Umap_donor_tissue
#calculate tissue composition table

# CD4CD8$Cell_cluster <- Idents(CD4CD8)

comp_tissue <-CD4CD8@meta.data %>% group_by(tissue, ID, Cell_pheno, cluster_no) %>%  
  summarise(n = n() ) %>% mutate(n = case_when(tissue == 'LLN'~ n,
                                               tissue == 'Lung' ~ -n)) %>% group_by(ID) %>% 
  mutate(percent = case_when(tissue == 'LLN'~ (n/sum(n)*100),
                             tissue == 'Lung' ~ -n/sum(n)*100  )) %>% as.data.frame() 

comp_tissue


cl_comp_flow <- ggplot(comp_tissue, aes(y = n, x = Cell_pheno, fill = ID, color = ID,
                                        stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = ID_cl  )+
  scale_color_manual(values = ID_cl)+
  theme_bw() + 
  ylab('Cell number of cluster per donor')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL,
                             keywidth = unit(2, 'mm'),
                             byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', 
        axis.text.x = element_text(angle = 315, size = 8, vjust = -0.6),
        
        legend.key.height  = unit(3, 'mm'),
        legend.key.width   = unit(5,'mm'),
        axis.line = element_blank())+
  NULL
cl_comp_flow


UMAPsCD4CD8 <- (Feature_rast(CD4CD8, sz = 0.3,
                            labelsize = 6,
                            noaxis = F
                            
                            )+ggtitle('abT cells from lung and luLN'))
UMAPsCD4CD8

figtissue <- PG(list(UMAPsCD4CD8, cl_comp_flow), rw = c(1, 1.1))
figtissue

figsave(figtissue, 'CD4CD8_UMAP_and_tissuedistribution.pdf', 200, 85)

# CD4CD8

CD4CD8$ID_CD <- paste0(CD4CD8$CD4CD8,'_',CD4CD8$patient)

Umap_donor_CD <- Feature_rast(subset(CD4CD8, CD4CD8 %in% c('CD4', 'CD8')), 'ID_CD', colorset = CD_cl, do.label = F, noaxis = F)+ ggtitle(NULL)+NoLegend() 
Umap_donor_CD

comp_CD <-CD4CD8@meta.data %>% filter(CD4CD8 %in% c('CD4', 'CD8') ) %>% 
  group_by(CD4CD8, patient, Cell_pheno, cluster_no) %>%  
  summarise(n = n() ) %>%group_by(patient, Cell_pheno,cluster_no)  %>% 
  mutate(percent = case_when(CD4CD8 == 'CD8'~ -(n/sum(n)*100),
                             CD4CD8 == 'CD4' ~ n/sum(n)*100  ))  %>%
  mutate(n = case_when(CD4CD8 == 'CD8'~ -n,CD4CD8 == 'CD4' ~ n) ) %>% ungroup() %>%  tidyr::complete(Cell_pheno,cluster_no, patient,CD4CD8,fill = list(n = 0, percent = 0))%>%group_by(patient, Cell_pheno,cluster_no)  %>% mutate(ID_CD = paste0(CD4CD8, '_', patient), total = sum(abs(n))) %>% filter(total>0)


comp_CD
library('tidyr')
library('broom')
library(rstatix)



meant_comp_CD<- comp_CD  %>% group_by(Cell_pheno, cluster_no,CD4CD8) %>% summarise( SD = sd(percent),percent = mean(percent))   %>% as.data.frame()
meant_comp_CD

meant_comp_CD

library(ggpubr)

cl_comp_flow_CD<-ggplot(meant_comp_CD, aes(y = percent, x = Cell_pheno  ))+
  geom_bar(stat = 'identity', aes(fill = CD4CD8), 
           color='black', size=0.3, width = 0.95)+
  # geom_errorbar(color = 'black', size = 0.2,aes(ymax = percent+SD, ymin=percent-SD))+
  theme_bw(base_line_size = 0)+
  geom_point(data=comp_CD, aes(color=ID_CD, y = percent, x = Cell_pheno, group = patient), position = position_dodge(width = 0.65) , size = 0.7 )+  
  scale_color_manual(values = c(CD_cl), labels = rep(sort(unique(CD4CD8$patient)),2)  )+
  scale_y_continuous(labels = abs ,breaks = seq(-100, 100, 20))+
  scale_fill_manual(values = alpha(c('purple','green' ), 0.2))+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 315, size = 8, vjust = -0.6),
  )+
  guides(  fill = guide_legend(nrow = 2, title = NULL),
    color = guide_legend(nrow = 2,  byrow = T, label.position = 'top', title = 'patient',override.aes = list(size = 2))
       )
cl_comp_flow_CD
figCD <- PG(list(Umap_donor_CD, cl_comp_flow_CD), rw = c(1, 1.3))
figCD%T>% figsave('CD4CD8_CD4CD8_distribution.pdf',210,110)

compositions <- PG(list(figtissue, figCD), ncol = 1) %T>% figsave('tissue_donor_CD4_8.pdf',200,220)

comp_CD2<- comp_CD %>% mutate(percent = abs(percent) )  %>% ungroup() %>% group_by(Cell_cluster)  %>% mutate(nrows = NROW(CD4CD8)) %>% filter(nrows >4 ) 

pvs <- comp_CD %>% ungroup() %>% group_by(Cell_pheno) %>% summarise(pv = t.test(abs(percent) ~ CD4CD8, paired = T)$p.value)

pvs$fdr <- p.adjust(pvs$pv, method = 'fdr', n = length(pvs$pv))

pvs$p.adjust <- p.adjust(pvs$pv, method = 'bonferroni', n = length(pvs$pv))
pvs

  






# DEGs --------------------------------------------------------------------
DefaultAssay(CD4CD8) <- 'RNA'


CD4CD8_DEGs <- 
  FindAllMarkers(CD4CD8, test.use = 'bimod', min.pct = 0.10,   only.pos = T )%>%
  filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))
CD4CD8_DEGs %>% filter(avg_log2FC >0) %>% count(cluster)

CD4CD8_DEGs  %<>% filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) 
write.xlsx(CD4CD8_DEGs, 'PulmabT_DEGs.xls')

top10deg <- CD4CD8_DEGs %>%
  # filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(5, avg_log2FC )


top10DEG_heat <- 
  DoHeatmap(subset(CD4CD8, downsample=500),features = top10deg$gene,raster = T,  assay = "RNA",
            group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>% figsave('top10DEG.pdf', 190, 370)
top10DEG_heat





 
DoHeatmap(subset(CD4CD8, Cell_cluster %in% paste0('C',13:18)), features = top15_TrmMarker$gene) %>% heat_theme()

# Tnv markers
Tnvmarkers <-FindAllMarkers(subset(CD4CD8, Cell_cluster %in% c(7,8,10)),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA') 
top15_Tnvmarkers <- Tnvmarkers %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(15, avg_log2FC ) %>% as.data.frame()

DoHeatmap(subset(CD4CD8, Cell_cluster %in% c(7,8,10)), features = top15_Tnvmarkers$gene) %>% heat_theme()

TFs <- readLines('HumanTFs_1639.txt')
TFs_T <- intersect(TFs, rownames(CD4CD8))


grep("GZM", rownames(CD4CD8), value = T)


Cytlist <- read.table("Cytokine_gene_list.txt", header = T) %>% pull(Symbol) %>% intersect(rownames(CD4CD8)) %>% c(grep("GZM", rownames(CD4CD8), value = T))

Cytlist


cd4cd828 <-ClusterCompare(CD4CD8, 'TRM_1_LG', 'TRM_3_LG', features = TFs_T)

ClusterCompare(CD4CD8, 'TRM_1_LG', 'TRM_3_LG')
cd4cd828$plot

ClusterCompare(CD4CD8, 'TRM_1_LG', 'TRM_3_LG', features = Cytlist,genetoshow = 20)


ClusterCompare(CD4CD8, 'TRM_1_LG', 'Th17_M', features = Cytlist,genetoshow = 20)



intersect(P2P8_TFs$table$gene[1:63], cd4cd828$table$gene[1:26])  %>% 
  Feature_density(CD4CD8, .)
intersect(P2P8_TFs$table$gene[1:63], cd4cd828$table$gene[1:26])  %>% 
  Feature_density(GDTlung_s, .)

Feature_rast(GDTlung_s, c('PDCD1', 'CTLA4'), sz= 0.2)
cd4cd828$table


inst_cls <- c('TRM_1_P',  'TRM_2_P', 'TRM_3_P','Temra_1_P', 'Temra_2_P', 'TEM_M')

TrmMarker <- FindAllMarkers(subset(CD4CD8, Cell_pheno %in% inst_cls),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA')%>%
  filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))

TrmMarker

top15_TrmMarker <- TrmMarker %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC ) %>% as.data.frame()



DoHeatmap(subset(CD4CD8, Cell_pheno %in% inst_cls, downsample=500),features = top15_TrmMarker$gene,raster = T, 
          group.colors = umap.colors,size = gs(8))%>%heat_theme() 


# psedubulk 


pseudobuklCD4CD8 <- AggregateExpression(CD4CD8, assays = "RNA", return.seurat = T, group.by = c(  "Cell_pheno", "patient"))

colnames(pseudobuklCD4CD8)

pseudobuklCD4CD8$orig.ident

bulkDEGCD4CD8 <- FindAllMarkers(pseudobuklCD4CD8, test.use = "DESeq2", logfc.threshold = 0.1)

cytokine_DEG <-  bulkDEGCD4CD8 %>% dplyr::filter(gene  %in% Cytlist & avg_log2FC > 0.5)
cytokine_DEG

top5deg_CD4CD8_bulk <- bulkDEGCD4CD8 %>%
  arrange(cluster, desc(avg_log2FC))%>%
  filter(!grepl('^RP|^MT|^TR', gene)) %>%
  group_by(cluster) %>% top_n(5, avg_log2FC )


top5deg_CD4CD8_bulk


(DoHeatmap(pseudobuklCD4CD8,features = top5deg_CD4CD8_bulk$gene,raster = T,  
          assay = "RNA", l
          group.colors = umap.colors,size = gs(8) )  )%>% heat_theme()

ClusterCompare(pseudobuklCD4CD8, "TRM-1-LG", "TRM-3-LG", test = "DESeq2", rm = "MT|RP|TR")


ClusterCompare(pseudobuklCD4CD8, "Temra-1-LG", "TRM-3-LG", test = "DESeq2", rm = "MT|RP|TR")


DoHeatmap()

# DEG TRM TEMRA TEM -------------------------------------------------------

trmcl <- c(
  
   "TEM_M",
  "Temra_1_LG",
  "TRM_1_LG",
  "TRM_2_LG",
  "TRM_3_LG"
  # "Temra_2_P"
  
)




TRM_DEGs <- 
  FindAllMarkers(CD4CD8 %>% subset(Cell_pheno %in% trmcl), test.use = 'bimod', min.pct = 0.10,   only.pos = FALSE )%>%
  filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))


TRM_DEGs_TF <- 
  FindAllMarkers(CD4CD8 , test.use = 'bimod', min.pct = 0.10,   only.pos = T, features = TFs_T )%>%
  filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))

TRM_DEGs_TF %>%  filter(cluster == "TRM_1_P")


# TRM_DEGs  %<>% filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) 
# write.xlsx(CD4CD8_DEGs, 'PulmabT_DEGs.xls')

top10deg_TRM <- TRM_DEGs %>%
  filter(!grepl('^RP|^MT|^TR|^HSP|^HIST', gene)) %>%
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC )


top10TF_TRM <- TRM_DEGs_TF %>%
  filter(!grepl('^RP|^MT|^TR|^HSP|^HIST', gene) &
         cluster %in% trmcl
         ) %>% group_by(gene) %>% 
  top_n(1, avg_log2FC) %>% 
  
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC )


view(TRM_DEGs_TF)
heatDOT_top10_TRM <-( DotPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl),dot.scale = 3.5,
                               features = rev(unique(top10deg_TRM$gene)))+mytheme+heattheme+
                         
                         theme(text = element_text(size = 8), axis.text.y = element_text(size = 8),
                               axis.line.y.right = element_line(),
                               axis.text.x = element_text(size = 8, angle= 90),
                               legend.box.margin = margin(5,0,0,15,unit = 'mm'),
                               legend.box = "horizontal",legend.position = 'bottom',
                               axis.title = element_blank())+coord_flip()+
                         scale_y_discrete(position = 'right')+
                         scale_x_discrete(position = 'top')+
                         xlab(NULL)+ylab(NULL)+
                         scale_color_gradient2(low = '#003399', mid = 'lightgreen',  high = "#990000")+
                         guides(
                           color = guide_colorbar(title.position = 'top',direction = 'horizontal',
                           ),
                           size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom'))) 
view(top10TF_TRM)
heatDOT_top10_TRM
heatDOT_top10_TF_TRM <-( DotPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl),dot.scale = 3.5,
                              features = rev(unique(top10TF_TRM$gene)))+mytheme+heattheme+
                        
                        theme(text = element_text(size = 8), axis.text.y = element_text(size = 8),
                              axis.line.y.right = element_line(),
                              axis.text.x = element_text(size = 8, angle= 90),
                              legend.box.margin = margin(5,0,0,15,unit = 'mm'),
                              legend.box = "horizontal",legend.position = 'bottom',
                              axis.title = element_blank())+coord_flip()+
                        scale_y_discrete(position = 'right')+
                        scale_x_discrete(position = 'top')+
                        xlab(NULL)+ylab(NULL)+
                          viridis::scale_color_viridis(discrete = F, option ="A")+
                        guides(
                          color = guide_colorbar(title.position = 'top',direction = 'horizontal',
                          ),
                          size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom')))  %T>% print()



DEGs_TRM_Temra <- FindAllMarkers(CD4CD8 %>% subset(subset =  Cell_pheno %in% c("TRM_1_P", "TRM_2_P", "TRM_3_P", "Temra_1_P", "Temra_2_P", "TEM_M")  ),   
                                 only.pos = T,assay = 'RNA' ) %>%     
  filter(p_val_adj < 0.05) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))
DEGs_TRM_Temra %<>%  filter(!grepl("TRA|TRB|TRG|MT|RP|HIST", gene) )

avggenel <- DEGs_TRM_Temra%>%
  pull(gene) %>% unique
DEGs_TRM_Temra$cluster
avggenel_2 <- CD4CD8_DEGs %>%  filter(cluster %in%  c("TRM_1_P", "TRM_2_P", "TRM_3_P", "Temra_1_P", "Temra_2_P", "TEM_M")) %>%   pull(gene) %>% unique


length(avggenel_2)

AVG_GENE <- AverageExpression(object = CD4CD8 %>% subset(subset =  Cell_pheno %in% c("TRM_1_P", "TRM_2_P", "TRM_3_P", "Temra_1_P", "Temra_2_P", "TEM_M")  ),  verbose = T, assays = 'RNA',
                              slot = 'counts', features = avggenel_2)
AVG_GENE_scaled <- t(AVG_GENE$RNA ) %>% 
  scale(center = T, scale = T) %>% as.data.frame()%>% t()
library(reticulate)
library(umap)
umap::umap()

Umap_gene <- umap(AVG_GENE_scaled, preserve.seed = 1 , method = "umap-learn") 
Umap_gene_gg <- Umap_gene$layout %>% as.data.frame() %>%
  `colnames<-`(c('UMAP_1', 'UMAP_2')) %>% rownames_to_column("gene")
nrow(Umap_gene_gg)

Feature_rast(Umap_gene_gg, g = "gene", do.label = F, othertheme = NoLegend(), colorset = rainbow(1000))
head(AVG_GENE_scaled)

write.xlsx(AVG_GENE_scaled, "AVG_GENE_scaled.xlsx")
write.csv(AVG_GENE_scaled, "AVG_GENE_scaled.csv")




DEG_TRM1_TRM3 <-  ClusterCompare(CD4CD8, 'TRM_1_P', "TRM_3_P")

DEG_CD4_CD8 <-  ClusterCompare(CD4CD8, 'CD4', "CD8", group.by = "CD4CD8")

TRM13unique <-  DEG_TRM1_TRM3$table %>% filter(!(gene %in% DEG_CD4_CD8$table$gene))

TRM13unique %>% filter(gene %in% TFs_T)


DoHeatmap(CD4CD8 %>% subset(Cell_pheno %in% c("TRM_1_P", "TRM_3_P")), features = TRM13unique$gene)


TRM13unique$gene

DEG_CD4_CD8$table$gene
# DEG CD4CD8 --------------------------------------------------------------

CD4CD8$CD_pheno <-   paste0(CD4CD8$CD4CD8, '_',CD4CD8$Cell_pheno)

Feature_rast(subset(CD4CD8, CD4CD8 %in% c('CD4', 'CD8')), 'CD_pheno', colorset = rainbow(40))



ClusterCompare(CD4CD8, 'CD4_Temra_1_P', 'CD8_Temra_1_P', group.by = 'CD_pheno')

# SCENIC result -----------------------------------------------------------
library(SCopeLoomR)
library(SCENIC)


loom <- open_loom('ABTlung_pyscenic/lungABT_output_tracks_nopara.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

AUCmat <- AUCell::getAUC(regulonAUC) 

rownames(AUCmat)
rownames(AUCmat)   %<>% gsub("\\(\\+\\)", "_REG", .)
CD4CD8[['AUC']] <- CreateAssayObject(data = AUCmat)

CD4CD8 <- ScaleData(CD4CD8, assay = 'AUC')

trmcl <- c(
  
  "TEM_M",
  "Th17_M",
  "Temra_1_LG",
  "Temra_2_LG",
  "TRM_1_LG",
  "TRM_2_LG",
  "TRM_3_LG"
  
)


DE_48_Reglon <- FindAllMarkers(CD4CD8 %>% subset(Cell_pheno %in% trmcl),, only.pos = T, assay = "AUC",min.pct = 0.25, test.use =  "bimod",slot = "scale.data")



DE_48_Reglon
DE_48_Reglon %>% filter(cluster == "TRM_1_LG")

top5_48_reg <-  DE_48_Reglon %>% group_by(cluster) %>%  filter(cluster  %in% trmcl) %>% top_n(5, avg_diff) %>%  pull(gene)


DoHeatmap(CD4CD8, features = top5_48_reg, assay = "AUC") %>% heat_theme()
DoHeatmap(CD4CD8 %>% subset(downsample = 500,Cell_pheno %in% trmcl), features = top5_48_reg, assay = "AUC") %>% heat_theme() %T>% figsave("top_reg_CD4CD8.pdf", 400, 400)


adjabt <- vroom::vroom("ABTlung_pyscenic/adj_abt.csv")

top20_ZNF559 <- adjabt %>% filter(TF == "ZNF559") %>% arrange(desc(importance)) %>% top_n(20, importance) %>%  pull(target)
Feature_rast(CD4CD8, top20_ZNF559)

adjabt %>% arrange(desc(importance))
adjabt %>% filter(TF == "ZNF683") %>% arrange(desc(importance))

top20_ZNF683 <- adjabt %>% filter(TF == "ZNF683") %>% arrange(desc(importance)) %>% top_n(20, importance) %>%  pull(target)
Feature_rast(CD4CD8, top20_ZNF683)


Motif <- vroom::vroom("GDTlung_pyscenic/auxilliaries/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

Motif <- vroom::vroom("GDTlung_pyscenic/auxilliaries/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")



glimpse(Motif)
Motif  %>% filter(gene_name == "ZNF559") %>% view()
view(Motif)


adjabt$TF %>% unique() 

regulon <- vroom::vroom("ABTlung_pyscenic/reg_2.csv")

regulon %>% filter(...1 == "ZNF559")
AUCcells <- vroom::vroom("ABTlung_pyscenic/lungABT_output_tracks_AUC_0.001.csv")  %>% `colnames<-`(
  gsub("\\(\\+\\)", "_REG", colnames(.) )  
)
  
colnames(AUCcells) 



colnames(AUCcells) 

head(regulons)
colnames(AUCcells) 
dim(AUCcells)


dim(UMAPmeta)

UMAPmeta$bc_backup == AUCcells$Cell

UMAPmeta <-  FetchData(CD4CD8, c("CD4CD8", "UMAP_1", "UMAP_2", "Cell_pheno", "bc_backup")) %>% cbind(AUCcells)
UMAPmeta


Feature_rast(CD4CD8, c( "EOMES-REG", "MAF-REG","TBX21-REG",   "RORC-REG", "RORA-REG"), colorgrd = "grd2", assay = "AUC")






ggplot(UMAPmeta, aes(x = Cell_pheno, y = EOMES_REG, fill = Cell_pheno))+geom_violin()+geom_boxplot()

data_frame(UMAPmeta, AUCcells)
colnames(UMAPmeta)

UMAPmeta  %<>% data_frame(AUCcells) %>% as.data.frame()

Feature_rast(UMAPmeta, "FOXP3_REG")


head(UMAPmeta)

class(UMAPmeta)

dim(UMAPmeta)
class(UMAPmeta %>% as.data.frame())

Feature_rast(UMAPmeta, "RUNX3_REG")

Feature_rast(CD4CD8, "RUNX3-REG", assay = "AUC")



UMAPmeta$`GATA3(+)`

regulons %>% filter(X == "GATA3") %>% view()

CD4CD8@assays$REG <- NULL

# scaleGM <- scale(t(CD4CD8@meta.data[,c(names(sigtable))]))
REGassay <- CreateAssayObject(data = AUCcells %>% column_to_rownames("Cell") %>% t())
CD4CD8@assays$REG <- REGassay
CD4CD8@assays$REG@key <- "REG_"

DEGreg <- FindAllMarkers(CD4CD8, assay = "REG", only.pos = T,min.pct = 0,
                         logfc.threshold = 0)
DEGreg %>% view()
DEGreg   %<>% filter(p_val_adj < 0.01)
view(DEGreg)

ClusterCompare(CD4CD8, "TRM_1_LG", "TRM_3_LG", assay = "AUC", log2fc = 0.01)

CD4CD8@assays$AUC


allregs_TRM1_3 <-  ClusterCompare(CD4CD8, "TRM_1_LG", "TRM_3_LG", assay = "AUC", log2fc = 0.00,p_cutoff = 0.8)


ggplot(allregs_TRM1_3$table, aes(x = avg_log2FC , y = -log10(p_val_adj ), color = abs(avg_log2FC)*-log10(p_val_adj )  ))+ geom_point()
regulons <- c("RORA-REG", "MAF-REG", "RUNX2-REG", 

              "EOMES-REG", "TBX21-REG", "RUNX3-REG")



ViolinPlot(CD4CD8 %>%  subset(Cell_pheno %in% c("TRM_1_LG" , "TRM_3_LG")), box = T,  sz = 0.2,
           x.angle = 330, ylabtext = "value",  ncol = 4, g = regulons,
           othertheme = list(stat_compare_means(  paired = F, method = 'wilcox.test',label = "p.format" ),
                             theme(axis.text.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 6))),
           assay = "AUC", colors = umap.colors[c(12,14)]) %T>%  print()


# # GSEA genesymobl instead --------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)



library(msigdbr)

Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") 

Hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") 

GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'BP')

Mc2 <-  msigdbr::msigdbr(species = "Homo sapiens", category = "C2") 





ALL_msigdb_G  <- rbind(
  # c7
  msigdbr::msigdbr(species = "Homo sapiens", category = "C7"), 
  # hallmarker
                       msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  # C2
                       msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                                        # GOBP
                                        subcategory = 'CP:KEGG'),msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')) %>% 
  dplyr::select(gs_name, gene_symbol)


TEMRA_TRM_3_G <- Genelist_generator(CD4CD8, id1 = "Temra_1_LG", id2 = "TRM_3_LG")






TEMRA_TRM_3_msigdb_G <- GSEA(geneList = TEMRA_TRM_3_G, TERM2GENE=ALL_msigdb_G, verbose = T,
                           pvalueCutoff = 0.05, pAdjustMethod = "BH") 
TEMRA_TRM_3_msigdb_G@result  <- TEMRA_TRM_3_msigdb_G@result %>% arrange(desc(NES)) 

view(TEMRA_TRM_3_msigdb_G@result )


TRM1_3 <- Genelist_generator(CD4CD8, 'TRM_1_LG','TRM_3_LG')
TRM1_3_Msigdb_GSEAc7 <- GSEA(geneList = TRM1_3, TERM2GENE=Mc7, verbose = T,
                             # minGSSize    = 15,
                             pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
TRM1_3_Msigdb_GSEAc7@result  <- TRM1_3_Msigdb_GSEAc7@result %>% arrange(desc(NES)) 
view(TRM1_3_Msigdb_GSEAc7@result)


TRM1_3_Msigdb_all_GSEA_G <- GSEA(geneList = TRM1_3, TERM2GENE=ALL_msigdb_G,  
                               pvalueCutoff = 0.05, pAdjustMethod = "BH") 
TRM1_3_Msigdb_all_GSEA_G@result  <- TRM1_3_Msigdb_all_GSEA_G@result %>% arrange(desc(NES)) 

view(TRM1_3_Msigdb_all_GSEA_G@result)




# TCR data frome processing & integration  --------------------------------------------------------------------

# read raw TCR files 
CD4CD8$bc_backup <- colnames(CD4CD8)
colnames(CD4CD8) %>% str_extract('_.') %>% unique()
CD4CD8$orig.ident %>% unique()
CD4CD8@meta.data %>% select(bc_backup, orig.ident) %>% group_by(orig.ident) %>% slice(1)

# [1] "_1" "_2" "_3" "_4" "_5" "_6"
proj
# [1] "p32"     "p45"     "p25_CD4" "p25_CD8" "p27"     "p71" 
pt


TCRdirs <- list.dirs(path = 'raw') %>% str_subset('VDJ.+outs$' )
# remove patient 31
TCRdirs <- TCRdirs[c(1:5,7:9)]




# note the number of the y argument is dependant the suffix of barcodes in the scRNAseq dataset

TCRs <- map2(TCRdirs, c(1,2,5:8, 3,4), ~   (
  do.call(rbind, lapply( list.files(path = .x, full.names = T,
                                    pattern = 'filtered_contig_annotations')  , read.csv))
    )  %>% 
               filter(productive == 'True'& is_cell == 'True') %>%
               dplyr::select(c(1, 5:10, 13,14))  %>%
               dplyr::mutate( patient = pt[[.y]],  
                              bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% 
  set_names(proj) %>% reduce(.f = rbind)

TCRs %>%  group_by(patient) %>%  select(patient, bc_backup) %>% slice(1)
TCRs %>%  group_by(patient) %>%  select(patient, bc_backup) %>%  filter(patient == 'p25') %>% 
  pull(bc_backup) %>% str_extract('_.') %>% unique()

TCRs %>% head()
  




# genereate TRA and TRB table

TRAs <- TCRs %>% filter(chain == 'TRA')%>% distinct(bc_backup, .keep_all = T) %>% 
  mutate(v_gene = str_remove(v_gene,'/DV\\d')) %>%
  rename_at(vars(-bc_backup,-patient), funs(sub('$','_TRA',.)))
TRBs <- TCRs %>% filter(chain == 'TRB')%>% distinct(bc_backup, .keep_all = T) %>% 
  rename_at(vars(-bc_backup, -patient), funs(sub('$','_TRB',.)))

#  then put it back
TCRs_paired <- full_join(TRAs,TRBs, by=c('bc_backup','patient')) %>%
  mutate(paired = case_when(!is.na(v_gene_TRA) & !is.na(v_gene_TRB) ~ paste0(str_remove(v_gene_TRA, 'TR'),
                                                                             ' ',str_remove(v_gene_TRB, 'TR')))) %>%
  mutate(cdr3_paired =   case_when(!is.na(v_gene_TRA) & !is.na(v_gene_TRB) ~ 
                       paste(paired,cdr3_TRA,cdr3_TRB)),
         paired_sp = case_when(!is.na(paired) ~ 'Paired TCR')) %>%
  filter(!duplicated(bc_backup))


TCRs_paired$paired_sp %>% table
##TCR CDR3 frequency
cdr3TRA_freq <- TCRs_paired %>% filter(chain_TRA == 'TRA' & cdr3_TRA != 'None') %>% group_by(patient)%>%
  dplyr::count(cdr3_TRA = cdr3_TRA) %>% arrange(desc(n)) %>% dplyr::rename(cdr3_TRA_freq = 'n') %>%
  mutate(cdr3_TRA_perc =  cdr3_TRA_freq/sum(cdr3_TRA_freq)*100) 
cdr3TRA_freq$cdr3_TRA_perc %>% sum()

cdr3TRB_freq <- TCRs_paired %>% filter(chain_TRB == 'TRB' & cdr3_TRB != 'None') %>% group_by(patient)%>%
  dplyr::count(cdr3_TRB = cdr3_TRB) %>% arrange(desc(n)) %>% dplyr::rename(cdr3_TRB_freq = 'n')%>%
  mutate(cdr3_TRB_perc =  cdr3_TRB_freq/sum(cdr3_TRB_freq)*100)
cdr3TRB_freq %>% top_n(5, cdr3_TRB_perc)
cdr3Paired_freq <- TCRs_paired %>% filter(!is.na(cdr3_paired)) %>% group_by(patient)%>%
  dplyr::count(cdr3_paired = cdr3_paired) %>% arrange(desc(n)) %>% 
  dplyr::rename(cdr3_paired_freq = 'n')%>%
  mutate(cdr3_paired_perc =  cdr3_paired_freq/sum(cdr3_paired_freq)*100)
cdr3Paired_freq$patient %>% unique()


TCRs_paired %<>%  left_join(cdr3TRA_freq, by = c('cdr3_TRA','patient')) %>% 
  left_join(cdr3TRB_freq, by = c('cdr3_TRB','patient')) %>%
  left_join(cdr3Paired_freq, by = c('cdr3_paired','patient'))

TCRs_paired %>% head()


TCRs_paired %>%  group_by(patient) %>%  select(patient, bc_backup) %>% slice(1)


CD4CD8@meta.data  %>%  group_by(patient) %>%  select(patient, bc_backup) %>% slice(1)

TCRs_paired %>%  group_by(patient)  %>% summarise(n = n())
saveRDS(TCRs_paired, 'abTCR_new.RDS')

# join the data 
CD4CD8$bc_backup <- rownames(CD4CD8@meta.data)


# in case there is duplicated data, remove old data at first

CD4CD8@meta.data  %<>%   select_at(.vars = vars(-contains(c('TRA','TRB','paired', '.x', '.y'))))
CD4CD8@meta.data  %<>%   select_at(.vars = vars(-contains(c('TRA','TRB','paired', '.x', '.y'))))

CD4CD8@meta.data  %<>%   select_at(.vars = vars(-contains(c('Age_'))))
colnames(CD4CD8@meta.data)

CD4CD8@meta.data %<>% left_join(TCRs_paired, by =c('bc_backup', 'patient'), suffix = c('', '') ) %>% `rownames<-`(CD4CD8$bc_backup)




Feature_rast(CD4CD8, 'cdr3_paired_perc',sz = 0.5)
Feature_rast(CD4CD8, c('cdr3_TRB_perc', 'tissue'), color_grd = 'grd', facets = 'patient')



saveRDS(CD4CD8, CD4CD8RDS)
dim(CD4CD8@assays$RNA@counts)

# TCR analysis -----------------------------------------------------------
# mapping rate
# to see how many T cells are paired with a TCR 


CD4CD8@meta.data %<>% mutate(TCR_summary = case_when(!is.na(paired) ~ 'paired TCR',
                                                    !is.na(chain_TRA)~'single TRA',
                                                    !is.na(chain_TRB)~ 'single TRB') )

Feature_rast(CD4CD8, 'TCR_summary', do.label = F)
table(CD4CD8$TCR_summary)

TCRmapping <- CD4CD8@meta.data %>% group_by(tissue, patient) %>%
  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)
 (ggplot(TCRmapping,aes(x = tissue,y= percent, group = TCR_summary, fill = TCR_summary))+geom_bar(stat = 'identity',position = position_stack(reverse = T))+facet_grid(~patient)+fill_m()+theme_bw()+mytheme) %T>% figsave('mappingrate_TCR.new.pdf', 150, 60)


# TCR clonal expansion ----------------------------------------------------

TCRabclonality_allcell <- Feature_rast(CD4CD8, 'cdr3_paired_perc', color_grd = 'grd', navalue = 'transparent')+ggtitle('TCRab clonality (%)')
patientID  <- pt %>% unique %>% sort

TCRabclonality_pt <-   map(patientID,~ CD4CD8 %>% subset(patient %in% .x) %>% 
                    Feature_rast('cdr3_paired_perc', color_grd = 'grd', sz = 0.5,navalue = 'transparent')+
                    ggtitle(paste0('TCRab clonality (%) in ', .x))) %>% 
  PG() 

TCRabclonality_fig <- PG(list(TCRabclonality_allcell, TCRabclonality_pt)) %T>% figsave('TCRabclonality_fig.pdf', 180, 100)

PG(list(TCRabclonality_allcell, TCRabclonality_pt))
saveRDS(CD4CD8, CD4CD8RDS)


CD4CD8$cdr3_paired_freq

CD4CD8@meta.data %<>% mutate(clonal_expansion =case_when(cdr3_paired_freq == 1 ~ 'monoclonal',
                                                       nr(cdr3_paired_freq, 2,9)~'low (2~9)',
                                                       nr(cdr3_paired_freq, 10,20)~'moderate (10~20)',
                                                       
                                                       cdr3_paired_freq >20 ~ 'high (>20)' ),
                             clonal_expansion = factor(clonal_expansion, levels = c('high (>20)', 'moderate (10~20)',
                                                                                    'low (2~9)', 'monoclonal'))   )

clone_exp_umap <- Feature_rast(CD4CD8, 'clonal_expansion',c('patient', 'tissue'),do.label = F, colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'), facetcol = 4) %T>% figsave('CD4CD8_clonalexpansion_tissue_patient.pdf',270,160) 

Feature_rast(CD4CD8, 'clonal_expansion',do.label = F, colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'))%T>% figsave('CD4CD8_TCR_clonalexpansion_UMAP.pdf',120,100) 

Feature_rast(CD4CD8 %>% subset(CD4CD8  %in% c('CD4', 'CD8')), 'clonal_expansion',c('tissue', 'CD4CD8'),do.label = F, colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2')) 
# top5 clones in CD4 and CD(   

top5 <- CD4CD8@meta.data %>% filter(CD4CD8 %in% c('CD4', 'CD8') & !is.na(paired)) %>% 
  group_by(CD4CD8, patient, cdr3_paired) %>% 
  summarise(pairedfreq = n()) %>% top_n(5, pairedfreq) %>% arrange(desc(pairedfreq))  %>%
  mutate(top5_paired_CD4 = case_when(CD4CD8 == 'CD4'~ cdr3_paired),
         top5_TCR_CD4 = case_when(CD4CD8 == 'CD4'~ paste0(patient, ' #', sprintf("%02d", 1:n()))),
         top5_TCR_CD8 = case_when(CD4CD8 == 'CD8'~ paste0(patient, ' #', sprintf("%02d", 1:n()))),
         
         top5_paired_CD8 = case_when(CD4CD8 == 'CD8'~ cdr3_paired) ) %>% ungroup() %>% 
  group_split(CD4CD8) %>% setNames(c('CD4', 'CD8'))

top5
  top10 <- CD4CD8@meta.data %>% filter(CD4CD8 %in% c('CD4', 'CD8') & !is.na(paired)) %>% 
  group_by(CD4CD8, patient, cdr3_paired) %>% 
  summarise(pairedfreq = n()) %>% filter(pairedfreq>3)  %>% top_n(10, pairedfreq) %>% arrange(desc(pairedfreq))  %>%
  mutate(top10_paired_CD4 = case_when(CD4CD8 == 'CD4'~ cdr3_paired),
         top10_TCR_CD4 = case_when(CD4CD8 == 'CD4'~ paste0(patient, ' #', sprintf("%02d", 1:n()) )),
         top10_TCR_CD8 = case_when(CD4CD8 == 'CD8'~ paste0(patient, ' #', sprintf("%02d", 1:n()) )),
         top10_paired_CD8 = case_when(CD4CD8 == 'CD8'~ cdr3_paired) ) %>% ungroup() %>% 
  group_split(CD4CD8) %>% setNames(c('CD4', 'CD8'))


top5$CD4

top10$CD4$top10_TCR_CD4 
CD4CD8@meta.data  %<>%  select(-contains('top'))

CD4CD8@meta.data  %<>%  left_join(select(top5$CD4, c(1:3, 5,6))) %>% left_join(select(top5$CD8, c(1:3, 7,8))) %>% `rownames<-`(CD4CD8$bc_backup)

CD4CD8@meta.data  %<>%  left_join(select(top10$CD4, c(1:3, 5,6))) %>% left_join(select(top10$CD8, c(1:3, 7,8))) %>% `rownames<-`(CD4CD8$bc_backup)
# CD4 top5 

top10CD4_umap <- map(patientID, ~
      Feature_rast(subset(CD4CD8, CD4CD8 == 'CD4' & patient == .x),  'top10_TCR_CD4', do.label = F, sz = 1)+
        ggtitle(.x )+ 
        guides(color = guide_legend(ncol =2, override.aes = list(size = 1.5))) ) %>% PG(ncol = 1)

top10CD8_umap <- map(patientID, ~
                      Feature_rast(subset(CD4CD8, CD4CD8 == 'CD8' & patient == .x),  'top10_TCR_CD8', do.label = F, sz  = 1)+
                       ggtitle(.x )+ 
                       guides(color = guide_legend(ncol =2, override.aes = list(size = 1.5))) ) %>% PG(ncol = 1)



top10expanded<- PG(list(top10CD4_umap, top10CD8_umap), ncol = 2, labels = c('Top 10 CD4 TCR', 'Top 10 CD8 TCR'),label_fontface = 'plain') %T>% figsave('top10expanded_CD4CD8.pdf', 200,250)

# Feature_rast(CD4CD8, 'cdr3_TRA_perc', color_grd = 'grd')



# TCR VDJ -----------------------------------------------------------------

CD4CD8@meta.data %<>% mutate(TRAVJ = case_when(!is.na(v_gene_TRA)  & !is.na(j_gene_TRA) ~ paste(v_gene_TRA, j_gene_TRA)  ) 
                            )

saveRDS(CD4CD8, CD4CD8RDS)
# mait cells 


# TCRJsharing Edge plot  --------------------------------------------------





# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)


# GDTlung_s@meta.data %<>%  mutate(pheno = str_replace(pheno, "memory", "circ") ) 

CD4TRB_data_meta <- CD4CD8@meta.data %>%  
  filter(cdr3_paired_freq > 1 &  CD4CD8 == "CD4" & 
           (Cell_pheno %in% c("TCM_1_M", "TCM_2_M", "TCM_3_M", "TEM_M", "Temra_1_LG", "Th17_M", "TRM_1_LG"))) %>% 
  # filter((grepl("LG", Cell_pheno))  & tissue == "Lung" | !(grepl("LG", Cell_pheno) ) & tissue == "LLN")  %>% 
  filter(!(str_detect(Cell_pheno, "TCM") & tissue == "LLN")) %>%
  mutate(Cell_pheno = str_replace(Cell_pheno, "TCM_\\d_M", "TCM_LN") )  

CD4TRB_data_meta$Cell_pheno

CD4_TRB_data <- CD4TRB_data_meta %>%  group_by(Cell_pheno) %>%  
  dplyr::count(cdr3_paired) %>% arrange(Cell_pheno, desc(n))%>% 
  # filter(n>1 )%>%
  mutate(Cell_pheno = factor(Cell_pheno, levels= c("TCM_LN", "Th17_M", "TRM_1_LG", "TEM_M", "Temra_1_LG")),
         uniquename = paste(Cell_pheno, cdr3_paired))

# GDTlung_s@meta.data %<>%  mutate(pheno = str_replace(pheno, "memory", "circ") ) 

CD8TRB_data_meta <- CD4CD8@meta.data %>%  
  filter(cdr3_paired_freq > 1 &  CD4CD8 == "CD8" & 
           (Cell_pheno %in% c("TCM_1_M", "TCM_2_M", "TCM_3_M", "TEM_M", "Temra_1_LG", "Th17_M", "TRM_1_LG", "TRM_3_LG"))) %>% 
  filter(!(str_detect(Cell_pheno, "TCM") & tissue == "LLN")) %>%
  mutate(Cell_pheno = str_replace(Cell_pheno, "TCM_\\d_M", "TCM_LN") )  


CD8TRB_data_meta$Cell_pheno

CD8_TRB_data <- CD8TRB_data_meta %>%  group_by(Cell_pheno) %>% 
  dplyr::count(cdr3_paired) %>% arrange(Cell_pheno, desc(n))%>% 
  # filter(n>1 )%>%
  mutate(Cell_pheno = factor(Cell_pheno, levels= c("TCM_LN", "Th17_M", "TRM_1_LG", "TEM_M", "Temra_1_LG", "TRM_3_LG")),
         uniquename = paste(Cell_pheno, cdr3_paired)) 

CD8_TRB_data
CD4_TCR_sharing <-  plot_tcr_sharing(CD4_TRB_data, title = "CD4+ TCR Sharing Across Clusters") %T>% print()


CD8_TCR_sharing <-  plot_tcr_sharing(CD8_TRB_data, title = "CD8+ TCR Sharing Across Clusters") %T>% print()





# TCR similarity ----------------------------------------------------------
# library(scRepertoire)
library(immunarch)




Total_list_CD4 <- CD4CD8@meta.data %>% mutate(aa = cdr3_TRB) %>% 
  filter(!is.na(cdr3_TRB) & CD4CD8 == 'CD4')   %>% 
  split(f = .$Cell_pheno)

Total_list_CD8 <- CD4CD8@meta.data %>% mutate(aa = cdr3_TRB) %>% 
  filter(!is.na(cdr3_TRB) & CD4CD8 == 'CD8')   %>% 
  split(f = .$Cell_pheno)



# 
# 



# library(conflicted)

conflicted::conflict_prefer_all('dplyr')

Total_list1_CD4 <- CD4CD8@meta.data %>% mutate(CDR3.nt = cdr3_nt_TRB) %>% filter(!is.na(CDR3.nt) & CD4CD8 == 'CD4') %>% 
  group_by(Cell_pheno) %>% 
  count(CDR3.nt, name = 'Clones') %>% 
  mutate(Proportion = Clones/sum(Clones)) %>%  as.data.frame %>% 
  
  split(f = .$Cell_pheno)



Total_list1_CD8 <- CD4CD8@meta.data %>% mutate(CDR3.nt = cdr3_nt_TRB) %>% 
  filter(!is.na(CDR3.nt) & CD4CD8 == 'CD8' & cdr3_TRB_freq > 1) %>% 
  group_by(Cell_pheno) %>% 
  count(CDR3.nt, name = 'Clones') %>% 
  mutate(Proportion = Clones/sum(Clones)) %>%  as.data.frame %>% 
  
  split(f = .$Cell_pheno)



Mori_result_CD4 <- repOverlap(Total_list1_CD4,.col = 'nt',
                              .method = "morisita", .verbose = F)

vis(Mori_result_CD4)+ 
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index CD4 TCRs')+xlab('Cluster')+ylab('Cluster')



Mori_result_CD8 <- repOverlap(Total_list1_CD8, .col='nt',
                              .method = "morisita", .verbose = F)

vis(Mori_result_CD8)+ 
  scale_fill_gradientn( na.value = 'white',
                         colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index CD8 TCRs')+xlab('Cluster')+ylab('Cluster')




  Mori_result_CD4 <- repOverlap(Total_list1_CD4, .col='nt',
                              .method = "morisita", .verbose = F)

vis(Mori_result_CD4)+ 
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index')+xlab('Cluster')+ylab('Cluster')


# Morisita per patient

Total_list1_CD8_pt <- CD4CD8@meta.data %>% mutate(CDR3.aa = cdr3_TRB) %>% filter(!is.na(CDR3.aa) & CD4CD8 == 'CD8') %>%    split(f = .$patient)  %>% map(
  ~  .x %>%  group_by(Cell_pheno) %>% 
    count(CDR3.aa, name = 'Clones') %>% 
    mutate(Proportion = Clones/sum(Clones)) %>%  as.data.frame %>% 
    
    split(f = .$Cell_pheno)
  
)


Total_list1_CD8_pt$p27

Mori_result_CD8_pt <- map(Total_list1_CD8_pt, ~ 
                            repOverlap(.x, .method = "morisita", .verbose = F) %>% 
                            vis()+ 
                            scale_fill_gradientn( na.value = 'white',
                                                  colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
                            ggtitle('Morisita index')+xlab('Cluster')+ylab('Cluster')     
                            
                            )

PG(Mori_result_CD8_pt)

Mori_result_CD8 <- repOverlap(Total_list1_CD8, .method = "morisita", .verbose = F)

vis(Mori_result_CD8)+ 
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index')+xlab('Cluster')+ylab('Cluster')







Morisita_GDTlung

clonalOverlap(Total_list1.1, cloneCall="cdr3_TRD", method="morisita")

clonalOverlap(Total_list[c(1,2,8)], cloneCall="cdr3_TRD", method="morisita")




Patient_list<- SplitObject(GDTlung_s, split.by = 'patient') %>%  map( ~ .x@meta.data %>% filter(!is.na(cdr3_TRD))    %>%  split(f = .$Cell_cluster) )

Morisita_patient <- map(Patient_list,~ 
                          clonalOverlap(.x, cloneCall="cdr3_TRD", method="morisita")  )


clonalOverlap(Patient_list$p27, cloneCall="cdr3_TRD", method="morisita") 


test <- subset(GDTlung_s, subset = is.na(cdr3_TRD) )




test<- clonalOverlap(combined3, cloneCall="cdr3_nt_TRD", method="morisita")


clonalOverlap(combined2, cloneCall="cdr3_TRD", method="overlap")


# gini index ---------------------------------------------------------

library(reldist)

ABTPairedfreq <-  CD4CD8@meta.data%>% 
  filter(!is.na(paired) & Cell_pheno != "unidentified_LG"  ) %>% 
  dplyr::count(CD4CD8, patient, Cell_pheno, cdr3_TRB) %>%
  arrange(desc(n))

ABTPairedfreq$Cell_pheno



gini_ABTpaire <-  ABTPairedfreq %>% 
  dplyr::group_by(CD4CD8, patient, Cell_pheno) %>%
  summarise(Gini_Index = gini(n), sum = sum(n))  %>%
  mutate(Gini_Index = replace(Gini_Index, sum < 10, NA))%>% 
  filter(Cell_pheno != "unidentified_LG"  ) %>% 
  dplyr::rename(CD4orCD8 = CD4CD8 ) %>% ungroup() %>% 
  tidyr::complete(CD4orCD8, patient, 
                  fill = list(Gini_Index = NA, sum = NA))





gini_ABTpaire$Cell_pheno



giniindex_ABT <- ggplot(gini_ABTpaire, aes(x = Cell_pheno , y = Gini_Index, 
                          color = patient, group= Cell_pheno))+geom_boxplot()+geom_point()+facet_wrap(~CD4orCD8,ncol = 1)+color_m(color = set_sample(umap.colors, s = 22))+theme_minimal()  +
  theme(axis.text.x = element_text(angle = 90))

figsave(giniindex_ABT, 'CD4CD8_giniindex.pdf', 100, 150)

giniindex_ABT

# MAIT cells by TCR -------------------------------------------------------




Feature_rast(subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20')) ,
            c( 'TRAVJ'), noaxis = F, sz = 1, facets = 'patient')

Feature_rast(subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20')) ,
             c( 'TRAVJ', 'v_gene_TRB',  'cdr3_TRA_perc'), noaxis = F, sz = 1)

CD4CD8$tissue


BG = Feature_rast(CD4CD8, 'bg', do.label = F)




BG +(geom_point_rast(data =    subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20' ))  %>% FetchData(c('UMAP_1', 'UMAP_2', 'TRAVJ', 'c_gene_TRB', 'patient', 'tissue')
), aes(x = UMAP_1, y = UMAP_2, color = TRAVJ, size  = 1) +scale_color_manual(values = ggplotColours(12)) )                  
                     )+ggtitle('TRAV1-2 J12 TRBV6')+ facet_grid('patient')


DimPlot(CD4CD8)


CD4CD8$v_gene_TRA %>%  unique()

Feature_rast(CD4CD8, 'patient')

Feature_rast(CD4CD8, c('tissue', 'KLRB1', 'DPP4', 'RORC' , 'CD4', 'CD8A'), color_grd = 'grd', ncol =3)



CD4CD8@meta.data %<>%  mutate( unconventional = case_when(
  v_gene_TRA == 'TRAV1-2' & 
    j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20') & 
       v_gene_TRB %in% c('TRBV20OR9-2', 'TRBV6-1', 'TRBV20-1')     ~  'MAIT',
  v_gene_TRA == 'TRAV10' & j_gene_TRA == 'TRAJ18'  & v_gene_TRB == 'TRBV25-1'  ~  'iNKT'   ) ) 

Feature_rast(CD4CD8, 'unconventional', navalue = alpha('lightgrey', 0.05), do.label = F, facets = 'tissue')


ViolinPlot(CD4CD8, 'CD103.protein', group.by = 'unconventional', assay = 'CITE')

CD4CD8@meta.data %>%  count(unconventional, patient, tissue)

Feature_rast(CD4CD8, 'TRAV1-2', assay = 'RNA')
FeaturePlot(CD4CD8, 'TRAV1-2',)

grep('TRAV1', rownames(CD4CD8@assays$RNA@counts), value = T)
rownames(CD4CD8)

CD4CD8$cdr3_paired_perc

Feature_rast(CD4CD8,'cdr3_TRB_perc', facets = 'CD4CD8', navalue = 'transparent')




#TCR most expanded clones in each cluster and each patients ------------------





totalTRBcounts<- CD4CD8@meta.data %>% group_by(patient, CD4CD8, Cell_pheno) %>%  
  filter(!is.na(cdr3_nt_TRB)) %>%  summarise(totalTRB = n() ) 
totalTRBcounts

top1TRBperclperpt <- CD4CD8@meta.data %>% group_by(patient, CD4CD8, Cell_pheno, cdr3_nt_TRB) %>%  
  filter(!is.na(cdr3_nt_TRB)) %>%  summarise(TRBcount = n() ) %>% filter(TRBcount> 2) %>% 
  group_by(patient, CD4CD8, Cell_pheno)%>% top_n(1, TRBcount) %>% 
  left_join(totalTRBcounts)  %>%  mutate(MostExp = paste0(
    'Most_expanded_',CD4CD8,'_', Cell_pheno, '_clone'
    
  )) 

top1TRBperclperpt_s <- spread(top1TRBperclperpt, key = 'MostExp', value = TRBcount )  %>% 
  mutate(across(starts_with("Most_expanded_"), ~ !is.na(.))) %>% 
  group_by(patient, CD4CD8, cdr3_nt_TRB) %>% 
  summarize(across(starts_with("Most_expanded_"), any, na.rm = TRUE))


rowSums(top1TRBperclperpt_s[,4:24])

spread(top1TRBperclperpt, key = 'MostExp', value = TRBcount )  

CD4CD8@meta.data  %<>% select_at(vars(-contains('Most_expanded')))
  
CD4CD8@meta.data %<>% 
  left_join( top1TRBperclperpt_s,
                                 by = c('cdr3_nt_TRB', 'CD4CD8', 'patient') ,
             suffix = c('', '')) %>%
  mutate(across(starts_with("Most_expanded_"), ~ ifelse(., ., NA), .ptype = logical())) %>% 

  `rownames<-`(CD4CD8$bc_backup)

Feature_rast(CD4CD8, 'Most_expanded_CD8_TMRA_1_P_clone',
              do.label = F)



CD4CD8$bg <- NA

CD4CD8_meta <- CD4CD8@meta.data


Feature_rast(CD4CD8)


CD4CD8_meta  <- cbind(FetchData(CD4CD8, c('UMAP_1', 'UMAP_2')), CD4CD8_meta)

mostcolnames <- grep('Most_expanded_', colnames(CD4CD8_meta), value = T)

allsubdata <- map(mostcolnames, ~CD4CD8_meta[!is.na(CD4CD8_meta[[.]]), ])%>%  setNames(mostcolnames)


phenocolors  <- setNames(umap.colors[1: length(levels(CD4CD8_meta$Cell_pheno))], 
                         levels(CD4CD8_meta$Cell_pheno))


CD4CD8_bg <- Feature_rast(CD4CD8, 'bg', sz = 0.1) + ggtitle(NULL)
CD4CD8_meta$patient
ALLMOSTEXPFIGS <- 
  map(mostcolnames, ~ 
        CD4CD8_bg+geom_point_rast(data = CD4CD8_meta[!is.na(CD4CD8_meta[[.]]), ],
                                  aes(x = UMAP_1, y = UMAP_2, color = Cell_pheno) , size = 0.5
        ) + scale_color_manual(values = phenocolors, na.value = alpha('lightgrey', 0.4)) +ggtitle(.)  + NoLegend  ()  ) %>%  setNames(mostcolnames)

ALLMOSTEXPFIGS <- 
  map(mostcolnames, ~ 
        CD4CD8_bg+geom_point_rast(data = CD4CD8_meta[!is.na(CD4CD8_meta[[.]]), ],
                                  aes(x = UMAP_1, y = UMAP_2, color = patient) , size = 0.5
        ) + color_m() +ggtitle(.)  + NoLegend  ()  ) %>%  setNames(mostcolnames)




ALLMOSTEXPFIGS$Most_expanded_CD4_TCM_2_M_clone

lg <- cowplot::get_legend(Feature_rast(CD4CD8, "patient"))


 PG(ALLMOSTEXPFIGS[c(3, 4,7, 13, 16, 18, 19, 20)], ncol = 4)
 
 
 PG(ALLMOSTEXPFIGS[c(3, 4,7, 13, 16, 18, 19, 20)], ncol = 4)
 

 EXP_memory_clones <-   PG(ALLMOSTEXPFIGS[c(3, 4,7, 13, 16, 18, 19, 20)], ncol = 4) %>% 
   list(lg) %>% 
   PG( rw = c(5.5,1)) %T>% figsave('CD4CD8_Expanded_memeory_T_clones.pdf', 210, 80) 
 
   
 PG(ALLMOSTEXPFIGS[c(13, 16, 18,  20)], ncol = 4) %>% 
   list(lg) %>% 
   PG( rw = c(6,1))    
   
 ALLMOSTEXPFIGS$Most_expanded_CD8_TMRA_1_P_clone +facet_wrap(~patient)



ALLMOSTEXPFIGS_tissue <- 
  map(mostcolnames, ~ 
        CD4CD8_bg+geom_point_rast(data = CD4CD8_meta[!is.na(CD4CD8_meta[[.]]), ],
                                  aes(x = UMAP_1, y = UMAP_2, color = tissue)
        ) + 
color_m()+
  ggtitle(.)   ) %>%  setNames(mostcolnames)
ALLMOSTEXPFIGS_tissue$Most_expanded_CD4_TRM_1_P_CDR3_TRB
PG(ALLMOSTEXPFIGS_tissue[c(3, 4,7, 13, 16, 18, 19, 20)], ncol = 4)



ALLMOSTEXPFIGS_tissue$Most_expanded_CD8_TMRA_1_P_clone+ facet_wrap(~patient) 

ALLMOSTEXPFIGS_tissue$Most_expanded_CD8_TMRA_1_P_CDR3_TRB + facet_wrap(~patient) 


ALLMOSTEXPFIGS_tissue$Most_expanded_CD8_TRM_3_P_CDR3_TRB + facet_wrap(~patient) 

# TCRsharing --------------------------------------------------------------


phenocolors

# TCR sharing between CD4 CD8 &  Pulm and LN

C14TCR <- CD4CD8@meta.data %>% 
  filter(!is.na(cdr3_paired) & 
           CD4CD8 %in% c('CD4', 'CD8')) %>% 
  select(cdr3_paired,CD4CD8) %>% group_split(CD4CD8)


intersect(C14TCR[[1]]$cdr3_paired,C14TCR[[2]]$cdr3_paired)



tissueTC <- CD4CD8@meta.data %>% filter(!is.na(cdr3_paired) &CD4CD8 %in% c('CD8')) %>% select(cdr3_paired,tissue) %>% group_split(tissue)


TCRpatient <- CD4CD8@meta.data %>% filter(!is.na(cdr3_paired) ) %>% select(cdr3_paired,patient) %>% group_split(patient) %>%set_names(patientID) %>% map(~ .x %>% pull(cdr3_paired) %>% unique)
  
(TCRpatient %>% unlist() %>% table() %>% sort(decreasing = T) )[1:100]

tissueTC[[2]]$cdr3_paired
 

intersect(tissueTC[[1]]$cdr3_paired,tissueTC[[2]]$cdr3_paired)

TCRbyCD4CD8 <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1  & CD4CD8 %in% c('CD4', 'CD8')) %>% select(cdr3_paired,CD4CD8) %>%group_by(CD4CD8,cdr3_paired) %>%  summarise(pairedfreq = n()) %>% arrange(CD4CD8) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )
  arrange(pairedfreq)
# nrow(TCRbyCD4CD8 %>% filter())
  
  
  



TCRBbyCD4CD8 <-CD4CD8@meta.data %>% filter(cdr3_TRB_freq >1  & CD4CD8 %in% c('CD4', 'CD8')) %>% select(cdr3_TRB,CD4CD8) %>%group_by(CD4CD8,cdr3_TRB) %>%  summarise(TRBfreq = n()) %>% arrange(CD4CD8) %>% mutate(cdr3_TRB = factor(cdr3_TRB, unique(cdr3_TRB)) )




CD4CD8$cdr3_nt_TRB

TCRB_shared_by_CD4CD8 <- CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRB)) %>% 
  group_by(CD4CD8, patient) %>%  distinct(cdr3_TRB)  %>%  ungroup() %>%  group_by(cdr3_TRB,patient) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) %>%  pull(cdr3_TRB)


TCRB_shared_by_CD4CD8 <- CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRB)) %>% 
  group_by(CD4CD8, patient) %>%  distinct(cdr3_nt_TRB, v_gene_TRB)  %>%  ungroup() %>%  group_by(v_gene_TRB,cdr3_nt_TRB,patient) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) %>%  pull(cdr3_nt_TRB)


CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRB)) %>% 
  group_by(CD4CD8, patient) %>%  distinct(cdr3_nt_TRB, v_gene_TRB)  %>%  ungroup() %>%  group_by(v_gene_TRB,cdr3_nt_TRB,patient) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) %>%  pull(v_gene_TRB)%>% table()

TCRB_shared_by_CD4CD8$TCRn %>%  unique()

TCRA_shared_by_CD4CD8 <- CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRA)) %>% 
  group_by(CD4CD8) %>%  distinct(cdr3_TRA)  %>%  ungroup() %>% 
  group_by(cdr3_TRA) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) %>%  pull(cdr3_TRA) 



TCRA_shared_by_CD4CD8 <- CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRA)) %>% 
  group_by(CD4CD8, patient) %>%  distinct(cdr3_nt_TRA, v_gene_TRB)  %>%  ungroup() %>%  group_by(v_gene_TRB,cdr3_nt_TRA,patient) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) %>%  pull(cdr3_nt_TRA)


TCRB_shared_by_patient <- CD4CD8@meta.data %>%  filter(!is.na(cdr3_TRB)) %>% 
  group_by(patient) %>%  distinct(cdr3_TRB)  %>%  ungroup() %>%  group_by(cdr3_TRB,patient) %>% 
  summarise(TCRn = n()) %>%  filter(TCRn >1) 

TCRB_shared_by_patient


CD4CD8@meta.data %<>% 
  mutate(TCR_shared_CD4CD8  = 
           case_when(
             cdr3_nt_TRB %in% TCRB_shared_by_CD4CD8 & cdr3_TRB_freq >= 10 ~  'expanded shared_TRB',
             cdr3_nt_TRB %in% TCRB_shared_by_CD4CD8 ~  'shared_TRB') ,
        TCRA_shared_CD4CD8  = 
           case_when(
            
             cdr3_nt_TRA %in% TCRA_shared_by_CD4CD8 ~  'shared_TRA') 
         
         
         )
CD4CD8$TCR_shared_CD4CD8 %>%  unique()

Feature_rast(CD4CD8, 'TCR_shared_CD4CD8', facets =  c('patient', 'CD4CD8'), colorset = 'gg', do.label = F, facetcol = 4)
Feature_rast(CD4CD8, 'TCRA_shared_CD4CD8', facets =  c('patient', 'CD4CD8'), colorset = 'gg', do.label = F, facetcol = 4)


saveRDS(CD4CD8, CD4CD8RDS)
CD4CD8@reductions$umap@cell.embeddings[,1]
CD4CD8$umap1_CD4CD8 <- CD4CD8@reductions$umap@cell.embeddings[,1]

Feature_rast(CD4CD8, noaxis = F, axis.number = T)

CD4CD8@meta.data  %<>%  mutate(umap1_CD4CD8 =
                                 case_when(CD4CD8 == 'CD8' ~umap1_CD4CD8 +8,
                                           CD4CD8 == 'CD4' ~umap1_CD4CD8 -8 ))


Feature_rast(CD4CD8, 'TCR_shared_CD4CD8',d1 = 'umap1_CD4CD8',colorset = 'gg',do.label = F)+NoLegend() +geom_line(aes(group =TCR_shared_CD4CD8))




# Paired TCR sharing between CD4 and CD8 ----------------------------------




library(ggalluvial)


TCRsharingCD4CD8 <- (TCRbyCD4CD8 %>% 
  ggplot(
    aes( x = CD4CD8 , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    ggtitle("TCR sharing between CD4 and CD8")+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
    theme_minimal_hgrid()+
    # scale_fill_manual(values = rainbow(884))+
    geom_stratum(size = 0.01, stroke = 0.1,  color = alpha('black', 0.2))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)
TCRsharingCD4CD8 %T>%  figsave('TCRsharing_CD4_CD8.pdf', 100, 100)


TCRbyCluster <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1  ) %>% select(cdr3_paired,Cell_pheno) %>%group_by(Cell_pheno,cdr3_paired) %>%  summarise(pairedfreq = n()) %>% arrange(Cell_pheno) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )


Lin_TRM3 <- c("TEM_M")

(TCRBbyCD4CD8 %>% 
    ggplot(
      aes( x = CD4CD8 , y = TRBfreq, fill = cdr3_TRB,  stratum= cdr3_TRB, alluvium  = cdr3_TRB))+
    ggtitle("TCR sharing between CD4 and CD8")+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
    theme_minimal_hgrid()+
    # scale_fill_manual(values = rainbow(1484))+
    geom_stratum(size = 0.1, color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRB frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)


# TCR shared by tissue ----------------------------------------------------



TCRbytissue <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >2 &CD4CD8 %in% c('CD8', 'CD4')) %>% select(cdr3_paired,tissue,CD4CD8) %>%group_by(tissue,cdr3_paired,CD4CD8) %>%  summarise(pairedfreq = n()) %>% arrange(CD4CD8, tissue) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )

nrow(TCRbytissue %>% filter())

TCRsharingtissue <- (TCRbytissue %>% 
                       ggplot(
                         aes( x = tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
                       ggtitle("TCR sharing between Lung and LLN")+
                       geom_flow(stat = "alluvium",
                                 color = "darkgray") +
                       # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
                       theme_minimal_hgrid()+
                       # scale_fill_manual(values = rainbow(1014))+
                       facet_wrap(~CD4CD8)+
                       geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
                       xlab(NULL) +ylab("TCRab frequencies")+
                       theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
                       guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)


TCRsharingtissue

TCRsharingplot <- PG(list(TCRsharingCD4CD8, TCRsharingtissue), rw = c(1,1.8)) %T>% figsave('TCRsharingplot_CD4_8_tissue.pdf', 150,120)

  

Feature_rast(CD4CD8)
# how TCRs shared between lung and LN?
CD4CD8_meta <- CD4CD8@meta.data
levels(CD4CD8_meta$Cell_pheno)
CD4CD8_meta$cluster_no


CD4CD8_meta$pheno_tissue <- paste0(CD4CD8_meta$tissue, '_', CD4CD8_meta$Cell_cluster)
CD4CD8_meta  %<>%  mutate(tissue = factor(tissue, c('Lung', 'LLN')))
TCR_paired_CD4CD8 <- CD4CD8_meta %>% filter(cdr3_paired_freq >1 ) %>% 
    mutate(pheno_tissue = case_when(tissue == 'LLN' ~ 'LLN',
                                    tissue == 'Lung' ~ pheno_tissue)) %>%
  group_split(CD4CD8)  %>%  setNames(c('CD4', 'CD8'))
  
TCR_paired_CD4CD8$CD8$Cell_cluster
p_phenos <- c('TRM_1_Lu', 'TRM_2_Lu', 'TRM_3_Lu', 'Temra_1_Lu', 'TEM_M', 'TCM_3_M')

TCR_paired_CD4CD8$CD8$pheno_tissue %>% unique()

CD8_sharedTCR_lung <- map(p_phenos, ~  

TCR_paired_CD4CD8$CD8 %>%  
  filter(tissue == 'LN' | grepl(.x,pheno_tissue )   ) %>% 
  group_by(pheno_tissue,cdr3_paired)  %>%  summarise(pairedfreq = n()) %>% ungroup() %>%
  arrange( pairedfreq)  %>% 
  mutate(cdr3_paired = factor(cdr3_paired, levels = unique(cdr3_paired))) %>%  
  ggplot(   aes( x = pheno_tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  theme_minimal_hgrid()+
  ggtitle(paste0('TCR clone sharing between\nCD8+ Pulm ',.x,' and CD8+ LN T cells'))+
  geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
  xlab(NULL) +
  ylab("TCRab frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"), 
        axis.text.x = element_text(angle = 60))+
  guides(fill = guide_legend(ncol = 1, title = NULL)) 
) %>% setNames(p_phenos)

PG(CD8_sharedTCR_lung)


CD4_sharedTCR_lung <- map(p_phenos, ~  
                            
                            TCR_paired_CD4CD8$CD4 %>%  
                            filter(tissue == 'LN' | grepl(.x,pheno_tissue )   ) %>% 
                            group_by(pheno_tissue,cdr3_paired)  %>%  summarise(pairedfreq = n()) %>% ungroup() %>%
                            arrange( cdr3_paired)  %>% 
                            ggplot(   aes( x = pheno_tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
                            geom_flow(stat = "alluvium",
                                      color = "darkgray") +
                            theme_minimal_hgrid()+
                            ggtitle(paste0('TCR clone sharing between\nCD4+ Pulm ',.x,' and CD4+ LN T cells'))+
                            geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
                            xlab(NULL) +
                            ylab("TCRab frequencies")+
                            theme(legend.position = 'none', legend.key.size = unit(2, "mm"), 
                                  axis.text.x = element_text(angle = 60))+
                            guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme
) %>% setNames(p_phenos)

PG(CD4_sharedTCR_lung)


# TCRbytissue_cluster <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 &CD4CD8 %in% c('CD8', 'CD4')) %>%
#   select(cdr3_paired,tissue,CD4CD8, Cell_pheno) %>%group_by(tissue,cdr3_paired,CD4CD8,Cell_pheno) %>%  summarise(pairedfreq = n()) %>% ungroup() %>% arrange(CD4CD8, tissue) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )
# TCRbytissue_cluster
# 
# Feature_rast(CD4CD8, 'cluster_no')
# 
# clusters <- c('c1', 'c2', 'c4',  'c6', 'c7', 'c12', 'c14')
# 
# clusters2 <- c('C1', 'C2', 'C4',  'C6', 'C7', 'C12', 'C14', 'C15', 'C17', 'C18')


TCRsharing_withincluster_tissue<-(TCRbytissue_cluster %>% filter(Cell_cluster %in% clusters) %>% 
    ggplot(
      aes( x = tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  geom_flow(stat = "alluvium",
              color = "darkgray") +
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(1100))+facet_wrap(~Cell_cluster+CD4CD8, scales = "free", ncol = 10)+
    geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 45))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) %T>% 
  figsave('TCRsharing_withincluster_tissue.pdf', 200,100)

TCRsharing_withincluster_tissue

(TCRbytissue_cluster %>% filter(Cell_cluster %in% clusters & CD4CD8 == 'CD8') %>% 
    ggplot(
      aes( x = tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(1100))+facet_wrap(~Cell_cluster, scales = "free", ncol = 5)+
    geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 45))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 



# TCR sharing between clusters
CD4CD8@meta.data  %<>% mutate(Cell_group = case_when(
  Cell_cluster %in% c('C8', 'C9', 'C10') ~'group 1: naive',
  Cell_cluster %in% c('C3', 'C6', 'C11', 'C5') ~'group 2: CD4 helpers',
  Cell_cluster %in% c('C1', 'C2', 'C4') ~'group 3: Tcm and Th1',
  Cell_cluster %in% c('C7', 'C12') ~'group 4: Tem',
  Cell_cluster %in% c('C15', 'C16', 'C17','C18', 'C19') ~'group 5: Trm effector',
  Cell_cluster %in% c('C14') ~'group 6: Th1_17 Trm'
)) %>% `rownames<-`(CD4CD8$bc_backup)

Feature_rast(CD4CD8, 'Cell_group')


TCRby_cluster <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 ) %>%
  # select(cdr3_paired,Cell_pheno) %>%
  group_by(cdr3_paired,Cell_pheno,) %>%  summarise(pairedfreq = n()) %>% ungroup() %>%
  # arrange(CD4CD8) %>%
  mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )

# CD4groups <- c('group 2: CD4 helpers', 'group 3: Tcm and Th1','group 4: Tem','group 6: Th1_17 Trm')

TRM3group <- c("TRM_1_P", "TRM_3_P", "Temra_1_P", "TEM_M",  "TCM_3_M")

TRM1group <- c("TRM_3_P", "TRM_1_P", "Th17_M",  "TCM_2_M")


TCRsharing_CD8_TRM3 <-(TCRby_cluster %>% 
    filter(Cell_pheno %in% TRM3group  ) %>% mutate(Cell_pheno = factor(Cell_pheno, TRM3group)) %>% 
    ggplot(
      aes( x = Cell_pheno , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    theme_minimal_hgrid()+
    # scale_fill_manual(values = rainbow(1100)[280:810])+
    # facet_wrap(~patient, scales = "free", ncol = 5)+
    geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+ggtitle('CD8 TCRab')+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 
TCRsharing_CD8_TRM3



# Top TCRs and gene expression  -------------------------------------------

# top3 CD8 in each donor in Lung 

top3_TRB_lung <-  CD4CD8@meta.data %>%  filter(CD4CD8 == "CD8" & 
                               !is.na(v_gene_TRB)&
                               tissue == "Lung") %>% 
  group_by(patient) %>% count(cdr3_TRB) %>% top_n(3, n) %>% 
  pull(cdr3_TRB)
  
top3_TRB_lung


Feature_rast(CD4CD8 %>% subset(cdr3_TRB %in% top3_TRB_lung), 
             g= "tissue", d1 = "Tissue.resident", d2 = "KLRG1", facets = "patient", noaxis = F, facetcol = 4, do.label = F) 



CD8_tissue <- (CD4CD8@meta.data %>% filter(cdr3_TRB %in% top3_TRB_lung) %>% 
  ggplot(aes(x = tissue, y = TRM_score, color = tissue))+
  geom_jitter_rast(size = 0.5)+
    geom_boxplot(color = "black", fill = "transparent", size = 0.3)+facet_wrap(~patient)+color_m()) %T>% print() 

CD8_tissue

Feature_rast(CD4CD8 %>% subset(cdr3_TRB %in% top3_TRB_lung), 
             g= "Tissue.resident",  facets = c("patient", "tissue"), noaxis = F, facetcol = 4, do.label = F) 


top3_TRB_lung_CD4 <-  CD4CD8@meta.data %>%  filter(CD4CD8 == "CD4" & 
                                                 !is.na(v_gene_TRB)&
                                                 tissue == "Lung") %>% 
  group_by(patient) %>% count(cdr3_TRB) %>% top_n(3, n) %>% 
  pull(cdr3_TRB)


Feature_rast(CD4CD8 %>% subset(cdr3_TRB %in% top3_TRB_lung_CD4), 
             g= "tissue", d1 = "Tissue.resident", d2 = "KLRG1", facets = "patient", noaxis = F, facetcol = 4, do.label = F) 



CD4_tissue <- (CD4CD8@meta.data %>% filter(cdr3_TRB %in% top3_TRB_lung_CD4) %>% 
                 ggplot(aes(x = tissue, y = Tissue.resident, color = tissue))+
                 geom_jitter_rast(size = 0.5)+
                 geom_boxplot(color = "black", fill = "transparent", size = 0.3)+facet_wrap(~patient)+color_m()) %T>% print() 

CD4_tissue


top3_TRD_lung <-  GDTlung_s@meta.data %>%  filter(v_gene_TRD != "TRDV2",
                                                     !is.na(v_gene_TRD)&
                                                     tissue == "Lung") %>% 
  group_by(patient) %>% count(cdr3_TRD) %>% top_n(3, n) %>% 
  pull(cdr3_TRD)



top3_TRDV2_lung <-  GDTlung_s@meta.data %>%  filter(v_gene_TRD == "TRDV2",
                                                  !is.na(v_gene_TRD)&
                                                    tissue == "Lung") %>% 
  group_by(patient) %>% count(cdr3_TRD) %>% top_n(3, n) %>% 
  pull(cdr3_TRD)


TRD_tissue <- (GDTlung_s@meta.data %>% filter(cdr3_TRD %in% top3_TRD_lung) %>% 
                 ggplot(aes(x = tissue, y = Tissue.resident, color = tissue))+
                 geom_jitter_rast(size = 0.5)+
                 geom_boxplot(color = "black", fill = "transparent", size = 0.3)+facet_wrap(~patient)+color_m()) %T>% print() 

TRD_tissue



TRDV2_tissue <- (GDTlung_s@meta.data %>% filter(cdr3_TRD %in% top3_TRDV2_lung) %>% 
                 ggplot(aes(x = tissue, y = Tissue.resident, color = tissue))+
                 geom_jitter_rast(size = 0.5)+
                   geom_boxplot(color = "black", fill = "transparent", size = 0.3)+facet_wrap(~patient)+color_m()) %T>% print() 

Feature_rast(GDTlung_s %>% subset(cdr3_TRD %in% top3_TRD_lung), 
             g= "tissue", d1 = "Tissue.resident", d2 = "KLRG1", facets = "patient", noaxis = F, facetcol = 4, do.label = F) 



top3intissue <- PG(list(CD8_tissue, CD4_tissue,
        TRD_tissue, TRDV2_tissue), ncol = 2, labels = c("CD8", "CD4", "nonVD2", "VD2")) 
figsave(top3intissue, "top3_TCR_of_lung.pdf", path = figpath_ni, 300, 300)



Feature_rast(GDTlung_s %>% subset(cdr3_TRD %in% top3_TRD_lung),             g= "Tissue.resident",  facets = c("patient", "tissue"), noaxis = F, facetcol = 4, do.label = F) 


expanded_TRM_CD8clone <-  CD4CD8@meta.data %>%  filter(CD4CD8 == "CD8" & 
                                                         Cell_pheno %in% c("TRM_1_LG", "TRM_2_LG", "TRM_3_LG"),
                                                         !is.na(v_gene_TRB)&
                                                         tissue == "Lung") %>% 
  group_by(patient) %>% count(cdr3_TRB) %>% filter(n >= 5) %>% 
  pull(cdr3_TRB)

expanded_TRM_CD8clone


CD8_tissue_TRM <- (CD4CD8@meta.data %>% filter(cdr3_TRB %in% expanded_TRM_CD8clone) %>% 
                 ggplot(aes(x = tissue, y = Tissue.resident, color = tissue))+
                 geom_jitter_rast(size = 0.5)+
                 geom_boxplot(color = "black", fill = "transparent", size = 0.3)+facet_wrap(~patient)+color_m()) %T>% print() 


# public data satija ------------------------------------------------------

Pulm_Satija <- readRDS('/home/big/tanlikai/lung/public/Azimuth.Pulm.Satija.rds')



CD8diffcopd <- ClusterCompare(Pulm_Satija %>% subset(annotation.l1 == 'CD8 T'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

Mfdiffcopd <- ClusterCompare(Pulm_Satija %>% subset(annotation.l1 == 'Macrophage'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

Basaldiffcopd <- ClusterCompare(Pulm_Satija %>% subset(annotation.l1 == 'Basal'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

CD14monodiffcopd <-
  ClusterCompare(Pulm_Satija %>% subset(annotation.l1 == 'CD14+ Monocyte'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

batcheff_copd_normal <- intersect(CD8diffcopd$table$gene, Basaldiffcopd$table$gene) %>% intersect(CD14monodiffcopd$table$gene) 
intersect(CD8diffcopd$table$gene, CD14monodiffcopd$table$gene)



intersect(CD8diffcopd$table$gene, Basaldiffcopd$table$gene) %>% intersect(CD14monodiffcopd$table$gene)  %>% sort()

Feature_rast(Pulm_Satija, c('tissue','TRDC', 'TRGC1', 'TRGC2', 'CD3D', 'CD3E') ,colorset = 'gg')
Feature_rast(Pulm_Satija, c('TRAC','TRDC', 'TRGC1', 'TRGC2', 'CD3D', 'CD3E') ,sz = 0.2)

Pulm_Satija$tissue %>% table()
tcells <-  Pulm_Satija$cell_type %>% unique() %>% str_subset('T cell|thy')

diseases  <-  c('chronic obstructive pulmonary disease', 'normal', 'COVID-19')

PulmT_Satija <-  readRDS('/home/big/tanlikai/lung/public/Azimuth.Pulm.Satija.rds') %>% 
  subset(cell_type %in% c(.$cell_type %>% unique() %>% str_subset('T cell|thy')) & disease %in% diseases)
dim(PulmT_Satija)

colnames(PulmT_Satija@meta.data)


Feature_rast(PulmT_Satija, c( 'CD3D', 'CD3E'), colorset = 'gg')
Feature_rast(PulmT_Satija, c( 'disease', 'health_status', 'tissue'), colorset = 'gg', sz = 0.1)
PulmT_Satija$disease %>% table
table(PulmT_Satija$disease, PulmT_Satija$health_status)


Feature_rast(PulmT_Satija, c('CD4', 'CD8A', 'RORC', 'CCR6', 'FOXP3', 'TOX', 'ITGAE' , 'ITGA1'))



load('/home/big/tanlikai/lung/public/GSE162498_NSCLC_CD3_4tumors_4Juxta_2Juxta.Rds')
Feature_rast(eleven.tils.cd3.integrated, 'IL17A')
Feature_rast(eleven.tils.cd3.integrated)



three.tissues.cd3.integrated$tissue %>% unique()
dim(three.tissues.cd3.integrated)

Feature_rast(three.tissues.cd3.integrated, c('CD4', 'CD8A', 'CD8B', 'ITGAE', 'ITGA1'))


three.tissues.cd3.integrated@assays$RNA@counts %>% dim()


table(three.tissues.cd3.integrated$tissue)





# integration with satija T celss -----------------------------------------

# clean the data  

PulmT_Satija  %<>% 
PercentageFeatureSet( '^MT', col.name =  'percent.mito') %>% 
  PercentageFeatureSet('^RP', col.name = 'percent.ribo')    %>% 
  PercentageFeatureSet('^HSPA', col.name =  'percent.hspa')



ViolinPlot(PulmT_Satija, 'nCount_RNA', box = T) +ylim(0, 10000)
ViolinPlot(PulmT_Satija, 'nFeature_RNA', box = T) +ylim(0, 5000)
ViolinPlot(PulmT_Satija, 'percent.mito', box = T) 

Feature_rast(PulmT_Satija, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', color_grd = 'grd',
             
             noaxis = F, axis.number = T)+
  geom_smooth(method = "lm")+
  
  scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
  scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
  geom_hline(yintercept = c(500,2500))+geom_vline(xintercept = c(1000,6500))



PulmT_Satija  %<>% subset( nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA > 1000 & nCount_RNA < 6500& percent.mito < 16)

dim(PulmT_Satija)
PulmT_Satija$donor %>%  unique()
table(PulmT_Satija$donor) %>% sort(decreasing = T)

table(PulmT_Satija$donor, PulmT_Satija$disease)


TCELLgenes <- c('PTPRC','CD3D', 'CD3G', 'CD3E', 'TRAC', 'TRBC1', 'TRBC2','TRDC', 'TRGC1', 'TRGC2', 'CD4', 'CD8A', 'CD8B')

Feature_rast(PulmT_Satija, c('donor', 'disease', 'tissue'), colorset = 'gg', ncol =1)
Feature_rast(PulmT_Satija, TCELLgenes, noaxis = F)
ViolinPlot(PulmT_Satija, c('CD3D', 'CD3G', 'CD3E', 'TRAC', 'TRBC1', 'TRBC2','TRDC', 'TRGC1', 'TRGC2'),colors = umap.colors)


# select real T cells 
PulmT_Satija  %<>%  AddModuleScore(features = list(c('CD3D', 'CD3G', 'CD3E')),name = 'CD3score')

ViolinPlot(PulmT_Satija, c('CD3score1'),colors = umap.colors)+geom_hline(yintercept = c(-0.1,0.6))

Feature_rast(PulmT_Satija, 'CD3score1', color_grd = 'grd', sz = 0.1)

PulmT_Satija$bc_backup <- rownames(PulmT_Satija@meta.data)  

PulmT_Satija@meta.data  %<>% mutate(CD3exp = case_when(CD3score1 < 0 ~ 'CD3neg',
                                                       nr(CD3score1, 0, 0.6) ~ 'CD3lo',
                                                       CD3score1 > 0.6 ~ 'CD3hi')) %>% `rownames<-`(PulmT_Satija$bc_backup)


Feature_rast(PulmT_Satija, 'CD3exp',  sz = 0.5, colorset = 'gg', noaxis = F, axis.number = T) +
  geom_hline(yintercept = c(0,9))+
  geom_vline(xintercept = c(-11,0))

table(PulmT_Satija$CD3exp)



PulmT_Satija  %<>%  subset(CD3exp %in% c('CD3lo', 'CD3hi') & UMAP_1 > -11 & UMAP_1 < 0 & UMAP_2 > 0 & UMAP_2 < 9 )

dim(PulmT_Satija)

Feature_rast(PulmT_Satija, c('donor', 'tissue', 'disease'), colorset = 'gg', ncol =1)


PulmT_Satija  %<>%  FindVariableFeatures(selection.method = 'vst', nfeatures = 5000)

PulmT_Satija@assays$RNA@data


vf1 <- PulmT_Satija@assays$RNA@var.features
vf2 <- PulmT_Satija@assays$RNA@var.features
intersect(vf1, vf2)

PulmT_Satija@assays$RNA@var.features %<>%str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T) 
VariableFeaturePlot(PulmT_Satija) +ylim(-1, 10)
VariableFeaturePlot(CD4CD8)

PulmT_Satija@assays$RNA@var.features
multicores(mem = 200)
Feature_rast(PulmT_Satija, 'annotation.l1')
PulmT_Satija %<>%  ScaleData(verbose = FALSE, assay = 'RNA',feature = PulmT_Satija@assays$RNA@var.features,

                       vars.to.regress = c('percent.mito',  'nCount_RNA','nFeature_RNA', 'donor' ) )


table(PulmT_Satija$dataset_origin, PulmT_Satija$disease)


fwPulmT_Satija  %<>%  RunPCA(npcs =  100) %>% RunUMAP(dims = 1:60, 
                                                    reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:60)) %>% 
  FindClusters(resulution = 0.6)

Feature_rast(PulmT_Satija, c('ident', 'tissue', 'disease'), colorset = 'gg')

saveRDS(PulmT_Satija, 'PulmT_Satija_cleaned_processed.rds')

table(PulmT_Satija$donor, PulmT_Satija$disease)


ClusterCompare(PulmT_Satija %>% subset(annotation.l1 == 'CD8 T'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease')

PulmT_Satija@assays$RNA@scale.data

table(PulmT_Satija$donor)

table(PulmT_Satija$tissue, PulmT_Satija$disease)

PulmT_Satija$project <- 'satija'



# integration -------------------------------------------------------------

CD4CD8@assays
CD4CD8$donor <- CD4CD8$patient
CD4CD8$project <- 'Prinz&Falk'
CD4CD8$dataset_origin <- 'Prinz&Falk'


CD4CD8$pj_tissue <- paste(CD4CD8$project, CD4CD8$tissue)
PulmT_Satija$pj_tissue <- paste(PulmT_Satija$project, PulmT_Satija$tissue)

# CD4CD8.l <- CD4CD8 
# CD4CD8.l  <- AddMetaData(CD4CD8.l, t(CD4CD8.l[['CITE']]@data)) 

CD4CD8  <- AddMetaData(CD4CD8, t(CD4CD8[['CITE']]@data)) 



CD4CD8.l <- SplitObject(CD4CD8, split.by = "patient")
# CD4CD8.l  %<>%  map( ~ NormalizeData(.x) %>% FindVariableFeatures)

CD4CD8.l <-   lapply(X = CD4CD8.l, FUN = function(x) {
  # x  <- AddMetaData(x, t(x[['CITE']]@data)) 
  x[['CITE']] <- NULL
  x[['HTO']] <- NULL
  x[['GM']] <- NULL
  x[['integrated']] <- NULL
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
  
for (x in names(CD4CD8.l))      {
  
CD4CD8.l[[x]]@assays$RNA@var.features %<>% str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T) 
  
}


PulmT_Satija.l  <- SplitObject(PulmT_Satija, split.by = "dataset_origin")

map(PulmT_Satija.l, ~ dim(.x))

PulmT_Satija.l$mayr_2020 <- NULL
PulmT_Satija.l$lukassen_2020 <- NULL


PulmT_Satija.l <-   lapply(X = PulmT_Satija.l, FUN = function(x) {
  # x  <- AddMetaData(x, t(x[['CITE']]@data)) 

  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

for (x in names(PulmT_Satija.l))      {
  
  PulmT_Satija.l[[x]]@assays$RNA@var.features %<>% str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T)  %>% 
    setdiff(batcheff_copd_normal)
  
}


multicores(mem = 400)

anchors <- FindIntegrationAnchors(object.list = append(CD4CD8.l,PulmT_Satija.l ), 
                                  dims = 1:50)
Pulm.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)


Pulm.integrated$project


DefaultAssay(Pulm.integrated) <- "integrated"

Pulm.integrated@meta.data  %<>% mutate(dataset_origin = case_when(project == 'Prinz&Falk' ~ 'Prinz&Falk',
                                                                  project == 'satija' ~ dataset_origin)) 
# scale data and run PCA
Pulm.integrated %<>% ScaleData( vars.to.regress = c('donor', 
                                                    'dataset_origin',
                                           "percent.mito",
                                           # 'project',
                                           "percent.ribo",
                                           'nCount_RNA',
                                           'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)



ClusterCompare(Pulm.integrated, 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease')

ClusterCompare(PulmT_Satija, 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', rm = 'NO-NAME|RP|MT|HIST|HSPA', genetoshow = 100)

ElbowPlot(Pulm.integrated,ndims = 100)


Feature_rast(PulmT_Satija, 'disease')


Pulm.integrated <- RunUMAP(object = Pulm.integrated, dims = 1:70, 
                  reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:70))
for (i in seq(0.6,1,0.1) %>% rev()) {
  Pulm.integrated <- FindClusters(Pulm.integrated, resolution = i)
}

ClusterCompare(Pulm.integrated, '3', '7', do.plot = F)


Feature_rast(Pulm.integrated, c('tissue', 'CD4CD8','pj_tissue', 'disease', 'project'), ncol = 2, sz  = 0.2)

Feature_rast(Pulm.integrated, 'ident', colorset = 'gg')

Pulm.integrated$donor %>% unique()


Pulm.integrated$tissue %>% unique()

Pulm.integrated$CD4CD8


Feature_rast(Pulm.integrated, c('CD4', 'CD8A', 'CD8B', 'KLRB1', 'KLRG1', 'ITGA1',  'TRGC1', 'TRDC', 'FOXP3',
                                'ITGAE', 'CD69' ,'RORC', 'CCR6','PDCD1', 'TOX') , assay = 'RNA')



ClusterCompare(Pulm.integrated, '3', '4')
saveRDS(Pulm.integrated, 'PulmT.integrated.prinz.satija.rds')

# harmony -----------------------------------------------------------------
library(harmony)
Pulm.integrated@meta.data  %<>% mutate(dataset_origin = case_when(project == 'Prinz&Falk' ~ 'Prinz&Falk',
                                                                  project == 'satija' ~ dataset_origin)) %>% 
  `rownames<-`(Pulm.integrated$bc_backup)
DefaultAssay(Pulm.integrated) <- 'RNA'

Pulm.integrated  %<>%  FindVariableFeatures(selection.method = 'vst', nfeatures = 5000)



Pulm.integrated@assays$RNA@var.features %<>%str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T)  %>% 
  setdiff(batcheff_copd_normal)



# export to H5ad ----------------------------------------------------------
library(SeuratData)
library(SeuratDisk)



 



# GDTlung_trimmed$Cell_cluster <- as.vector(GDTlung_trimmed$Cell_cluster )
CD4CD8@meta.data  %<>%  mutate_if(is.factor, as.character) %>%  `rownames<-`(CD4CD8$bc_backup)


CD48trim <- CD4CD8 %>% DietSeurat(dimreducs = "umap", assays = "RNA")


CD48trim@assays$RNA$scale.data <-  NULL




SaveH5Seurat(CD48trim, 'ABTlung_pyscenic/CD4CD8', overwrite = T)
SaveH5Seurat(subset(CD48trim, CD4CD8 == "CD4")%>% DietSeurat(dimreducs = "umap"), 'abt/CD4', overwrite = T)
SaveH5Seurat(subset(CD48trim, CD4CD8 == "CD8")%>% DietSeurat(dimreducs = "umap"), 'abt/CD8', overwrite = T)

Convert('ABTlung_pyscenic/CD4CD8.h5seurat', dest = 'h5ad', overwrite = T)
write.csv(CD4CD8@meta.data, "ABTlung_pyscenic/ABTmeta.csv")
write.csv(FetchData(CD4CD8, c("UMAP_1", "UMAP_2")), "ABTlung_pyscenic/ABTumap.csv")



Convert('abt/CD4.h5seurat', dest = 'h5ad', overwrite = T)
Convert('abt/CD8.h5seurat', dest = 'h5ad', overwrite = T)

rm(CD48trim)

# Convert('GDTlung.trimmed.h5seurat', dest = 'h5ad')



dim(CD4CD8@assays$RNA@scale.data)



Pulm.integrated %<>% ScaleData( vars.to.regress = c(
  # 'donor', 
                                                    "percent.mito",
                                                    'dataset_origin',
                                                    # 'G2M.Score',
                                                    "percent.ribo",
                                                    'nCount_RNA',
                                                    'nFeature_RNA' ), assay = 'RNA', feature = Pulm.integrated@assays$RNA@var.features) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40) 


Pulm.integrated %<>%  RunHarmony(group.by.vars = c('dataset_origin', 'donor')) %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(preduction = "harmony", dims = 1:40) %>% FindClusters(resolution = 0.6)
ElbowPlot(Pulm.integrated, reduction = 'harmony',ndims = 100)

Feature_rast(Pulm.integrated, c('tissue', 'CD4CD8','pj_tissue', 'disease', 'project'), ncol = 2, sz  = 0.2)

  
Feature_rast(Pulm.integrated, c('disease', 'dataset_origin', 'Cell_pheno'))


ClusterCompare(Pulm.integrated, '1', '6')


GDTlung_s <- readRDS('../GDTlung280622.rds')


ClusterCompare(GDTlung_s, 'L1', 'P8')



# figures for slides ----------------------------------------------

install.packages("gtools")

library(gtools)

CD4CD8@meta.data %<>% mutate(
  Cluster.no=case_when(
    grepl("_M", Cell_cluster) ~ paste0("M",  str_extract(cluster_no, "(?<=c)\\d+") ),
    grepl("_L", Cell_cluster) ~ paste0("L",  str_extract(cluster_no, "(?<=c)\\d+") ),
    grepl("_P", Cell_cluster) ~ paste0("P",  str_extract(cluster_no, "(?<=c)\\d+") )
    )
) %>% mutate(Cluster.no = factor(Cluster.no, mixedsort(unique(Cluster.no)) ))

Feature_rast(CD4CD8, "Cluster.no",noaxis = F)
Feature_rast(GDTlung_s, "Cluster.no",noaxis = F)

Feature_rast(CD4CD8, "tissue",noaxis = F)

UMAP_ID_CD4CD8<- Feature_rast(CD4CD8, 'ID', colorset = alpha((ID_cl), 0.8), 
                       
                       do.label = F, noaxis = T, sz = 0.2)+NoLegend()+ggtitle('patient and tissue') 

UMAP_ID_CD4CD8

UMAP_ID<- Feature_rast(GDTlung_s, 'ID', colorset = alpha(c(brewer.pal(9,'Blues')[2:9], brewer.pal(9,'Reds')[2:9]), 0.8), 
                       
                       do.label = F, noaxis = T, sz = 0.2)+NoLegend()+ggtitle('patient and tissue') 

# ggsave(plot=UMAP_ID,filename =  "figs/GDlung_id.png", bg ="transparent", width = 40, height = 40)




# Accurate identification of tissue resident T cells (Trm) 










(Feature_rast(CD4CD8_cite,
              g = 'ID',
              # sz=0.2,
              mythe =F,
              # facets = c( 'tissue'),
              d1 ='CD49a.protein', d2 =  'CD103.protein',
              colorset = c(brewer.pal(9,'Blues')[6:9], brewer.pal(9,'Reds')[6:9]),
              # colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 6) %>% rev(),
              do.label = F,
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 1.2, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
  # theme(axis.text = element_text(size = 10))
)



(Feature_rast(CD4CD8_cite,
              g = 'ID',
              # sz=0.2,
              mythe =F,
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = c(brewer.pal(9,'Blues')[6:9], brewer.pal(9,'Reds')[6:9]),
              
              do.label=F,
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.6, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
)
(Feature_rast(CD4CD8_cite %>% subset(patient != 'p25'),  g='ID',  d1='CD4.protein', d2='CD8.protein', assay = 'CITE', mythe = F,
              colorset = c(brewer.pal(9,'Blues')[6:9], brewer.pal(9,'Reds')[6:9]),
              do.label = F,
              noaxis = F,  axis.number = T) +
    # facet_wrap(~ patient,ncol = 2)+
    ggtitle('CITEseq on CD4 and CD8')+geom_hline(yintercept = 1)+geom_vline(xintercept = 1) )


Feature_rast(CD4CD8, g = "Cell_pheno",sz = 0.5, noaxis = F, mythe =F,othertheme = NoLegend())

Idents(CD4CD8) <- CD4CD8$Cell_pheno

# Gene modules
CD4CD8$Cell_pheno %>%  unique()

CD4CD8 %>%  subset(subset = c("TRM_1_P", "TRM_2_P") %in% Cell_pheno ) %>% dim()

CD4CD8 %>%  subset(subset =  Cell_pheno != "unidentified_P" ) %>% dim()



ViolinPlot(CD4CD8 %>%  subset(Cell_pheno != "unidentified_P"), c( "Effectors", "Tissue.resident","CD8.Cytotoxictiy", "Th17"), colors = umap.colors, ncol = 2, sz = 0.2, box = T ,mythe = F,x.angle = 90, size = 10, ylabtext = "Score")


Feature_rast(CD4CD8, c("ITGAE", "ZNF683", "KLRD1", "PRF1",
                       "KLRB1", "CCR6", "DPP4"))

heatDOT_top5TF_TRM <-( DotPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl),dot.scale = 3.5,
                              features = rev(unique(top5TF_TRM$gene)))+mytheme+heattheme+
                        
                        theme(text = element_text(size = 8), axis.text.y = element_text(size = 8),
                              axis.line.y.right = element_line(),
                              axis.text.x = element_text(size = 8, angle= 90),
                              legend.box.margin = margin(5,0,0,15,unit = 'mm'),
                              legend.box = "horizontal",legend.position = 'bottom',
                              axis.title = element_blank())+coord_flip()+
                        scale_y_discrete(position = 'right')+
                        scale_x_discrete(position = 'top')+
                        xlab(NULL)+ylab(NULL)+
                        scale_color_gradient2(low = '#003399', mid = '#ffccff',  high = "#990000")+
                        guides(
                          color = guide_colorbar(title.position = 'top',direction = 'horizontal',
                          ),
                          size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom'))) 

heatDOT_top10_TRM

heatDOT_top5TF_TRM

Feature_rast(CD4CD8, c("KLRD1","KLRB1", "FCGR3A", "CCR6",  "PRF1",
                        "DPP4"), ncol = 2, sz = 0.2, mythe = F, 
             color_grd = "threecolor",
             othertheme = NoLegend())



Feature_rast(CD4CD8, "cdr3_TRB_perc", 
             color_grd = "threecolor")




Feature_rast(CD4CD8, 'clonal_expansion',do.label = F,
             mythe = F,sz = 0.5,
             colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'))


M_cls <- c("TEM_M", "TCM_1_M", "TCM_2_M","Temra_1_P", "TRM_1_P", "TRM_2_P", "TRM_3_P"
           )


Mori_result_CD8 <- repOverlap(Total_list1_CD8[!(names(Total_list1_CD8) %in% c("unidentified_P",
                "Naive_2_L",                                                 "Treg_L", "Tfh_L"))], .col='nt',
                              .method = "morisita", .verbose = F)

vis(Mori_result_CD8)+ 
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index of CD8 TCR')+xlab('Cluster')+ylab('Cluster')

Total_list1_CD4

Mori_result_CD4 <- repOverlap(Total_list1_CD4[!(names(Total_list1_CD8) %in% c("unidentified_P", "TRM_2_P", "TRM_3_P",  "Naive_1_L", "Treg_L",    "Temra_2_P"))], .col='nt',
                              .method = "morisita", .verbose = F)

vis(Mori_result_CD4)+ 
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index of CD4 TCR')+xlab('Cluster')+ylab('Cluster')




GDTlung_s$pheno

GDTlung_s$clonal_expansion


Feature_rast(GDTlung_s, "clonal_expansion",
             colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'),
             do.label = F,
             noaxis = F, mythe=F, navalue = "transparent")


paste0(c("CD103", "CD49a", "KLRG1",  "CD8"), '.protein') %>% 
  Feature_rast(GDTlung_cite, ., assay = 'CITE',
               titlesize = 12, titleface = "plain",
               sz = 0.2, ncol =2, mythe = F)


Feature_rast(GDTlung_s, "pheno")

c("KLRB1",  "CCR6",  "DPP4") %>% 
Feature_rast(GDTlung_s,., sz = 0.5,color_grd = "threecolor",
             ncol = 1, mythe = F, titlesize = 12)

c( "NKG7", "FCGR3A", "EOMES") %>% 
  Feature_rast(GDTlung_s,., sz = 0.5,color_grd = "threecolor",
               ncol = 1, mythe = F, titlesize = 12)

c( "CTLA4","CSF1",  "ENTPD1","AREG", "TNFRSF18","RBPJ") %>% 
  Feature_rast(GDTlung_s,., sz = 0.5,color_grd = "threecolor",
               ncol = 2, mythe = F, titlesize = 12,
               othertheme = NoLegend())



  Feature_rast(GDTlung_s, "cdr3_TRD_perc", facets = "patient",,
             navalue = "transparent")



Feature_rast(CD4CD8, "cdr3_paired_perc", facets = "patient",,
             navalue = "transparent")


TRM_1_3_DEG <- ClusterCompare(CD4CD8, "TRM_1_P", "TRM_3_P")


TRM_1_3_DEG$table %>%  filter(grepl("GZ", gene) & avg_log2FC > 0)


TRM_1_Th17_DEG <-  ClusterCompare(CD4CD8,
                                  "Th17_M", "TRM_1_P") %T>% print(.$plot) 

TRM_1_Th17_DEG$plot




# Figure for manuscript --------------------------------------------------------


# Fig 5 and S1 ------------------------------------------------------------
# A Workflow
# B and D Umap abT and Umap CD4CD8 

# F5A T cells sorting

# F5C
F5A <- (Feature_rast(CD4CD8, sz = 0.3,
                             labelsize = 4,
                             noaxis = F, othertheme = list(coord_fixed(),theme(legend.spacing.x = unit(1, 'mm'),
                                                                               legend.margin = margin(0,0,5,0, "mm")))
                             
)+ggtitle('abT cells')+guides(color = guide_legend(ncol = 3,position = "bottom",
                                                      
                                                   override.aes = list(size = 1.5)))
) %T>% print() 
F5A



theme = theme(legend.key.size = unit(0.5, "mm") ,legend.spacing.y = unit(1, 'mm'))

CD4CD8_cite@assays$CITE@counts %>%  rownames()

# F5B15 <-  map( c("CD103.protein", "CD49a.protein", "CD45RA.protein", "CD27.protein" , "CD26.protein"), ~ 
#                  Feature_rast(CD4CD8_cite, ., colorgrd = "grd2",     sz = 0.2,othertheme = list(theme(
#                    legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())   )        
# )
# 
#  map( c("CD103.protein", "CD45RA.protein", "CD27.protein" , "CD26.protein", "KLRG1.protein", "CD8.protein"), ~ 
#                  Feature_rast(CD4CD8_cite, ., colorgrd = "grd2",     sz = 0.2,othertheme = list(theme(
#                    legend.margin = margin(0,0,0,-10, "pt")), 
#                    scale_m("grd2",  limits = c(0,3)),
#                    coord_fixed())   )        
# ) %>% PG()
# 
# F5B15[[5]]
# 
# 
# F5B6 <-  Feature_rast(CD4CD8_cite, "KLRG1.protein", colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
#   legend.margin = margin(0,0,0,-10, "pt")),
#   scale_m("grd2", c(0,1, 2), c(0,2)),
#   coord_fixed())   )  
# F5B6
# F5B15[[6]] <- F5B6
# 
# F5B <- PG(F5B15, ncol = 3)
# F5B
# 
# 
# Feature_rast(CD4CD8_cite, c("CD4.protein", "CD8.protein"), colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
#   legend.margin = margin(0,0,0,-10, "pt")),
#   scale_m("grd2", c(0,1, 2), c(0,2)),
#   coord_fixed())   )  


Feature_density(CD4CD8,c("AREG", "IFNG", "GZMB", "GZMA"))


# F5C <-  Feature_rast(CD4CD8, "ID", colorset = ID_cl, do.label = F, sz = 0.3,  noaxis = F,
#                      othertheme = list( NoLegend(),coord_fixed())) %T>% print()


F5B <- ggplot(comp_tissue, aes(y = n, x = Cell_pheno, fill = ID, color = ID,
                                stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = ID_cl  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell no. of\ncluster per donor')+
  # xlab('cluster')+ 
  xlab(NULL)+
  ggtitle("ab T cell tissue distribution")+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  mytheme+
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    
    axis.line = element_blank())+
  # NoLegend()
NULL


F5C1 <-  map( c("CD103.protein", "CD49a.protein",  "CD45RA.protein", "CD27.protein" , "CD26.protein","CD94.protein"), ~ 
               Feature_rast(CD4CD8_cite, ., colorgrd = "grd2",     sz = 0.1,othertheme = list(theme(
                 legend.position = "none"), 
                 coord_fixed())   )        
) %>% PG(ncol = 3) %T>% print()


F5C2 <-  map( c("KLRG1.protein", "CD4.protein" , "CD8.protein"), ~ 
                Feature_rast(CD4CD8_cite, ., colorgrd = "grd2",     sz = 0.1,othertheme = list(theme(
                  legend.margin = margin(0,0,0,-10, "pt")), 
                  scale_m("grd2",  limits = c(0,3.5)),
                  coord_fixed())   )        
) %>% PG(ncol = 3) 
F5C2


F5C <-  PG(list(F5C1, F5C2), ncol = 1, rh = c(2,1))


Feature_rast(CD4CD8_cite,c("CD103.protein", "CD49a.protein",  "CD45RA.protein"), sz = 0.2, colorgrd = "grd2" )  %>% figsave("F5test.pdf", 150,50)

figsave(F5C, "F5C.pdf", 150, 150)


F5ABC <-  PG(list(F5A, F5B, F5C), ncol = 3  , labels = "AUTO", rw = c(1, 1.4, 1.2) )  %T>% print()


F5D <- (Feature_rast(CD4CD8, sz = 0.3,
                     g= "CD4CD8", do.label = F,
                     othertheme = coord_fixed(),
                     noaxis = F
                     
)+ggtitle('CD4+ and CD8+ T cell distribution')+
  theme(legend.position = c(0.05, 00.05), legend.justification = c(0.05, 0.05))


) %T>% print() 







GMS_long <- FetchData(CD4CD8, c("CD4CD8", "Cell_pheno", names(sigtable))) %>% 
  reshape2::melt(value.name = "score" ) %>% 
  mutate(variable = if_else(variable =="CD8.Cytotoxictiy","Cytotoxicity",variable ) )

names(sigtable)

sigtable$Th1

head(GMS_long)
GMS_long$variable %>% unique()

nrow(GMS_long)

GMS_long %>% group_by(variable) %>% count()

GMS_long %>% group_by(CD4CD8, Cell_pheno, variable) %>% count() %>% filter(n < 150) %>% ungroup() %>% 
  distinct(CD4CD8, Cell_pheno)

GMS_long  %<>% 
  filter(!((CD4CD8 == 'CD4' &  Cell_pheno  %in% c("TRM_2_LG", "TRM_3_LG")) | (CD4CD8 == "CD8" & Cell_pheno == "Treg_LN")))

GMS_long  %>% 
  filter((CD4CD8 == 'CD4' &  Cell_pheno  %in% c("TRM_2_LG", "TRM_3_LG")) | (CD4CD8 == "CD8" & Cell_pheno == "Treg_LN"))


F5E <- (ggplot(GMS_long %>% filter(variable %in% c( "Tissue.resident", "Th17", "Cytotoxicity") & Cell_pheno != "unidentified_LG" ) %>% 
                 mutate(variable = factor(variable, level = c( "Tissue.resident", "Th17", "Cytotoxicity"))), 
               aes(x = Cell_pheno, y = score, fill = CD4CD8))+
          geom_violin_rast(size = 0.1) +
          facet_wrap(~variable, ncol = 1, strip.position = "right", scales = "free_y")+
          ylab("module score")+xlab(NULL)+
          geom_boxplot( alpha = 0.5, size = 0.1,notch = F,
                        outlier.alpha = 0,show.legend = FALSE)+
          
          fill_m()+
          theme_minimal()+
          mytheme+
          # guides(fill = guide_legend(title = "type", keyheight = unit(4,"mm"), keywidth  = unit(3,"mm"), label.position = "bottom")
          #   
          # )+
          NoLegend()+
          theme(axis.text.x = element_text(size = 8 , angle = 90))
        
) %T>% print() 

Feature_density(CD4CD8, c("IFNG", "GZMB","AREG"))





F5F <- ggplot(gini_ABTpaire, aes(x = Cell_pheno , y = Gini_Index, 
                          color = patient, group= Cell_pheno))+geom_boxplot(show.legend = FALSE, outlier.colour = "transparent")+
  geom_point(size = 0.5)+facet_wrap(~CD4orCD8,ncol = 1)+
  color_m(color = RColorBrewer::brewer.pal(n = 7,name = "Accent"))+theme_minimal()  +
  theme(axis.text.x = element_text(angle = 90),
        legend.margin = margin(0,0,0,-10, "pt"))+mytheme+
  guides(color = guide_legend(override.aes = list(size = 1)))


F5DEF <-  PG(list(F5D, F5E, F5F), ncol = 3  , labels = c("D", "E", "F"), rw = c(1, 1.2, 1) )  %T>% print() 



# PG(list(F2A, F2B), ncol = 2, rw =c(1, 0.6)) %T>%  print()

# F5G TFs   



PG(list(F5ABC, F5DEF), ncol = 1)


F5G_velo <-  NA

CD4_TCR_sharing <-  plot_tcr_sharing(CD4_TRB_data, title = "CD4+ TCR Sharing Across Clusters",size_range = c(0.1, 3)) + coord_fixed()+gglp(p = "b")


CD8_TCR_sharing <-  plot_tcr_sharing(CD8_TRB_data, title = "CD8+ TCR Sharing Across Clusters", size_range = c(0.1, 3)) + coord_fixed()+gglp(p = "b")

F5H_TCR <-  PG(list(CD4_TCR_sharing,CD8_TCR_sharing ), align = "hv") %T>% print()

F5GH <-  PG(list(F5G_velo, F5H_TCR), labels = c("G", "H"), rw = c(1,1)) %T>% print()

# F5I_reg 
regulons <- c("RORA-REG", "MAF-REG", "RUNX2-REG", 
              
              "EOMES-REG", "TBX21-REG", "RUNX3-REG")



F5I_reg <-  ViolinPlot(CD4CD8 %>%  subset(Cell_pheno %in% c("TRM_1_LG" , "TRM_3_LG")), box = T,  sz = 0.2,
           x.angle = 330, ylabtext = "value",  ncol = 3, g = regulons,
           othertheme = list(stat_compare_means(  paired = F, method = 'wilcox.test',label = "p.format" ),
                             theme(axis.text.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 6))),
           assay = "AUC", colors = umap.colors[c(12,14)]) %T>%  print()


F5J_genes  <-  Feature_density(CD4CD8, c("RORA", "RUNX2",  "EOMES", "RUNX3" , "GATA3", "ZNF683",
                                         "CD40LG", "IL17A",  "CTLA4", "GZMA",  "GZMB", "AREG" )  , ncol = 4,
                               othertheme =   list(theme(
                                 legend.position = "none"),
                                 coord_fixed()
                                 
                               ))  %T>% print()
  

F5IJ <-  PG(list(F5I_reg, F5J_genes), labels = c("I", "J"), rw = c(1,1)) %T>% print()
F5 <-  PG(list(F5ABC, F5DEF, F5GH, F5IJ), ncol = 1, rh = c(2.4,2.2,2,2)) %T>% print()

figsave(F5, "figure5_new_AbT_2025_v2.pdf",200, 280)




F5 <-  PG(list(F5ABC, F5DEF, F5G, F5I), ncol = 1, rh = c(2.4,2.2,1,1.5)) %T>% print()

figsave(F5, "figure5_new_AbT_2025_v2.pdf",200, 230, path = figpath_ni)



F5G <-  Feature_density(CD4CD8, c("TBX21", "EOMES", "GATA3",  "ZNF683", "RORC", "RUNX2", "RUNX3")  , ncol = 7,
                     othertheme =   list(theme(
                       legend.position = "none"),
                       coord_fixed()
                       
                     ))  %T>% print()

trmcl2 <- c(
  "TEM_M",
  "Th17_M",
  "Temra_1_LG",
  "TRM_1_LG",
  # "TRM_2_LG",
  "TRM_3_LG"
  # "Temra_2_P"
  
)
genef22  <-  c( "CD40LG", "IL17A",  "CTLA4", "GZMA", "GZMK", "GZMB", "AREG")
Feature_density(CD4CD8, genef22  , ncol = 7,
                othertheme =   list(theme(
                  legend.position = "none"),
                  coord_fixed()
                  
                ))  %T>% print()

F5I <- ViolinPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl2), genef22, colors = umap.colors, 
                   box = F, ncol =7, sz = 0.2, othertheme = theme(axis.text.x = element_text(angle = 70, vjust = 0))
)  %T>% print()


F5 <-  PG(list(F5ABC, F5DEF, F5G, F5I), ncol = 1, rh = c(2.4,2.2,1,1.5)) %T>% print()

figsave(F5, "figure5_new_AbT_2025.pdf",200, 230, path = figpath_ni)



trmcl <- c(
  
  "TEM_M",
  "Temra_1_LG",
  "TRM_1_LG",
  "TRM_2_LG",
  "TRM_3_LG"
  # "Temra_2_P"
  
)



# Feature_rast(CD4CD8)
# 
# F2D <-( DotPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl),dot.scale = 3.5,
#                 features = rev(c(unique(top10TF_TRM$gene), "TOX", "RORC"))
# )+mytheme+heattheme+
#   
#   theme(text = element_text(size = 8),
#         axis.text.y = element_text(size = 8),
#         axis.line.y.left = element_line(),
#         axis.text.x = element_text(size = 8, angle= 90,face = "italic"),
#         legend.box.margin = margin(5,0,0,5,unit = 'mm'),
#         legend.box = "horizontal",legend.position = 'bottom',
#         axis.title = element_blank())+
#   # coord_flip()+
#   scale_y_discrete(position = 'right')+
#   scale_x_discrete(position = 'bottom')+
#   xlab(NULL)+ylab(NULL)+
#   ggtitle("Transcription factors")+
#   viridis::scale_color_viridis(discrete = F, option ="D")+
#   guides(
#     color = guide_colorbar(title.position = 'top',direction = 'horizontal',
#     ),
#     size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom'))+
#   NoLegend()
# )  %T>% print()
# 
# 
# 
# F2C <-  DoHeatmap(CD4CD8 %>% subset(downsample = 500,Cell_pheno %in% trmcl), features = top5_48_reg, assay = "AUC", size = gs(8)) %>% heat_theme() +mytheme+heattheme 
# 
# 
# F2C
# 
# CTG <- CD4CD8_DEGs %>%  filter(gene %in% Cytlist & cluster %in% trmcl) %>% 
#   group_by(gene) %>% 
#   top_n(1, avg_log2FC) %>% 
#   arrange(cluster, desc(avg_log2FC))%>%
#   pull(gene)
# 
# 
# 
# CD4CD8_DEGs %>%  filter(gene %in% Cytlist & cluster %in% "TRM_1_LG")
# 
# 
# F2D <- ( DotPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl),dot.scale = 3.5,
#                  features = CTG)+mytheme+heattheme+
#            
#            theme(text = element_text(size = 8),
#                  axis.text.y = element_text(size = 8),
#                  axis.line.y.left = element_line(),
#                  axis.text.x = element_text(size = 8, angle= 90,face = "italic"),
#                  legend.box.margin = margin(5,0,0,5,unit = 'mm'),
#                  legend.box = "horizontal",legend.position = 'bottom',
#                  axis.title = element_blank())+
#            ggtitle("Cytokines and Granzymes")+
#            
#            scale_y_discrete(position = 'right')+
#            scale_x_discrete(position = 'bottom')+
#            xlab(NULL)+ylab(NULL)+
#            viridis::scale_color_viridis(discrete = F, option ="D")+
#            guides(
#              color = guide_colorbar(title.position = 'top',direction = 'horizontal',
#              ),
#              size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom')))  %T>% print()
# 
# 
# 
# F2E <-  Feature_rast(CD4CD8,genef2d, sz = 0.2, ncol = 6,
#                      othertheme =theme(
#                        legend.margin = margin(0,0,0,-10, "pt"))) %T>%  print()
# 
# F2AB <- PG(list(F2A, F2B), ncol = 2, rw =c(1, 0.6), labels = "AUTO") 
# F2CDE <-  PG(list(F2C, F2E), ncol = 1, rh = c(0.8, 0.4),labels = c("C", "D") )
# F2 <- PG(list(F2AB,F2CDE), ncol = 1, rh=c(12, 16)) %T>% figsave("Fig2Test_2024.pdf", 200, 280) 
# 
# 
# trmcl2 <- c(
#   # "TEM_M",
#   "Th17_M",
#   "Temra_1_LG",
#   "TRM_1_LG",
#   # "TRM_2_LG",
#   "TRM_3_LG"
#   # "Temra_2_P"
#   
# )
# genef21 <- c("GZMA",   "GZMK",  "IFNG", "TNF") 
# genef22  <-  c( "CD40LG", "IL17A",  "CTLA4", "GZMA", "GZMB", "AREG", "TNF")
# 
# F2F5 <- ViolinPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl2), genef21, colors = umap.colors, box = F, ncol = 4,othertheme = theme(axis.text.x = element_blank()) ) %T>% print()
# 
# 
# ViolinPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl2), genef21, colors = umap.colors, box = F, ncol = 4,othertheme = theme(axis.text.x = element_blank()) ,sz =  0.2) %T>% print()
# 
# F2F2 <- ViolinPlot(CD4CD8 %>% subset(Cell_pheno %in% trmcl2), genef22, colors = umap.colors, 
#                    box = F, ncol =4, sz = 0.5, othertheme = theme(axis.text.x = element_text(angle = 70, vjust = 0))
# )  %T>% print()
# F2F <- PG(list(F2F5, F2F2), ncol = 1, rh = c(1,1.3))
# F2F
# 
# CD4CD8_DEGs %>%  filter(gene %in% Cytlist & cluster %in% "TRM_3_LG")
# 
# 
# F2 <- PG(list(F2ABC, F2D, NA, F2F), ncol = 1, labels = c(NA, "D", "E", "F"), rh = c(1, 0.35, 0.8, 1)) %T>% print() %T>% figsave("FIg2_2024.pdf",200, 270, path = figpath_ni) 




# PG(list(F5A, F5D, NA), ncol = 3)
# 
# PG(list(F5A, F5C, F5D), ncol = 3  , labels = "AUTO", rw = c(1.1, 1, 1.3) ) 
# 
# F5F
# 
# F5I <-  Feature_rast(GDTlung_s, "ID", colorset = ID_cl, do.label = F, sz = 0.3, othertheme = list( NoLegend(),coord_fixed()))  %T>% print()
# 
# 
# Feature_rast(GDTlung_s, "ID", colorset = ID_cl, do.label = F, sz = 0.3, othertheme = list( coord_fixed()))  %T>% print()
# 
# Feature_rast(CD4CD8, "ID", colorset = ID_cl, do.label = F, sz = 0.3) 
# 
# 
# F5J <- ggplot(cl_comp_gd, aes(y = n, x = Cell_cluster, fill = ID, 
#                                color = ID,
#                                stratum = ID  )) +
#   scale_x_discrete(expand = c(.1, .1)) +
#   geom_stratum(alpha=0.8 ,size = 0.2)+
#   scale_fill_manual(values = ID_cl,
#                     labels = rep(paste0('p', c(25,27,32,45,71,73,77)),2)
#   )+
#   scale_color_manual(values = ID_cl)+
#   theme_minimal() + 
#   ylab('Cell no. of\ncluster per donor')+
#   # xlab('cluster')+ 
#   xlab(NULL)+
#   ggtitle("gd T cell tissue distribution")+
#   scale_y_continuous(labels = abs)+
#   geom_hline(yintercept = 0,size = 0.5)+
#   guides(fill = guide_legend(nrow = 2, title = 'LLN\n\nLung',
#                              byrow = T, label.position = 'bottom',
#                              override.aes = list(color = "white")),
#          color = F )+
#   mytheme+
#   theme(legend.position = 'bottom', 
#         
#         legend.key.height  = unit(5, 'mm'), 
#         legend.key.width   = unit(5,'mm'), 
#         axis.text.x = element_text(size = 8, angle = 90),
#         legend.text = element_text(size=8),
#         axis.line = element_blank())+
#   NULL
# 
# F5J
# 
# lgd <-  get_legend(F5J) %>% ggplotify::as.ggplot()
# 
# F5EFIJ <- PG(list(F5E, F5F, F5I, (F5J+NoLegend()), lgd), labels = c("E", "F", "I", "J", NA),ncol = 2, rw = c(1,1.2), rh = c(1,1,0.15)) %T>% print()  
# 
# Feature_rast(GDTlung_s, sz = 0.3, "AREG",
#              labelsize = 8,
#              noaxis = F, legendcol = 2, othertheme = coord_fixed())
# 
# F5G <- (Feature_rast(GDTlung_s, sz = 0.3,
#                      labelsize = 8,
#                      noaxis = F, legendcol = 3, othertheme = coord_fixed()
#                      
# )+ggtitle('gdT cells')+gglp(p = "b"))%T>% print() 
# 
# F5H <-  Feature_rast(GDTlung_cite, c("CD103.protein", "CD49a.protein", "CD45RA.protein", "CD27.protein" , "IL7R.protein"), colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
#   legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())   ) %T>% print()
# 
# 
# 
# F5H15 <-  map( c("CD103.protein", "CD49a.protein", "CD45RA.protein", "CD27.protein" , "IL7R.protein"), ~ 
#                  Feature_rast(GDTlung_cite, ., colorgrd = "grd2",     sz = 0.2,othertheme = list(theme(
#                    legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())   )             
#                  )
# 
# F5H6 <-  Feature_rast(GDTlung_cite, "KLRG1.protein", colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
#   legend.margin = margin(0,0,0,-10, "pt")),
#   scale_m("grd2", c(0,0.5, 1,1.5), c(0,1.5)),
#   coord_fixed())   )
# 
# F5H15[[6]] <- F5H6
# 
# F5H <- PG(F5H15, ncol = 3)
# F5H
# 
# PG(F5H15, ncol = 2)
# 
# 
# F5GH <- PG(list(F5G, F5H), labels = c('G', 'H'), ncol = 1, rh = c(1.3,1))
# F5EFIJGH <- PG(list(F5EFIJ,F5GH), ncol = 2, rw = c(1.6,1) ) %T>% print()
# 
# 
# F5_new <- PG(list(F5AB, F5CD, F5EFIJGH), ncol = 1, rh = c(0.6,1,2)) %T>% print() %T>% figsave("Fig1_2024_new_design_1.pdf", 200, 260, path = figpath_ni)
# 
# F5_new
# supplementary

FS1A <-  
  (Feature_rast(CD4CD8_cite,
                g = "T_pheno", sz = 0.3,
                facets = c('tissue'),  do.label = F,
                othertheme = list(theme_minimal()),
                d1 ='CD45RA.protein', d2 =  'CD27.protein',
                colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
                
                # slot = 'scale.data',
                noaxis = F, assay = 'CITE',
                axis.number = T)+ggtitle("Memory phenotype of ab T cells")+
     geom_hline(yintercept = 0.5, linewidth = 0.2)+
     geom_vline(xintercept = 1.2,linewidth = 0.2)
  ) %T>%  print()

FS1A


FS1C <-  
  (Feature_rast(GDTlung_cite, 
                g = "T_pheno", sz = 0.3,
                othertheme = (theme_minimal()),
                
                facets = c('tissue'),  do.label = F,
                d1 ='CD45RA.protein', d2 =  'CD27.protein',
                colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
                
                # slot = 'scale.data',
                noaxis = F, assay = 'CITE',
                axis.number = T)+ggtitle("Memory phenotype of gd T cells")+
     geom_hline(yintercept = 0.6, linewidth = 0.2)+
     geom_vline(xintercept = 2.5,linewidth = 0.2)
  ) %T>%  print()

FS1B <- (DoHeatmap(pseudobuklCD4CD8,features = top5deg_CD4CD8_bulk$gene,raster = T,  
                   assay = "RNA", 
                   group.colors = umap.colors,size = gs(8) ) +
           ggtitle("top5 DEG abT cells") )%>% 
  heat_theme() %T>% print()


cols.use <- list(Cell_pheno = umap.colors, tissue = c("#003399", "#990000"), patient =RColorBrewer::brewer.pal(n = 7,name = "Accent")
)

 DoMultiBarHeatmap(subset(CD4CD8,downsample = 500),features = top10deg$gene,  group.by='Cell_pheno', size = gs(6),additional.group.sort.by = "tissue", angle = 20,
                          additional.group.by = c("tissue", "patient", "CD4CD8"),
                          cols.use =cols.use) %>% heat_theme(legend.position = "right") %T>% print()



FS1D <- (DoHeatmap(pseudobuklGD,features = top5degGDT_bulk$gene,raster = T,  
                   assay = "RNA", 
                   group.colors = umap.colors,size = gs(8) ) +
           ggtitle("top5 DEG gdT cells") )%>% heat_theme() %T>% print()

PG(list(FS1A,FS1C, FS1B, FS1D), ncol = 2, rh = c(1,4), labels = c("A", "C", "B", "D")) %T>% print() %T>% figsave("Fig_S1_2024_1.pdf", 200, 270,path = figpath_ni) 
figpath_ni

F5


Feature_rast(CD4CD8, c( "EOMES-REG", "MAF-REG","TBX21-REG",   "RORC-REG", "RORA-REG"), colorgrd = "grd2", assay = "AUC")


Feature_rast(CD4CD8,  "EOMES-REG", assay = "AUC", colorgrd = "grd2", sz = 0.2,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  scale_m("grd2", c(0,0.1, 0.2), c(0,0.2)),
  coord_fixed())   )



Feature_rast(CD4CD8,  "RORA-REG", assay = "AUC", colorgrd = "grd2", sz = 0.2,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  # scale_m("grd2", c(0,0.05, 0.1,0.15), c(0,0.15)),
  coord_fixed())   )
Feature_rast(CD4CD8, facets = "CD4CD8", facetcol = 1, ,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  # scale_m("grd2", c(0,0.05, 0.1,0.15), c(0,0.15)),
  coord_fixed()) )


# Fig6 ABT TCR ------------------------------------------------------------

Feature_rast(CD4CD8, c("ident", "CD4CD8"),
             othertheme =  list( 
                                 coord_fixed() ), ncol = 1, noaxis = F)



F3B <- Feature_rast(CD4CD8, 'clonal_expansion', facets = ('CD4CD8'),do.label = F,
                    sz = 0.5, facetcol = 1,
                    othertheme =  list( theme(legend.position = c(0.05, 00.05), 
                                              legend.justification = c(0.05, 0.05)),
                                        coord_fixed() ),
                                        
               colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2')
                    
                    ) %T>%  print()


F3A <-  ggplot(gini_ABTpaire, aes(x = Cell_pheno , y = Gini_Index, 
                                                  color = patient, group= Cell_pheno))+geom_boxplot(show.legend = FALSE, outlier.colour = "transparent")+
  geom_point(size = 0.5)+facet_wrap(~CD4orCD8,ncol = 1)+
  color_m(color = set_sample(umap.colors, s = 22))+theme_minimal()  +
  theme(axis.text.x = element_text(angle = 90),
        legend.margin = margin(0,0,0,-10, "pt"))+mytheme+
  guides(color = guide_legend(override.aes = list(size = 1)))






F3AB <-  PG(list(F3A,F3B), ncol = 2, rw = c(1, 1.2), labels = c("A", "B"))

F3C <-    
  PG(ALLMOSTEXPFIGS[c(13, 16, 18,  20)], ncol = 4) %>% 
  list(lg) %>% 
  PG( rw = c(6,1))    




# Plot the heatmap using ggplot2
MCD4 <-  ggplot(reshape2::melt(Mori_result_CD4, id.vars = "feature"), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + theme_nothing()+
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  mytheme+gglp("r")+
  theme(axis.text.x = element_text(size = 8 , angle  = 90))+
  xlab(NULL)+ylab(NULL)+
  # theme(text = element_text(size = 2))+ 
  ggtitle('Morisita index of CD4 TCRs') 


M_cls <- c("TEM_M", "TCM_1_M", "TCM_2_M","Temra_1_Lu", "TRM_1_Lu", "TRM_2_Lu", "TRM_3_Lu"
)


Mori_result_CD8 <- repOverlap(Total_list1_CD8[!(names(Total_list1_CD8) %in% c("unidentified_Lu",
                                                                              "Naive_2_L",                                                 "Treg_L", "Tfh_L"))], .col='nt',
                              .method = "morisita", .verbose = F)



MCD8 <-  ggplot(reshape2::melt(Mori_result_CD8, id.vars = "feature"), aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + theme_nothing()+
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  mytheme+gglp("r")+
  theme(axis.text.x = element_text(size = 8 , angle  = 90))+
  xlab(NULL)+ylab(NULL)+
  # theme(text = element_text(size = 2))+ 
  ggtitle('Morisita index of CD8 TCRs') 




F3E <-  PG(list(MCD4, MCD8), ncol = 1)
F3E


F3ABE <-  PG(list(F3AB, F3E), labels = c(NA, "E") , rw = c(3,1.5)) %T>% print()



TRM1group <- c("TCM_1_M", "Th17_M","TRM_1_LG",   "TCM_3_M")

TCRby_cluster$pairedfreq
TCRby_cluster4 <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 & CD4CD8 == "CD4") %>%
  # select(cdr3_paired,Cell_pheno) %>%
  group_by(cdr3_paired,Cell_pheno,) %>%  summarise(pairedfreq = n()) %>% ungroup() %>%
  arrange(Cell_pheno, pairedfreq) %>%
  mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )


TCRby_cluster8 <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 & CD4CD8 == "CD8") %>%
  # select(cdr3_paired,Cell_pheno) %>%
  group_by(cdr3_paired,Cell_pheno,) %>%  summarise(pairedfreq = n()) %>% ungroup() %>%
  arrange(Cell_pheno, pairedfreq) %>%
  mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )

F3F <-(TCRby_cluster4 %>% 
         filter(Cell_pheno %in% TRM1group  ) %>% mutate(Cell_pheno = factor(Cell_pheno, TRM1group)) %>% 
         ggplot(
           aes( x = Cell_pheno , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
         geom_flow(stat = "alluvium",
                   size = 0.1,
                   color = "darkgray") +
         theme_minimal_hgrid()+
         # scale_fill_manual(values = rainbow(1100)[280:810])+
         # facet_wrap(~patient, scales = "free", ncol = 5)+
         geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
         xlab(NULL) +ylab("TCRab frequencies")+ggtitle('TCR TRM_1 lineage')+
         theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
         guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 

F3F
TRM3group <- c( "TRM_3_LG","TCM_3_M",  "Temra_1_LG")

F3G <-(TCRby_cluster8 %>% 
                         filter(Cell_pheno %in% TRM3group  ) %>% mutate(Cell_pheno = factor(Cell_pheno, TRM3group)) %>% 
                         ggplot(
                           aes( x = Cell_pheno , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
                         geom_flow(stat = "alluvium",
                                   size = 0.1,
                                   color = "darkgray") +
                         theme_minimal_hgrid()+
                         # scale_fill_manual(values = rainbow(1100)[280:810])+
                         # facet_wrap(~patient, scales = "free", ncol = 5)+
                         geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
                         xlab(NULL) +ylab("TCRab frequencies")+ggtitle('TCR TRM_3 lineage')+
                         theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
                         guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 
F3G


F3FG <-  PG(list(F3F, F3G), ncol = 1, labels = c("F","G"))
F3CDFG <- PG(list(NA, F3FG), ncol = 2, rw = c(1.8,1))

F3 <- PG(list(F3ABE, F3CDFG), ncol = 1, rh = c(1, 1.3)) %T>%
  figsave("Fig3Test.pdf", 200, 240, path = figpath_ni)

CD4CD8$length_TRB
ViolinPlot(CD4CD8, "length_TRB", box = T)

Feature_rast(CD4CD8, d1 = "Cell_pheno", d2 = "length_TRB")


F3AB <- PG(list(F3A,F3B), ncol = 2, rw = c(2, 1), labels = "AUTO")
F3AB %T>% figsave("Fig3ABtest.pdf", 200, 80)

F3C
F3DE <- PG(list(F3E,F3F), ncol = 1,  labels = c("D", "E"))

F3CDE <- PG(list(F3C, F3DE, NA), ncol =3, rw= c(0.9, 1, 0.5), labels = c("C", NA,NA))

F3 <- PG(list(F3ABE, F3CDE), ncol = 1, rh = c(1, 1.5)) %T>%
  figsave("Fig3Test.pdf", 200, 200, path = figpath_ni)


CD4CD8_meta %>%  filter(Cell_pheno %in% c("TRM_1_Lu", "TRM_3_Lu", "Naive_1_L", "Naive_2_L" ) & !is.na(cdr3_TRB)) %>%  
  ggplot(aes(x = Cell_pheno, y = length_TRB))+ geom_boxplot_jitter() +facet_grid(CD4CD8 ~ patient)


library(data.table)

fwrite(CD4CD8[["RNA"]]@data %>% as.data.frame(), "abt/Normalized_count_CD4CD8.csv", row.names = T)

fwrite(data.frame(FetchData(CD4CD8, c("UMAP_1", "UMAP_2")), CD4CD8@meta.data),
          "abt/CD4CD8lung_meta.csv", row.names = T)


fwrite(data.frame(FetchData(GDTlung_s, c("UMAP_1", "UMAP_2")), GDTlung_s@meta.data),
       "GDTlung_meta.csv", row.names = T)


CD4CD8@meta.data %>% filter(!is.na(cdr3_TRB)) %>% count(patient, CD4CD8)

Feature_rast(CD4CD8 %>% subset(CD4CD8 == "CD4"), "cdr3_TRB_perc", facets="patient", navalue = "transparent")

Feature_rast(CD4CD8 %>% subset(CD4CD8 == "CD4" & patient == "p25"),
             "cdr3_TRB_freq", navalue = "transparent", 
             colorgrd = rev(brewer.pal(n = 11, name = "Spectral")))

Feature_rast(CD4CD8 %>% subset(CD4CD8 == "CD4" & patient == "p25"),
             "cdr3_TRB_freq", navalue = "transparent", 
             colorgrd = rev(brewer.pal(n = 11, name = "Spectral")))




CD4CD8@meta.data %>% dplyr::filter(CD4CD8 == "CD8" & !is.na(cdr3_TRB)) %>%  count(patient, CD4CD8, tissue)



Feature_rast(GDTlung_s, c("GATA3", "EOMES", "TNFRSF58", "AREG", "CSF2", "GZMA"),
              ncol = 3, sz = 0.4)

Feature_rast(CD4CD8, c("GATA3", "EOMES", "TNFRSF58", "AREG", "CSF2", "GZMA"), 
             ncol = 3, sz = 0.4)

GDTlung_s@meta.data %>% mutate(JC = paste0(j_gene_TRG, c_gene_TRG)) %>% 
  count(JC)



