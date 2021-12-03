
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

library(harmony)
# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')


file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

#paralell computation
options(future.fork.enable = TRUE)
options(future.globals.maxSize= 400*1024^3)
future::plan(strategy = "multicore", workers = 40)
# future::plan(strategy = "sequential")
# ?parallelly::supportsMulticore
#single core
setwd("/home/big/tanlikai/Lung")

GDTlung_s <- readRDS('GDTlung_6p_Seurat_07062021.rds')

# #set the working directory   --------------------------------------------

orig_cl <-c("#33cc33", #CB1
            "#339933", #CB2
            "#660066", #PB1
            "#cc00cc"  #PB2
) 

Feature_rast(GDTlung_s)

# work dir ----------------------------------------------------------------
# GDTlung.rds <- '/home/big/tanlikai/Lung/GDTlung_3p_20210213.RDS'
# GDTlung <- readRDS(GDTlung.rds)
# GDT_2020AUG.rds <- '/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds'
# GDT_2020AUG <- readRDS(GDT_2020AUG.rds)



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
rawdata_4$`Antibody Capture` %>% glimpse()
rownames(rawdata_4$`Antibody Capture`)


rawdata_4$`Antibody Capture`[1:6, ] 

GDTlung_4_6 <- CreateSeuratObject(rawdata_4$`Gene Expression`, project = 'p71&p31&p27',
                                meta.data = t(rawdata_4$`Antibody Capture`))

GDTlung_4_6@meta.data

GDTlung_4_6[['HTO']] <- CreateAssayObject(rawdata_4$`Antibody Capture`[1:6,])
GDTlung_4_6[['CITE']] <- CreateAssayObject(rawdata_4$`Antibody Capture`[7:23,])



GDTlung_4_6$bc_backup <- rownames(GDTlung_4_6@meta.data)


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
  `rownames<-`(GDTlung_1$bc_backup)
GDTlung_1$patient <- 'p25'
dim(GDTlung_1)
GDTlung_1 <-GDTlung_1 %>%  subset( subset = tissue %in% c('lung', 'luLN'))

table(GDTlung_1$tissue,useNA = 'always')

table(GDTlung_1$tissue,useNA = 'always')


tissuehashtag <- Feature_rast(GDTlung_1, g = 'tissue', d1 = "X1_TotalSeqC_LuPa" , d2 = "X9_TotalSeqC_LN",
             noaxis = F, axis.number = T)+xlim(0,5000)+ylim(0, 10000)+
  geom_vline(xintercept = 100)+geom_hline(yintercept = 180)
figsave(tissuehashtag,'tissuehashtag.pdf', 150,150)

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


# Feature_rast(GDTlung_2_3, 'hash.ID', d1 = 'X3_TotalSeqC_Lu_453',d2 = 'X9_TotalSeqC_Lu_322', do.label = F,
#              noaxis = F, axis.number = T )+geom_hline(yintercept = 200)
# 
# Feature_rast(GDTlung_2_3, 'hash.ID',d2 = 'X1_TotalSeqC_LN_325', d1 = 'X4_TotalSeqC_LN_456', do.label = F,
#              noaxis = F, axis.number = T )

# 
# 
# # adjust 
# GDTlung_2_3$hash.ID <- as.character(GDTlung_2_3$hash.ID)
# GDTlung_2_3@meta.data <- GDTlung_2_3@meta.data %>% mutate(hash.ID = case_when(
#   hash.ID == 'Negative' &  X9_TotalSeqC_Lu_322 > 200 ~ '9-TotalSeqC-Lu-322',
#   is.character(hash.ID) ~ hash.ID
# )) %>%
#   `rownames<-`(GDTlung_2_3$bc_backup)
# table(GDTlung_2_3$hash.ID)

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

# QC ----------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


GDTlung <- list(GDTlung_1, GDTlung_2_3, GDTlung_4_6) %>% setNames(c('p1', 'p23', 'p456'))

GDTlung %<>% map(., ~ PercentageFeatureSet(.x, '^MT', col.name =  'percent.mito') %>% 
                   PercentageFeatureSet('^RP', col.name = 'percent.ribo') #%>%
                   # CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
)

GDTlung$p1$percent.mito

GDTlung$p23$Phase
QCvio_T <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                         colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>% 
  PG(ncol = 3)

QCvio_T


QCvio_P <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'patient',
                                           colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>% PG(ncol = 3)

QCvio_P


figsave(QCvio_T,'beforeQC_violin_3patient.pdf', 200, 200)

QC_scatter <- GDTlung %>%  map( ., ~Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', facet = 'tissue',
                                           noaxis = F, axis.number = T)+grd+facet_wrap(~tissue)+
                    scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
                    scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
                    geom_hline(yintercept = c(400,2000))+geom_vline(xintercept = c(1000,7000))) %>%
  PG(ncol = 3)
  
  
QC_scatter
figsave(QC_scatter,'beforeQC_scatter.pdf',600,250)

# clean 
GDTlung %<>% map(.,~ subset(.x, subset =  nCount_RNA %in% 800:7000 &
                    nFeature_RNA %in% 400:2000 &  percent.mito <25&
                  tissue %in% c('lung','luLN')) 
                 )
dim(GDTlung$p1)
dim(GDTlung$p23)
# cell number
map(GDTlung, ~ dim(.x))
#mean counts 
map(GDTlung,~ .x$nCount_RNA %>% mean)


QCvio_T_clean <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                           colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>%
  PG(ncol = 3) %T>% figsave('afterQC_violin_6patient.pdf',150,150)

QCvio_T_clean

# figsave(QCvio_T_clean,'afterQC_violin_6patient.pdf',150,150)

saveRDS(GDTlung,'GDTlung_6p_beforeintegration.rds')

L_GDTlung <- readRDS('GDTlung_6p_beforeintegration.rds')

# integration by harmony --------------------------------------------------
GDTlung<- L_GDTlung %>% reduce(.f = merge, project = 'lungGDT')


test <- list(GDTlung, GDTlung_s, GDT_2020AUG)%>% reduce(.f = merge, project = 'lungGDT')


L_GDTlung$p23@assays$RNA@key


GDTlung%<>% NormalizeData(verbose = FALSE, assay = 'RNA') %>% 
  FindVariableFeatures(selection.method = 'vst', nfeatures = 3500) 


GDTlung@assays$RNA@var.features %<>%str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T) #%T>% print(length())
GDTlung@assays$RNA@var.features

#%>% 
GDTlung%<>%   ScaleData(verbose = FALSE, assay = 'RNA',    
                        vars.to.regress = c('percent.mito',  'nCount_RNA','nFeature_RNA' )
                        ) %>% 
  RunPCA(npcs = 100, verbose=T)
GDTlung%<>% RunHarmony(c("patient", 'orig.ident', 'percent.mito',  'nCount_RNA','nFeature_RNA'), plot_convergence = TRUE)

ElbowPlot(GDTlung, ndims = 50,reduction = 'harmony')

#chose 35 PCs 





# dim reduce via harmony --------------------------------------------------------------


GDTlung %<>% RunUMAP(reduction = 'harmony', dims = 1:30) %>% FindNeighbors(redcution = 'harmony', dims = 1:30)

for (i in seq(0.6,1.5,0.1) %>% rev()) {
  GDTlung <- FindClusters(GDTlung, resolution = i)
}

GDTlung %<>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


Feature_rast(GDTlung, c('ident', 'tissue', 'patient'))
 
ViolinPlot(GDTlung, c('nCount_RNA','percent.mito', 'percent.ribo'), colors = umap.colors)

ClusterCompare(GDTlung, '3', '7')


# Since by harmoony  we failed to achieve clear UMAP, we decide to stay with Seurat protocol
saveRDS(GDTlung, 'GDTlung.harmony:integration.rds')
rm(GDTlung )



# integration by Seurat ---------------------------------------------------



# normalization -----------------------------------------------------------
L_GDTlung%<>% map(.,~ NormalizeData(.x, normalization.method = 'LogNormalize',
                          scale.factor = 10000, assay = 'RNA') %>%
  ScaleData( assay = 'RNA',
             vars.to.regress = c("percent.mito",
                                 # "S.Score", 'G2M.Score',
                                 "percent.ribo",
                                 'nCount_RNA','nFeature_RNA' )) %>%
  FindVariableFeatures(assay = 'RNA',nfeatures = 2500, selection.method = 'vst')
)
L_GDTlung$p456

for (i in c('p1','p23', 'p456')) {
L_GDTlung[[i]]@assays$RNA@var.features%<>%  
    str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T)
}





# data integration --------------------------------------------------------
anchors <- FindIntegrationAnchors(L_GDTlung, dims = 1:60)

GDTlung_s <- IntegrateData(anchors, dims = 1:60)

DefaultAssay(GDTlung_s) <- "integrated"


GDTlung_s %<>% ScaleData(vars.to.regress = c("percent.mito", 
                                             'patient',
                                                    # "S.Score", 'G2M.Score',
                                                    "percent.ribo",
                                                    'nCount_RNA','nFeature_RNA' ))

GDTlung_s %<>% ScaleData(assay = 'RNA', features = rownames(GDTlung_s),
                       vars.to.regress = c("percent.mito", 'patient',
                                           # "S.Score", 'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA','nFeature_RNA' )
                        )


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


GDTlung_s %<>%  CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


# PCA  --------------------------------------------------------------------
dim(GDTlung_s)

GDTlung_s %<>% RunPCA( npcs = 100, verbose =
                     T, nfeatures.print = 40) %>% 
  JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
JackStrawPlot(GDTlung_s, dims = 1:50 ) %T>%
  figsave('GDTlung_6p.jackstraw.pdf' , w = 400, h = 400)
# saveRDS(GDTlung,GDTlung.rds)



# dimensional reduction & clustering --------------------------------------

GDTlung_s  %<>%  RunUMAP( dims = 1:31, 
                    reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:31))

for (i in seq(0.5,1.2,0.1) %>% rev()) {
  GDTlung_s  %<>% FindClusters( resolution = i)
}
GDTlung_s$RNA_snn_res.1.2

GDTlung_s$ID <- paste0(GDTlung_s$tissue, GDTlung_s$patient)

Feature_rast(GDTlung_s,c('ident','tissue','patient', 'ID'), ncol = 2, noaxis = F)
ViolinPlot(GDTlung_s,'nFeature_RNA', group.by = 'integrated_snn_res.0.9')
DefaultAssay(GDTlung_s) <- "RNA"


Feature_rast(GDTlung_s,c('TRDV2','TRDV1','TRDV3','TRDC','TRGC1','TRGV9', "CD3D",'CD3E'), ncol =4)
Feature_rast(GDTlung_s, axis.number = T, noaxis = F)+geom_vline(xintercept = c(-5.3,7))
# trim non-T cells 
GDTlung_s %<>% subset( UMAP_1 > -5.3 )
Feature_rast(GDTlung_s,c('ident','tissue','patient', 'ID'), ncol = 2, noaxis = F)


# 2nd dim-reduc &clustering -----------------------------------------------

DefaultAssay(GDTlung_s) <- "integrated"



GDTlung_s %<>% ScaleData(vars.to.regress = c("percent.mito", 'patient',
                                           # "S.Score", 'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA','nFeature_RNA' ))

GDTlung_s %<>% ScaleData(assay = 'RNA', features = rownames(GDTlung),
                       vars.to.regress = c("percent.mito", 'patient',
                                           # "S.Score", 'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA','nFeature_RNA' )
)


GDTlung_s %<>%RunPCA( npcs = 100, verbose =
                    T, nfeatures.print = 40)


ElbowPlot(GDTlung_s, ndims = 100)
GDTlung_s %<>% JackStraw(num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 
JackStrawPlot(GDTlung_s, dims = 1:50 )
JackStrawPlot(GDTlung_s, dims = 1:50 ) %T>%
  figsave('GDTlung_6p.jackstraw.pdf' , w = 400, h = 400)

GDTlung_s  %<>%  RunUMAP( dims = 1:36, 
                   reduction = 'pca', min.dist = 0.15) %>%
  FindNeighbors(dims = c(1:36))

for (i in seq(0.5,1.2,0.1) %>% rev()) {
  GDTlung_s <- FindClusters(GDTlung_s, resolution = i)
}
Feature_rast(GDTlung_s)


allres <- Feature_rast(GDTlung_s, paste0("integrated_snn_res.", seq(0.5,1.2,0.1)), ncol =4) %T>% 
figsave('lunggdt_6p_allres.pdf', 400,200)
 allres




table(GDTlung_s$seurat_clusters)

DefaultAssay(GDTlung_s) <- "RNA"


Feature_rast(GDTlung_s, 'Phase')

Feature_rast(GDTlung_s,c('TRDV2','TRDV1','TRDV3','TRDC','TRGC1','TRGV9', "CD3D",'CD3E'), ncol =4)
ClusterCompare(GDTlung_s,'4','6', group.by = 'integrated_snn_res.0.7' , log2fc = 0.5)
ClusterCompare(GDTlung_s,'1','6',group.by = 'integrated_snn_res.0.9' , log2fc = 0.75)


Idents(GDTlung_s) <- paste0('c',(as.numeric(as.vector(GDTlung_s$integrated_snn_res.0.7))+1)) %>%
  factor(levels = paste0('c',1:12))

Feature_rast(GDTlung_s, c('ident', 'tissue', 'patient')) %T>% figsave('gdt_lung_6p_UMAP_aftertrimming.pdf', 200, 80)


GDTlung_s$Cell_cluster <- Idents(GDTlung_s) 

GDTlung_s$ID <- paste0(GDTlung_s$tissue,'_',GDTlung_s$patient)

Feature_rast(GDTlung_s,c('TRAC','TRBC1','TRBC2'))

 




cl_comp <- GDTlung_s@meta.data %>% group_by(tissue, ID, Cell_cluster) %>%  
  summarise(n = n() ) %>% group_by(ID) %>% 
  mutate(percent = n/sum(n)*100) %>% as.data.frame()
# cl_comp$ID <- factor(cl_comp$ID, levels = c(paste0('luLN#',1:),paste0('lung#',1:3) ))
# cl_comp$group <- factor(cl_comp$group, levels = c(, "NCB"))

cl_comp[cl_comp$tissue =='lung',]$n <-  -cl_comp[cl_comp$tissue =='lung',]$n 
cl_comp[cl_comp$tissue =='lung',]$percent <-  -cl_comp[cl_comp$tissue =='lung',]$percent



cl_comp$percent %>% sum
cl_comp
library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[4:9], brewer.pal(9,'Reds')[4:9])

cl_comp_flow <- ggplot(cl_comp, aes(y = n, x = Cell_cluster, fill = ID, color = ID,
                                     stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = ID_cl#,
                    # labels = c('luLN #1: 1135','#2: 203','#3: 822',
                    #            'lung #1: 1756','#2: 632','#3: 964')
                    )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+
  NULL
cl_comp_flow


cl_comp_flow +theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12))


Umap_GDTlung<- Feature_rast(GDTlung_s,c('Cell_cluster','tissue'), sz = 0.5,
                             noaxis = F, labels = c('AUTO')) %>%
  list(cl_comp_flow) %>% PG(nrow = 1, rw = c(1.7,1), labels = c(NA,'C'))
Umap_GDTlung

genes <- c('TRDV1','TRDV2','TRDV3','TRGV9','TRBC1',
          'ZBTB16','TCF7', 'CCR7','TBX21', 'ITGAE','ITGB2',  
       'GZMA','FCGR3A','IFNG',"XCL1", "AREG",
           'KLRB1', 'KLRC1', 'KLRD1','KLRG1',
          "RORC", 'CCR6','DPP4','CXCR3','CXCR6')

gene_UMAP <- Feature_rast(GDTlung_s,genes,  sz = 0.4)
gene_UMAP


fig1 <- PG(list(Umap_GDTlung, gene_UMAP), nrow =2, rh = c(1,3), labels = c(NA,'D'))
fig1


# figsave(gene_UMAP2, 'test_geneumap_gdt2020_2.pdf', 400, 400)

figsave(fig1, 'fig1_lunggdt_6p_UMAP.pdf', 210,270)

saveRDS(GDTlung_s, 'GDTlung_6p_Seurat_07062021.rds')


# DEGs --------------------------------------------------------------------
hgdTmarkers <- FindAllMarkers(GDTlung_s, test.use = 'bimod', min.pct = 0.10,   only.pos = FALSE )%>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))

hgdTmarkers %>% group_by(cluster) %>%   filter(avg_log2FC >0) %>% summarise(n = n())

top10_hgdT <- hgdTmarkers %>% filter(!grepl('^RP|^MT|^TRG|^TRD|^AC', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC)


Heat_GDTlung <- DoHeatmap(GDTlung_s,features = top10_hgdT$gene,raster = T, 
                          group.colors = umap.colors,size = gs(8))%>%heat_theme() %>%
  figsave('heatmap_GDTlung_6p.pdf',150,250)
Heat_GDTlung


surfacemakers <- c('PTPRC','PDCD1','CCR6','ITGAE','CD69','CCR7','KLRB1','CD27','KLRG1','IL7R',
                   'DPP4','ITGA1','NCR1','NCR3','FCGR3A')


sf <- Feature_rast(GDTlung_s, c('ident','tissue',  surfacemakers), ncol = 4,sz = 0.2) %T>%
figsave('sfmarkers_gdt.pdf', 200,200)
sf
library(Nebulosa)

plot_density(GDTlung_s, surfacemakers)   %T>%
  figsave('density_test.pdf', 400,400)
plot_density(GDTlung_s, 'TRGC2')
Feature_rast(GDTlung_s, 'TRGC2')

Feature_rast(GDTlung_s, c('ident','orig.ident')  )
saveRDS(GDTlung_s,'GDTlung_6p_Seurat_07062021.rds')


Cytokines <- c('IL17A','IL22',  'AREG','IFNG','CSF2', 'TNF',
               'GZMA', 'GZMB','GZMK', 'PRF1')


Cytokines %>% Feature_rast(GDTlung_s, .)
# Gene modules ------------------------------------------------------------
library(openxlsx)


gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_c_list_ent$GM
gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)
gene_cluster$GM_B

for (i in gc_name) {
  GDTlung_s %<>%  AddModuleScore( features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(GDTlung_s@meta.data) %<>% str_replace("(?<=\\w)1", '')
  
}
colnames(GDTlung_s@meta.data)

Feature_rast(GDTlung_s, c(gc_name), color_grd = 'grd') 





ViolinPlot(GDTlung_s, gc_name, colors = umap.colors)

GMS  <- GDTlung_s@meta.data[,c(gc_name)] %>% as.data.frame() %>% t()
scaleGM <- scale(t(GDTlung_s@meta.data[,c(gc_name)]))
scaleGMassay <- CreateAssayObject(data = scaleGM)
GDTlung_s@assays$GM <- scaleGMassay
gms <- paste0('GM-', LETTERS[1:8])

gcanno <- as.vector(c('GM_A: Naive or immature T cell',
                      "GM_B: Innate T cell differentiation",
                      'GM_C: Proliferating',
                      "GM_D: Type-3 immunity",
                      "GM_E: Interferon induced",
                      "GM_F: CTL response (VD2)",
                      "GM_G: CTL response (VD1)",
                      
                      "GM_H: Acute activation"))


c(1,2,4,6,7) %>% map(~
                       Feature_rast(GDTlung_s, gc_name[.x], color_grd = 'grd', mythe = F)+
                       ggtitle(gcanno[.x])
) %>% PG(ncol = 3) %T>% figsave('GM_score.pdf', 200, 150)


GM_score_heat <- DoHeatmap(subset(GDTlung_s, downsample = 1000), 
                           raster =T, draw.lines = T, angle = 45,
                           lines.width = 10,group.colors = umap.colors,
                           
                           assay = 'GM', features = gms, slot = 'data', size = gs(8)) +hmp2 + mytheme+
  theme(legend.position = 'bottom',
        legend.key.height = unit(2,'mm'))+
  guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))
GM_score_heat

# TCR  --------------------------------------------------------------------

# here we have TCR data from two different datasets: 1 and 2_3
# and we will combine them as one
# for each sample  we have two sequence run.

TCRdirs <- list.dirs(path = 'raw') %>% str_subset('VDJ.+outs$' )

read.csv('raw/VDJ_322_325_453_456_gd/outs/filtered_contig_annotations*.csv')




TCRs <- map2(TCRdirs, c(2,1,3), ~ rbind(read.csv(paste0(.x, '/filtered_contig_annotations.csv' )),
                                        read.csv(paste0(.x, '/filtered_contig_annotations_p2.csv' )))%>% 
               filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene)) %>%
               dplyr::select(c(1, 5:10, 13,14))  %>%
               mutate(  
                       bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% reduce(.f = rbind)
TCRs
# 
# lung1TCR <- read.csv('raw/VDJ_Falk1_gd/outs/all_contig_annotations.csv') %>% 
#   filter(productive == 'True' & grepl('GV|DV', v_gene))%>%
#   dplyr::select(c(1, 5:10, 13,14)) %>% mutate(bc_backup = paste0(barcode, "_1"))
# 
# lung2_3TCR <- read.csv('raw/VDJ_322_325_453_456_gd/outs/all_contig_annotations.csv') %>% 
#   filter(productive == 'True' & grepl('GV|DV', v_gene))%>%
#   dplyr::select(c(1, 5:10, 13,14)) %>% mutate(bc_backup = paste0(barcode, "_2"))
# 
# lung4_6TCR <- read.csv('raw/VDJ_falk3_gd/outs/all_contig_annotations.csv') %>% 
#   filter(productive == 'True' & grepl('GV|DV', v_gene))%>%
#   dplyr::select(c(1, 5:10, 13,14)) %>% mutate(bc_backup = paste0(barcode, "_3"))
# 
# 
# lungTCR <- rbind(lung1TCR, lung2_3TCR,lung4_6TCR) %>% select(-barcode)


TRGs <- lungTCR %>% filter(grepl('GV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
  # mutate(v_gene = str_remove(v_gene,'DV\\d')) %>%
  rename_at(vars(-bc_backup), funs(sub('$','_TRG',.)))
TRGs

TRDs <- lungTCR %>% filter(grepl('DV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
  rename_at(vars(-bc_backup), funs(sub('$','_TRD',.))) %>%
  mutate(v_gene_TRD = str_remove(v_gene_TRD, 'AV\\d\\d'))%>% 
  mutate(v_gene_TRD = str_remove(v_gene_TRD, '-2'))
TRDs


TCRs <- full_join(TRGs,TRDs, by = c('bc_backup'))
 TCRs  %<>%  mutate(paired = 
                           case_when(!is.na(v_gene_TRD) & !is.na(v_gene_TRG) ~ 
                                       paste0(str_remove(v_gene_TRG,'TR'),' ',str_remove(v_gene_TRD,'TR'))))
 TCRs  %<>% mutate(cdr3_paired = 
                          case_when(!is.na(v_gene_TRD) & !is.na(v_gene_TRG) ~ 
                                      paste(str_remove(v_gene_TRG,'TR'),cdr3_TRG,str_remove(v_gene_TRD,'TR'),cdr3_TRD))) 

top9paired <- TCRs %>% filter(!is.na(paired))  %>%count(paired)   %>% 
  arrange(desc(n)) %>% top_n(9,n) %>% pull(paired)
unique(TCRs$v_gene_TRD)

TCRs<- TCRs %>% 
  mutate(paired_sp =
           case_when(paired %in% top9paired ~ paired, 
                     !is.na(paired) ~ "Other paired")    ) %>%
  mutate(paired_sp = factor(paired_sp, levels = c(top9paired, "Other paired")))
# paired TRG and J 
TCRs<- TCRs %>% mutate(TRGJP = case_when(!is.na(chain_TRG) ~ paste(v_gene_TRG, j_gene_TRG)))%>%
  mutate(TRGJP = str_replace(TRGJP, ' TRG',' '))

TCRs$bc_backup
TCRs
# # joinTCR to seurat --------------------------------------------------------
# GDTlung is the name of seurat object


GDTlung_s$bc_backup <- rownames(GDTlung_s@meta.data)

metabackup <- GDTlung_s@meta.data


GDTlung_s@meta.data  %<>% left_join(TCRs, by = 'bc_backup') %>% `rownames<-`(GDTlung_s@meta.data$bc_backup )

GDTlung_s@meta.data

# this line to remove repeated join
GDTlung_s@meta.data %<>% `colnames<-`(str_remove(colnames(GDTlung_s@meta.data), '.y')) %>%
  select(!ends_with('.x')) 
# mapped TCR
TRDmapp <- GDTlung_s@meta.data %>% count(ID, v_gene_TRD) %>%  group_by(ID)  %>%
     mutate(percent = n/sum(n)*100) 
TRGmapp <- GDTlung_s@meta.data %>% count(ID,v_gene_TRG) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 


pairedmapp <-GDTlung_s@meta.data %>% count(ID,paired_sp) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 



GDTlung_s@meta.data %<>% mutate(TCR_summary = case_when(!is.na(paired) ~ 'paired TCR',
                                                     !is.na(chain_TRG)~'single TRG',
                                                     !is.na(chain_TRD)~ 'single TRD') )


mappedTCR <- map2(list(TRDmapp, TRGmapp, pairedmapp), list('v_gene_TRD','v_gene_TRG', 'paired_sp'),~
                    ggplot(.x, aes_string('ID', y = "percent", color = .y, fill = .y))+
  geom_bar(stat = 'identity',  width = 0.5,position = position_stack(reverse = T))+
    color_m()+
    fill_m()+
    xlab(NULL)+
    ylab('%')+
    ggtitle(.y)+
    guides(color = F, fill = guide_legend(nrow =4, title = NULL))+
    theme_minimal()+
    mytheme+gglp('b')
  
) %>% PG(nrow =1, labels = 'AUTO')
mappedTCR
ViolinPlot(subset(GDTlung, length_TRD >0 & patient == '#3'),'length_TRD')


TCRmapping <- GDTlung_s@meta.data %>% group_by(tissue, patient) %>%  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)
(ggplot(TCRmapping,aes(x = tissue,y= percent, group = TCR_summary, fill = TCR_summary))+geom_bar(stat = 'identity',position = position_stack(reverse = T))+facet_grid(~patient)+fill_m()+theme_bw()+
    scale_y_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))+
    mytheme) %T>% figsave('mappingrate_TCR.pdf', 150, 60)


figsave(mappedTCR,'lunggdt_mappedTCR.pdf', 200,100)


# TCR frequency -----------------------------------------------------------


cdr3TRD_freq <- GDTlung_s@meta.data %>% filter(c_gene_TRD == 'TRDC' & cdr3_TRD != 'None') %>%  group_by(patient) %>%
  dplyr::count(cdr3_TRD = cdr3_TRD) %>% arrange(desc(n)) %>% dplyr::rename(cdr3_TRD_freq = 'n') %>%
  mutate(cdr3_TRD_perc = cdr3_TRD_freq/sum(cdr3_TRD_freq)*100)

cdr3TRG_freq <- GDTlung_s@meta.data %>% filter(grepl('TRGC', c_gene_TRG) & cdr3_TRG != 'None') %>%  group_by(patient) %>%
  dplyr::count(cdr3_TRG = cdr3_TRG) %>% arrange(desc(n)) %>% dplyr::rename(cdr3_TRG_freq = 'n') %>%
  mutate(cdr3_TRG_perc = cdr3_TRG_freq/sum(cdr3_TRG_freq)*100)
cdr3TRG_freq

cdr3TRG_freq$cdr3_TRG_perc %>% sum
cdr3Paired_freq <- GDTlung_s@meta.data  %>% filter(!is.na(cdr3_paired)) %>% group_by(patient) %>%
  dplyr::count(cdr3_paired = cdr3_paired) %>% arrange(desc(n)) %>% 
  dplyr::rename(cdr3_paired_freq = 'n') %>%
  mutate(cdr3_paired_perc = cdr3_paired_freq/sum(cdr3_paired_freq)*100)

cdr3Paired_freq


GDTlung_s@meta.data %<>% left_join(cdr3TRD_freq, by = c('cdr3_TRD','patient')) %>% 
  left_join(cdr3Paired_freq, by = c('cdr3_paired','patient'))%>% 
  left_join(cdr3TRG_freq, by = c('cdr3_TRG','patient')) %>% `rownames<-`(GDTlung_s$bc_backup)

# rownames(GDTlung@meta.data) <-
##paired_simplified


saveRDS(GDTlung_s,'GDTlung_6p_Seurat_07062021.rds')

patient <- c(GDTlung_s$patient %>% unique()) %>% sort()

TCRD_FREQ <-map(patient,~
      Feature_rast(subset(GDTlung_s, patient == .x),'cdr3_TRD_perc',color_grd = 'grd')+
        ggtitle(paste('patient',.x,'\nTRD_clone (%)'))
) %>% PG(nrow =2)

TCRD_FREQ


Feature_rast(GDTlung_s, 'cdr3_TRD_perc', color_grd = 'grd', mythe = F)+
  ggtitle('TRD_clone (%)')+
  theme(plot.title = element_text(size = 12))




# analysis and make figures -----------------------------------------------



# TRD on umap
TCRumap <- Feature_rast(GDTlung_s,c('Cell_cluster','v_gene_TRD','paired_sp'),do.label = F, sz =1, labels =LETTERS[4:6]   )
TCRumap
Feature_rast(GDTlung_s,'ident', c('ID'))+facet_grid(~ID)


GDTlung_s$length_TRD

ViolinPlot(GDTlung_s %>% subset(v_gene_TRD == 'TRDV1') , 'length_TRD', box = T)+ ylim(350,700)


Feature_rast(GDTlung, facet = 'patient')+facet_grid(~patient)

# top TRD and paired

Top5TRD_lung <- GDTlung_s@meta.data %>% filter( !is.na(cdr3_TRD) & tissue == 'lung') %>%
  group_by(patient,v_gene_TRD, cdr3_TRD) %>%
  dplyr::count(patient,cdr3_TRD) %>% arrange(desc(n)) %>%
   group_by(patient) %>%  dplyr::slice(1:5) %>% filter(n>=5) %>%
  dplyr::mutate(Top5TRD_lung = paste0(patient, ' ',v_gene_TRD," ",cdr3_TRD,' ',  n)) %>% 
  arrange(patient,desc(n)) %>%
  select(patient, cdr3_TRD, Top5TRD_lung) %>%
   ungroup()




Top5TRD_luLN <- GDTlung_s@meta.data %>% filter( !is.na(cdr3_TRD) & tissue == 'luLN') %>%
  group_by(patient,v_gene_TRD, cdr3_TRD) %>%
  dplyr::count(patient,cdr3_TRD) %>% arrange(desc(n)) %>%
 group_by(patient) %>%  dplyr::slice(1:5) %>% filter(n>=5) %>%
  dplyr::mutate(Top5TRD_luLN = paste0(patient, ' ',v_gene_TRD," ",cdr3_TRD,' ',  n)) %>% 
  arrange(patient,desc(n)) %>%
  select(patient, cdr3_TRD, Top5TRD_luLN) %>%
  ungroup()


# top10TRD <- GDTlung@meta.data %>% filter( !is.na(cdr3_TRD)) %>% 
#   group_by(tissue,v_gene_TRD, cdr3_TRD) %>%dplyr::count(tissue,patient,cdr3_TRD) %>%
#   arrange(desc(n)) %>% group_by(tissue,patient) %>%  dplyr::slice(1:10) %>%
#   dplyr::mutate(top10_TRD = paste0(tissue, " ",v_gene_TRD," ",cdr3_TRD)) %>%
#   arrange(tissue,desc(n)) %>% ungroup()

# top10TRD %>% group_split(patient) %>% map(.,~ .x%>% pull(cdr3_TRD) %>% unique() %>% length)
# 
# 
# 
# top10TRD %>% filter(patient == '#3')

GDTlung_s$Top5TRD_lung <- NULL
GDTlung_s$Top5TRD_luLN <- NULL
 GDTlung_s@meta.data %<>% 
  left_join(Top5TRD_lung,by = c('cdr3_TRD','patient') )%>%
   left_join(Top5TRD_luLN,by = c('cdr3_TRD','patient') )%>%
  mutate(Top5TRD_lung = case_when( !is.na(Top5TRD_lung) ~ Top5TRD_lung,
                                  is.na(Top5TRD_lung) & !is.na(cdr3_TRD)  ~ 'mapped TCR')) %>%
   mutate(Top5TRD_luLN = case_when( !is.na(Top5TRD_luLN) ~ Top5TRD_luLN,
                                    is.na(Top5TRD_luLN) & !is.na(cdr3_TRD)  ~ 'mapped TCR')) %>%
  `rownames<-`(GDTlung_s$bc_backup)
 GDTlung_s$Top5TRD_lung %<>% factor(levels = GDTlung_s$Top5TRD_lung %>%unique() %>%  str_sort())
 GDTlung_s$Top5TRD_luLN %<>% factor(levels = GDTlung_s$Top5TRD_luLN %>%unique() %>%  str_sort())
 

 Feature_rast(GDTlung_s)

 # GDTlung@meta.data %<>% `colnames<-`(str_remove(colnames(GDTlung@meta.data), '.y')) %>% select(!ends_with('.x')) 
 
 Feature_rast(GDTlung_s, c('tissue', 'Top5TRD_lung','Top5TRD_luLN'), do.label = F, 
              colorset = ggplotColours(28))
 
 top10trdcl <-
   c('#993333','#ff3300','#ff0066','#cc0099','#cc3300','#391313','#ff6666','#ffa64d','#ff8080','#660066',
     '#000066','#0000ff', '#336699','#006666','#6666ff', '#004d3d','#00b3b3','#003399','#39ac73','#1ab2ff',
     '#ffffcc'
   )
 
 luLN_expanded<-Feature_rast(GDTlung_s, "Top5TRD_luLN", do.label = F, facet ='patient')+
     scale_color_manual(values = c(  
       top10trdcl[1:(length(unique(GDTlung_s$Top5TRD_luLN))-2)],
                          'mapped TCR' = alpha('yellow',0.5)), 
       na.value =alpha('lightgrey',0.5) )+facet_grid(~patient)+
   ggtitle("Most Expanded TRD in luLN")
 
 
 lung_expanded<- Feature_rast(GDTlung, "Top5TRD_lung", do.label = F, facet ='patient')+
   scale_color_manual(values = c(  
     top10trdcl[1:(length(unique(GDTlung$Top5TRD_lung))-2)],
     'mapped TCR' = alpha('yellow',0.5)), 
     na.value =alpha('lightgrey',0.5) )+facet_grid(~patient)+
   ggtitle("Most Expanded TRD in lung")
 
 

# venny nonVD2  TRD  ----------------------------------------------------


Feature_rast(GDTlung_s)
GDTlung_s@meta.data %<>% mutate(big_cluster = case_when(Cell_cluster %in% c('c1','c2','c8', 'c9')~ 'nonVD2_lung_1',
                                                      Cell_cluster %in% c('c5','c7' ,'12')~ 'nonVD2_lung_2',
                                                      Cell_cluster %in% c('c3', 'c6', 'c10')~ 'nonVD2_luLN_1',
                              Cell_cluster %in% c('c11')~ 'nonVD2_luLN_2',
                              Cell_cluster %in% c('c4')~ 'VD2')
)%>%
  `rownames<-`(GDTlung_s$bc_backup)
GDTlung_s$big_cluster

UMAP_bigcluster <- Feature_rast(GDTlung_s, 'big_cluster')

ClusterCompare(GDTlung_s,'nonVD2_lung_1','nonVD2_lung_2', group.by = 'big_cluster', log2fc = 1)

patient

TRD_uniq <- c()
  for (i in patient) {
    TRD_uniq[[i]] <-    map(c('nonVD2_lung_1','nonVD2_lung_2','nonVD2_luLN_1','nonVD2_luLN_2'), ~ 
                              GDTlung_s@meta.data %>% 
          dplyr::filter( big_cluster==.x & v_gene_TRD != 'TRDV2' & !is.na(cdr3_TRD) &patient ==i) %>% 
          pull(cdr3_TRD) %>% unique ) %>% 
      setNames(c('nonVD2_lung_1','nonVD2_lung_2','nonVD2_luLN_1','nonVD2_luLN_2'))
  }


TRD_uniq$p25


library(eulerr)
library(ggplotify)

set.seed(8)
euler(TRD_uniq$`#3`, shape = 'ellipse')

set.seed(8)
Venny_nonVD2 <- map2(TRD_uniq, list(F,F,list(cex = .66),F,F,list(cex = .66)),~ plot(euler(.x, shape = 'ellipse'), 
                        quantities = list(cex = 0.66),  prob = T,
                        labels = NULL,legend =.y,
                        fill = alpha(umap.colors,0.8), 
) %>% ggplotify::as.ggplot() + ggtitle('nonVD2 unique CDR3')+ 
  theme(plot.title = element_text(size = 8, hjust = 0.5))) %>% 
   PG(nrow = 2, rw = c(1,1,1.8), labels = patient, label_fontface = 'plain',label_size = 8 )

# venny VD2  TRD  ----------------------------------------------------
TRD_VD2_uniq <- c()
for (i in patient) {
  TRD_VD2_uniq[[i]] <-    map(c('lung','luLN'), function(x) GDTlung_s@meta.data %>% 
                            dplyr::filter( Cell_cluster== 'c4' & v_gene_TRD == 'TRDV2' & patient ==i & tissue == x) %>% 
                            pull(cdr3_TRD) %>% unique ) %>% 
    setNames(c('lung','luLN'))
}

set.seed(8)
Venny_VD2 <- map2(TRD_VD2_uniq, list(F,F,list(cex = .66),F,F,list(cex = .66)),~ plot(euler(.x, shape = 'ellipse'), 
                                                                quantities = list(cex = 0.66),  prob = T,
                                                                labels = NULL,legend =.y,
                                                                fill = alpha(umap.colors,0.8), 
) %>% ggplotify::as.ggplot() + ggtitle('VD2 unique CDR3 in c5')+ 
  theme(plot.title = element_text(size = 8, hjust = 0.5))) %>% 
  PG(nrow = 2, rw = c(1,1,1.8), labels = patient, label_fontface = 'plain',label_size = 8 )


# TCR_comp <- PG(list(Venn_TRD_uniq, Venn_paired_uniq), labels = c('h','i'))
Fig2test <- PG(list(mappedTCR,TCRumap,TCRD_FREQ), 
               ncol = 1)
figsave(Fig2test,'fig2_lunggdt_TCRcomp.pdf', 195, 230)


UMAP_tissue <- Feature_rast(GDTlung, c('tissue', 'big_cluster'))

fig3_TCRsharing <- PG(list(   PG(list(UMAP_tissue, NA), rw = c(2,1)),
                         lung_expanded, luLN_expanded, 
                           Venny_nonVD2  ), ncol = 1, labels = 'AUTO', rh = c(1.2,1.5,1.5,1.2))

fig3_TCRsharing

figsave(fig3_TCRsharing,'fig3_lunggdt_TCRsharing.pdf', 195, 270)
saveRDS(GDTlung,GDTlung.rds)


# # TCR sharing by barplot ------------------------------------------------

library(ggalluvial)
# non VD2 

GDTlung_s@meta.data %>% ungroup()%>% filter(cdr3_TRD_freq >1 ) %>% top_n(1000, cdr3_TRD_freq) %>% pull(cdr3_TRG) %>% unique()

GDTlung_s@meta.data %>% ungroup()%>% filter(cdr3_TRD_freq >0 ) %>% top_n(100, cdr3_TRD_freq)

TCRfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_paired_freq >1 )  %>% group_by(Cell_cluster, cdr3_paired, patient) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq)

TCRfreqtable %>% filter(Cell_cluster %in% c('c1', 'c2' ,'c5', 'c7', 'c6', 'c3')) %>%  
  mutate(Cell_cluster = factor(Cell_cluster, levels = c('c1', 'c2' ,'c3','c6', 'c5', 'c7' ))) %>% 


ggplot(
       aes( x = Cell_cluster, y = TCRfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  ggtitle("TCRD sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  scale_fill_manual(values = rainbow(514))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+facet_wrap(~patient)


# VD2

VD2freqtable <- GDTlung_s@meta.data %>% filter(cdr3_TRD_freq >0 &v_gene_TRD == 'TRDV2' & Cell_cluster == 'c4')  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq)
VD2freqtable



ggplot(
  aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("VD2 TCRD sharing in c4")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  scale_fill_manual(values = rainbow(200))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq
# CITESEQ analysis --------------------------------------------------------




GDTlung_s@assays$CITE

GDTlung_s$orig.ident %>% unique



GDTlung_cite <- subset(GDTlung_s, orig.ident == "p71&p31&p27")


Feature_rast(GDTlung_cite)
DefaultAssay(GDTlung_cite) <- 'CITE'

# antibodies <- GDTlung_cite@assays$CITE@counts 
rownames(GDTlung_cite@assays$CITE@counts ) 
abnames <- paste0(c('CD4', 'CD8','CD45RA', 'CD45RO', 'PD1', 'CXCR3', 'CCR6','CD103', 'CD69', 'CCR7',
                   'KLRB1', 'CD27', 'KLRG1', 'IL7R', 'CD26', 'CD49a', 'KLRD1'), '.protein')


rownames(GDTlung_cite@assays$CITE@counts ) <- abnames
rownames(GDTlung_cite@assays$CITE@data ) <- abnames

rownames(GDTlung_s@assays$CITE@counts ) <- abnames
rownames(GDTlung_s@assays$CITE@data ) <- abnames

saveRDS(GDTlung_s, 'GDTlung_6p_Seurat_07062021.rds')
saveRDS(GDTlung_cite, 'GDTlung_CITEseqvisual.rds')


GDTlung_cite %<>% NormalizeData(normalization.method ='CLR', margin = 2)
GDTlung_s %<>% NormalizeData(normalization.method ='CLR', margin = 2)

GDTlung_cite

CITEprotein <- Feature_rast(GDTlung_cite, abnames,color_grd = 'grd', sz = 0.4, ncol =6) %T>% figsave('CITEseq.pdf', 200, 120)


genenames <- c('ident', 'CD4', 'CD8A','PTPRC', 'PDCD1', 'CXCR3', 'CCR6','ITGAE', 'CD69', 'CCR7',
               'KLRB1', 'CD27', 'KLRG1', 'IL7R', 'DPP4', 'ITGA1', 'KLRD1')
DefaultAssay(GDTlung_cite) <- 'RNA'

CITEgene <- Feature_rast(GDTlung_cite, genenames, labelsize = 6, sz = 0.4, ncol =6) 


CITEgene

PG(list(CITEgene, CITEprotein), nrow = 2 ) %T>% figsave('genevsprotein.pdf', 270, 200)



Feature_rast(GDTlung_s, 'MKI67')




GDTlung_s@meta.data  %<>%  mutate(Vg9Vd2 = case_when(paired_sp == 'GV9 DV2' ~ 'GV9 DV2',
                                                     !is.na(paired_sp)  ~ 'other γδTCR'
                                                     )) %>% `row.names<-`(GDTlung_s$bc_backup)



saveRDS(GDTlung_s, 'GDTlung_6p_Seurat_07062021.rds')



# integration with gdT from adult PBMC ------------------------------------

Adult_GDT_2020AUG <- readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds') %>% 
  subset(group == 'Adult') 
Adult_GDT_2020AUG%<>%   ScaleData(assay = 'RNA', 
                                      vars.to.regress = c("percent.mito", 'orig.ident',
                                                          "S.Score", 'G2M.Score',
                                                          "percent.ribo",
                                                          'nCount_RNA','nFeature_RNA' ) )

DefaultAssay(Adult_GDT_2020AUG) <- 'RNA'

GDTlung_s$tissue

# annotaion pehotypes 
Adult_GDT_2020AUG@meta.data %<>% mutate(
  project = 'PBMC',
  tissue = 'blood',
  pheno = case_when( 
  Cell_cluster %in% c('c1', 'c2', 'c3', 'c4' ) ~ 'Naive',
  Cell_cluster == 'c5' ~ 'Vg9Vd2 type2',
  Cell_cluster == 'c6' ~ 'Vg9Vd2 type3',
  Cell_cluster == 'c7' ~ 'Vg9Vd2 type1 Tmem',
  Cell_cluster == 'c8' ~ 'Vg9Vd2 type1 Teff',
  Cell_cluster == 'c9' ~ 'Vg9Vd2 IFN induced',
  Cell_cluster == 'c10' ~ 'Vd1 Tem',
  Cell_cluster == 'c11' ~ 'Vd1 Teff',
  Cell_cluster == 'c12' ~ 'early activation')
           ) %>% `rownames<-`(Adult_GDT_2020AUG$bc_backup)

Feature_rast(Adult_GDT_2020AUG, 'pheno')

Adult_GDT_2020AUG[['integrated']] <- NULL
Adult_GDT_2020AUG[['GM']] <- NULL

GDTlung_inte <- GDTlung_s

GDTlung_inte  %<>% AddMetaData(GDTlung_inte@assays$CITE@data %>% t())
GDTlung_inte@meta.data

GDTlung_inte[['integrated']]<- NULL
GDTlung_inte[['HTO']]<- NULL
GDTlung_inte[['CITE']]<- NULL
GDTlung_inte[['GM']]<- NULL
# integration
# GDT_lung_PB$tissue %>% unique



# anchors <- FindIntegrationAnchors(list(GDTlung_inte, Adult_GDT_2020AUG), dims = 1:60)

# anchors
# GDT_lung_PB <- IntegrateData(anchors, dims = 1:20)
# 
# 
# DefaultAssay(GDT_lung_PB) <- "integrated"
# 
# 
# GDT_lung_PB %<>% ScaleData(vars.to.regress = c("percent.mito", 
#                                              # 'patient',
#                                              "S.Score", 'G2M.Score',
#                                              "percent.ribo",
#                                              'nCount_RNA','nFeature_RNA' )) %>% 
#   RunPCA( npcs = 100, verbose =   T, nfeatures.print = 40) %>% 
#   JackStraw( num.replicate = 100, dims = 80)%>%
#   ScoreJackStraw(dims = 1:80) 
# JackStrawPlot(GDTlung_s, dims = 1:50 ) %T>%
#   figsave('GDT_lung_PB.jackstraw.pdf' , w = 400, h = 400)

# integration by harmony --------------------------------------------------
library(harmony)

GDT_lung_PB<- merge(x = GDTlung_inte, y = Adult_GDT_2020AUG,add.cell.ids = c("lung", "bl" ), project = 'lb' )
# GDT_lung_PB<- merge(x = Adult_GDT_2020AUG, y = GDT_2020AUG)
dim(GDT_lung_PB)
  

GDT_lung_PB%<>% NormalizeData(verbose = FALSE, assay = 'RNA') %>% 
  FindVariableFeatures(selection.method = 'vst', nfeatures = 3500) 


GDT_lung_PB@assays$RNA@var.features %<>%str_subset('^TRG|^TRD|^TRB|^TRA|^MT|^IG|^MIR|^HSP', negate = T) #%T>% print(length())
GDT_lung_PB@assays$RNA@var.features

#%>% 
GDT_lung_PB%<>%   ScaleData(verbose = FALSE, assay = 'RNA',    
                        vars.to.regress = c('percent.mito',  'nCount_RNA','nFeature_RNA' )
) %>% 
  RunPCA(npcs = 100, verbose=T)
GDT_lung_PB%<>% RunHarmony(c(  'patient', 'percent.mito',  'nCount_RNA','nFeature_RNA'), plot_convergence = TRUE)



ElbowPlot(GDT_lung_PB, ndims = 50,reduction = 'harmony')

saveRDS(GDT_lung_PB ,'GDT_lung_PB_harmony.rds')

#chose 35 PCs 





# dim reduce via harmony --------------------------------------------------------------




GDT_lung_PB %<>% RunUMAP(reduction = 'harmony', dims = 1:30) %>% FindNeighbors(redcution = 'harmony', dims = 1:30)

for (i in seq(0.6,1.5,0.1) %>% rev()) {
  GDTlung <- FindClusters(GDTlung, resolution = i)
}

GDTlung %<>% CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


Feature_rast(GDTlung, c('ident', 'tissue', 'patient'))

ViolinPlot(GDTlung, c('nCount_RNA','percent.mito', 'percent.ribo'), colors = umap.colors)

ClusterCompare(GDTlung, '3', '7')


# Since by harmoony  we failed to achieve clear UMAP, we decide to stay with Seurat protocol
saveRDS(GDTlung, 'GDTlung.harmony:integration.rds')#


# Figs --------------------------------------------------------------------




Feature_rast(GDTlung_s, 'paired_sp', colorset = c(umap.colors[1:9], 'lightblue'), sz = 1.2)


umap_idnet_lunggdt <- Feature_rast(GDTlung_s, c('ident'), noaxis = F, labelsize = 12)+ggtitle('Phenotypes of lung γδT cells')+
  theme(plot.title = element_text(size = 12, face = 'plain'))



vg9d2 <- Feature_rast(GDTlung_s, 'Vg9Vd2', colorset = c('red', 'yellow'),labelsize = 12)+ggtitle('γδTCR')+
  theme(plot.title = element_text(size = 12, face = 'plain'))


GMD_vio <- ViolinPlot(GDTlung_s, 'GM_D', colors = umap.colors, box = T)+ggtitle('Gene module: type-3 immunity')+
  theme(plot.title = element_text(size = 12, face = 'plain'))


PG(list(umap_idnet_lunggdt, vg9d2, GMD_vio), ncol = 3) %T>% figsave('Vg9Vd2areType3', 200, 80) 

GDT_2020AUG <- readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds')


Feature_rast(GDT_2020AUG, 'ZBTB16')+ggtitle('ZBTB16 (PLZF)')+
  theme(plot.title = element_text(size = 16))



# VIPER and PISCES --------------------------------------------------------
library(PISCES)

sDefaultAssay(GDTlung_s) <- 'integrated'
 
GDTlung_s@assays$integrated@misc


GDTlung_s <- CorDist(GDTlung_s)

sil.list <- list()
for (i in seq(0.5, 1, by = 0.1)) {
  clust.name <- paste('integrated_snn_res.', i, sep = '')
  clust.vect <- as.numeric(GDTlung_s[[clust.name]][,1])
  sil.list[[clust.name]] <- mean(cluster::silhouette(clust.vect, GDTlung_s@assays$integrated@misc$dist.mat)[,3])
}

sil.list$integrated_snn_res.0.7
meta.mats <- MetaCells(as.matrix(GDTlung_s@assays$RNA@counts), GDTlung_s@assays$integrated@misc$dist.mat, 
                       GDTlung_s$Cell_cluster, min.samps = 300, num.neighbors = 10)

for (m.name in names(meta.mats)) {
  f.name <- paste('viper/GDTlung', m.name, '-meta.rds', sep = '')
  meta.mats[[m.name]] <-  as.data.frame(meta.mats[[m.name]]) %>% tibble::rownames_to_column('gene') %>% as.data.frame()
  write.table(meta.mats[[m.name]],paste0('viper/GDTlung', m.name, '-meta.txt'), sep = '\t',quote = FALSE, col.names = T, row.names = F)
  # saveRDS(meta.mats[[m.name]], f.name)
}

read.table('viper/GDTlungc1-meta.txt', header = T)

GDTlung_s  %<>% AddPISCESAssay()  %>% CPMTransform() %>% GESTransform()



Feature_rast(GDTlung_s, c('ident', 'paired_sp'))


cbmc_networks %>% glimpse()
cbmc_networks[[1]] %>%  head()


GDTlung_s %<>% PISCESViper( cbmc_networks)
GDTlung_s %<>% CorDist() %>% LouvainResRange(rmin =10, rmax =300) %>% MakeUMAP() %T>% PlotClusters()


PlotClusters(GDTlung_s)
Feature_rast(GDTlung_s)


GDTlung_s@assays$PISCES@data


Feature_rast(GDTlung_s, 'CD3D', color_grd = 'grd')


DefaultAssay(GDTlung_s) <- 'RNA'


c('CD3D','TBX21','EOMES','RORC','KLRG1','CD27','CCR7','CCR6') %>% 
  DoHeatmap(GDTlung_s,.) %>% heat_theme()


GDTlung_s@assays$PISCES <- NULL

GDTlung_s$patient


Feature_rast(GDTlung_s, 'tissue',  other = c('cdr3_TRD', 'patient'))+ geom_line(
  data = FetchData(subset(GDTlung_s, cdr3_TRD_freq >2), c('UMAP_1', 'UMAP_2', 'cdr3_TRD','patient')),
  aes(group = cdr3_TRD),color = alpha('black', 0.3))+facet_wrap(~patient)





Feature_rast(GDTlung_s, 'tissue', facets = 'patient')


GDTlung_s@meta.data$cdr3_paired

GDTlung_s@meta.data %>% top_n(1000,cdr3_TRD_freq) %>% pull(cdr3_TRG) %>% unique() %>% length()

GDTlung_s@meta.data %>% top_n(1000,cdr3_TRD_freq) %>% filter(!is.na(cdr3_paired)) %>% arrange(cdr3_TRD) %>% 
  select(cdr3_TRD, cdr3_TRD_freq, cdr3_TRG, cdr3_TRG_freq) %>% distinct()

