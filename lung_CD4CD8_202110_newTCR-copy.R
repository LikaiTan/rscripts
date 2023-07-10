
# lung CD4 and CD8 project
library(Seurat)
library(ggplot2)
library(purrr)
library(dplyr)
library(Matrix)
library(tibble)
library(cowplot)
library(ggrastr)
library(Cairo)

library(ggrepel)
library(openxlsx)
library(stringr)
library(magrittr)
# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')

# file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

multicores()


#paralell computation

# work dir ----------------------------------------------------------thx

setwd('/home/big/tanlikai/Lung/abt')

CD4CD8 <- readRDS('CD4CD8_integrated_2022_8p.rds')

Feature_rast(CD4CD8)
# read_rawdata ---------------------------------------------------------------

# directories to raw gene expression and CITEseq data

dirs <- list.dirs(path = 'raw') %>% str_subset('surface.+outs$' )

dirs

proj <- c('p32', 'p45','p25_CD4', 'p25_CD8','p27','p31','p71')
pt <- c('p32', 'p45','p25', 'p25','p27','p31','p71')

# read raw GEX and CITEseq data 
rawdata <- dirs %>% map(~ Read10X_h5(paste0(.x, '/filtered_feature_bc_matrix.h5')) ) %>% setNames(proj)
# The quality of p31 is very poor, so we need remove it.
rawdata$p31 <- NULL
proj <- c('p32', 'p45','p25_CD4', 'p25_CD8','p27','p71')
pt <- str_extract(proj, 'p\\d\\d')
# creat seurat project of gene expression datasepeately into a big list

lungTcell <- map2(rawdata,proj, ~
                   CreateSeuratObject(counts = .x$`Gene Expression`,project = .y,
                                       min.features = 20) %>% 
                    AddMetaData(str_extract(.y, 'p\\d\\d'),col.name = 'patient')
                 )




map(lungTcell, ~ .@meta.data %>% colnames )
map(lungTcell, ~ .@meta.data %>% ncol )

map(lungTcell, ~ dim(.))

# the data quality of p31 is very poor, thus need to be removed 
lungTcell$p31 <- NULL

# demultiplexing ----------------------------------------------------------



#  citeseq and hashtaq data

hash <- unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )) %>%  str_subset('^\\d_')
hash

cite <- setdiff(unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )), hash)
cite

# integrate HTO and CITESEQdata into seurat projects 

for (x in proj) {
  lungTcell[[x]][['HTO']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(hash, rownames(rawdata[[x]]$`Antibody Capture`)),])
  lungTcell[[x]][['CITE']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(cite, rownames(rawdata[[x]]$`Antibody Capture`)),])
  rownames(lungTcell[[x]][['HTO']]@counts) <- paste0('X',rownames(lungTcell[[x]][['HTO']]@counts) )
  rownames(lungTcell[[x]][['HTO']]@data) <- paste0('X',rownames(lungTcell[[x]][['HTO']]@data) )
  
  
}

# change names of HTO to more recognizable names 

for (x in proj) {
  rownames(lungTcell[[x]][['HTO']]@counts) %<>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+') %>% paste0('S',.)
  rownames(lungTcell[[x]][['HTO']]@data)%<>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+')%>% paste0('S',.)
  
  
}

hash_sample <- lungTcell %>% map(~ .[['HTO']] %>% rownames)
hash_sample
lungTcell$p32@assays$HTO@counts

# normalization
for (i in proj) {
  lungTcell[[i]] %<>%  NormalizeData( assay = "HTO", normalization.method = "CLR") %>% 
    ScaleData(assay='HTO', features = rownames( .[['HTO']])) %>% 
    HTODemux( assay = "HTO", positive.quantile = 0.9, nstarts = 20)
}


for (i in proj) {
  lungTcell[[i]] %<>%  NormalizeData( assay = "CITE", normalization.method = "CLR") 
}


map(lungTcell, ~ .@meta.data %>% colnames )
map(lungTcell, ~ .@meta.data %>% ncol )

for (x in proj) {
  DefaultAssay(lungTcell[[x]] ) <- 'HTO'
}

saveRDS(lungTcell,'lungTcell_pre_inte.rds')

lungTcell <- readRDS('lungTcell_pre_inte.rds')



lungTcell$p32$hash.ID
map2(lungTcell,hash_sample,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]], assay = 'HTO',
                                             g = 'hash.ID',noaxis =F, axis.number=T)+
       xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.5)+
       geom_vline(xintercept = 0.8)) %>% PG(ncol =6) %T>%
  figsave('hash.id.automatic.pdf', 400, 100)
Feature_rast(lungTcell$p32, 'hash.ID')

# the automiatic demultiplexing by seurat is not very accurate 

# add HTO to metadata

for (x in proj) {
 lungTcell[[x]]  %<>% AddMetaData(t(lungTcell[[x]]@assays$HTO@data), col.name = rownames(lungTcell[[x]]@assays$HTO@data))
}

for (x in proj) {
  lungTcell[[x]]  %<>% AddMetaData( rownames(lungTcell[[x]]@meta.data), 
                                    col.name = 'bc_backup')
}


#  demultiplexing manually

map(lungTcell,~ ncol(.x@meta.data)
    )

syms(hash_sample)
newmeta <- list()

# marker tissues based on HTO

newmeta <- map(lungTcell,  ~
                .x@meta.data  %<>%     
                  mutate(tissue = case_when(.[[15]] >= 0.8 & .[[16]] < 0.5 ~ 'lung',
                                            .[[15]]  < 0.8 & .[[16]] >= 0.5 ~ 'luLN',
                                            .[[15]]  >= 0.8 & .[[16]] >= 0.5 ~ 'DP'))%>%
  `rownames<-`(.x$bc_backup)
)
                
newmeta$p32  %<>%    
  mutate(tissue = case_when(.[[15]] >= 0.8 & .[[16]] < 0.5 ~ 'luLN',
                            .[[15]]  < 0.8 & .[[16]] >= 0.5 ~ 'lung',
                            .[[15]]  >= 0.8 & .[[16]] >= 0.5 ~ 'DP'))%>%
  `rownames<-`(.$bc_backup)
  
for (i in proj) {
  lungTcell[[i]]@meta.data <- newmeta[[i]]
}

map2(lungTcell,hash_sample,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]],
                                             g = 'tissue',noaxis =F, axis.number=T)+
       ggtitle(unique(.x$orig.ident))+
       xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.5)+
       geom_vline(xintercept = 0.8)) %>% PG(ncol =6) %T>%
  figsave('hash.id.mannual.pdf', 400, 100)


for (x in proj) {
  DefaultAssay(lungTcell[[x]] ) <- 'RNA'
}




# QC ----------------------------------------------------------------------

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
lungTcell$p32$percent.ribo

# calculate percent of mitochondrial  ribosam  and hspa genes

lungTcell <- map(lungTcell, ~ PercentageFeatureSet(.x, '^MT', col.name =  'percent.mito') %>% 
             PercentageFeatureSet('^RP', col.name = 'percent.ribo')    %>% 
  PercentageFeatureSet('^HSPA', col.name =  'percent.hspa') )
  
# QC plot

QCvio <- map(lungTcell, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                colors =umap.colors  ,box = T, jitter = T , ncol = 3  )) %>% 
    PG(labels = proj, ncol = 1) %T>%
  figsave('beforeQC_violin.pdf', 150, 270)

  QCvio

  figsave(QCvio,'beforeQC_violin_CD4CD8.pdf', 200, 200)
  

QC_scatter <- map(lungTcell, ~ Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', 
                                      noaxis = F, axis.number = T)+grd+
                    geom_smooth(method = "lm")+
                    
                    scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
                    scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
                    geom_hline(yintercept = c(400,2000))+geom_vline(xintercept = c(1200,7000)))  %>% 
  PG(labels = proj)%T>%
  figsave('beforeQC_scatter.pdf',600,400)
QC_scatter

figsave(QC_scatter,'beforeQC_scatter.pdf',600,300)


saveRDS(lungTcell, 'lungTcell_pre_inte.rds')
# data cleaning 
for (x in proj) {
  lungTcell[[x]] <- subset(lungTcell[[x]], subset =  nCount_RNA %in% 1200:7000 &
                       nFeature_RNA %in% 400:2000 &  percent.mito <15 &
                       tissue %in% c('lung','luLN'))
}
QCvio_clean <- map(lungTcell, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                              colors =umap.colors  ,box = T, jitter = T , ncol = 3  )) %>% 
  PG(labels = proj, ncol = 2) %T>%
  figsave('afterQC_violin_CD4CD8.pdf',200,270)
QCvio_clean

map(lungTcell,~ dim(.x))



# $CD4
# [1] 33538  2398
# 
# $CD8
# [1] 33538  3626


# normalization and scaling -----------------------------------------------
# scale and find high var genes 

lungTcell %<>% map(~   NormalizeData(.x, normalization.method = 'LogNormalize',
                               scale.factor = 10000, assay = 'RNA')%>% 
                     CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F) %>%
                     ScaleData( assay = 'RNA',
                                vars.to.regress = c("percent.mito",
                                                    "S.Score", 'G2M.Score',
                                                    "percent.ribo",
                                                    'nCount_RNA','nFeature_RNA' )) %>%
                     FindVariableFeatures(assay = 'RNA',nfeatures = 3000, selection.method = 'vst') )
lungTcell$p32@assays$RNA@var.features
for (x in proj) {
  lungTcell[[x]]@assays$RNA@var.features <- lungTcell[[x]]@assays$RNA@var.features%>%  
    str_subset('^RP|^MT|^HIST', negate = T)
}

saveRDS(lungTcell, 'lungTcell_pre_inte.rds')

map(lungTcell, ~ .x@assays$RNA@var.features %>% length() )

# integration -------------------------------------------------------------
anchors <- FindIntegrationAnchors(lungTcell, dims = 1:60)

CD4CD8 <- IntegrateData(anchors, dims = 1:60)

# mark CD4 and CD8 based on CITEseq ---------------------------------------

CD4CD8@assays$CITE@data
# change name of cite seq antibodies

CD4CD8 %<>% ScaleData( assay = 'CITE',vars.to.regress = c('patient'))
rownames(CD4CD8@assays$CITE@data) %<>% str_replace('-TotalSeqC', '.protein')
rownames(CD4CD8@assays$CITE@counts) %<>% str_replace('-TotalSeqC', '.protein')
rownames(CD4CD8@assays$CITE@scale.data) %<>% str_replace('-TotalSeqC', '.protein')



Feature_rast(CD4CD8, d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8, d1='CD4', d2='CD8B',  color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8,  g='patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T) +facet_grid(~ patient)+geom_hline(yintercept = 0.6)+geom_vline(xintercept = 0.8)

dim(CD4CD8)


Feature_rast(CD4CD8, c('CD4', 'CD8A', 'CD8B'))


CD4CD8 %<>% AddMetaData(FetchData(CD4CD8, c('CD4.protein', 'CD8.protein'), slot = 'data'))
# CD4 CD8 demultiplexing 

CD4CD8@meta.data %<>% mutate(CD4CD8= case_when(orig.ident == 'p25_CD4'~ 'CD4', orig.ident == 'p25_CD8'~ 'CD8', cite_CD4.protein >= 0.8 & cite_CD8.protein >= 0.6 ~ 'DP',cite_CD8.protein >= 0.6 ~ 'CD8',cite_CD4.protein >= 0.8 ~'CD4',cite_CD4.protein < 0.8 & cite_CD8.protein < 0.6 ~ 'DN'))

(Feature_rast(CD4CD8 %>% subset(patient != 'p25'),  g='CD4CD8', 'patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T) +facet_wrap(~ patient,ncol = 2)+ggtitle('CITEseq on CD4 and CD8')+geom_hline(yintercept = 0.6)+geom_vline(xintercept = 0.8) )%T>%  figsave('CD4andCD8staining_cite.pdf', 200, 200)

table(CD4CD8$CD4CD8)

Feature_rast(CD4CD8, 'CD4CD8','patient',  colorset = 'gg')
table(CD4CD8$CD4CD8)

ViolinPlot(CD4CD8, 'nCount_RNA', group.by = 'CD4CD8')

saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0716.rds')

# DP cells are significantly larger than other cells, thus we remove it as doublets
# CD4CD8  %<>% subset(subset = CD4CD8 %in% c('CD4', 'CD8', 'DN'))


Feature_rast(CD4CD8,c('CD103.protein', 'CD49a.protein', 'CD69.protein'), assay = 'CITE', color_grd = 'grd')


# dimensional reduction by PCA and UMAP -----------------------------------


DefaultAssay(CD4CD8) <- "integrated"
# scale data and run PCA
CD4CD8 %<>% ScaleData( vars.to.regress = c('patient', 
                                           "percent.mito",
                                                   "S.Score",
                                                   'G2M.Score',
                                                   "percent.ribo",
                                                   'nCount_RNA',
                                                   'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)




ElbowPlot(CD4CD8, ndims = 50)

# jackstraw is a way to choose significant PCs
CD4CD8 %<>% JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
CD4CD8 %>% JackStrawPlot( dims = 1:50 ) %T>%
  figsave('CD4CD8_lung_5p.jackstraw.pdf' , w = 400, h = 400)

CD4CD8 <- RunUMAP(object = CD4CD8, dims = 1:42, 
                   reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:42))
for (i in seq(0.6,1.5,0.1) %>% rev()) {
  CD4CD8 <- FindClusters(CD4CD8, resolution = i)
}

CD4CD8$integrated_snn_res.0.6

Feature_rast(CD4CD8, paste0('integrated_snn_res.',seq(0.6,1.3,0.1)), ncol =4) %T>% figsave('CD4CD8_different_reso.pdf', 400,200) 


Feature_rast(CD4CD8,c('patient','integrated_snn_res.1.4','tissue', 'CD4CD8','CCR6', 'RORC', 'FOXP3', 'DPP4', 'IL23R'), sz = 1)

Feature_rast(CD4CD8,'tissue',facet = 'patient')+facet_grid(~patient)

cite

CD4CD8@assays$CITE@counts

saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0716.rds')
Feature_rast(CD4CD8, c('CCR6', 'RORC', 'FOXP3', 'DPP4', 'IL23R', 'IL17A', 'IL17F', 'AREG', 'IFNG', 'IL22'), assay = 'RNA')


# clean and 2nd clustering ------------------------------------------------
# remove cells that not T cells 
Feature_rast(CD4CD8, noaxis = F, axis.number = T)
CD4CD8  %<>% subset(UMAP_1 >= -6)
DefaultAssay(CD4CD8) <- "RNA"

# CD4CD8<- PercentageFeatureSet(CD4CD8, '^HSPA', col.name =  'percent.hspa') 
DefaultAssay(CD4CD8) <- "integrated"
CD4CD8$percent.hspa
CD4CD8@assays$integrated@var.features
# scaling and 2nd round of clustering 

CD4CD8 %<>% ScaleData( vars.to.regress = c(
                                           "percent.mito",
                                           'percent.hspa',
                                           "S.Score",
                                           'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA',
                                           'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)

CD4CD8 %<>% JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
CD4CD8 %>% JackStrawPlot( dims = 1:50 ) %T>%
  figsave('CD4CD8_lung_5p.jackstraw.pdf' , w = 400, h = 400)
saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0716.rds')

# set.seed(123)
CD4CD8 <- RunUMAP(object = CD4CD8, dims = 1:45, seed.use = 1007,
                  reduction = 'pca', min.dist = 0.1) %>%
  FindNeighbors(dims = c(1:45))
for (i in seq(0.6,1.8,0.1) %>% rev()) {
  CD4CD8 <- FindClusters(CD4CD8, resolution = i, random.seed = 123)
}
Feature_rast(CD4CD8)

table(CD4CD8@active.ident)






ClusterCompare(CD4CD8, '17', '2')

Feature_rast(CD4CD8)
# decide to use resolution 1.4 for clustering 
Idents(CD4CD8)<- CD4CD8$integrated_snn_res.1.4

CD4CD8$Cell_cluster<- Idents(CD4CD8)

ClusterCompare(CD4CD8, '10', '7', group.by = 'integrated_snn_res.1.5')
ViFeature_rast(CD4CD8, c('CD4.protein', 'CD8.protein'), assay = 'CITE', color_grd = 'grd')
Feature_rast(CD4CD8, c('ident','CD4CD8'))


Feature_rast(CD4CD8)+facet_wrap('ident')


# # distribution between tissue and CD4CD9 --------------------------------


library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[5:9], brewer.pal(9,'Reds')[5:9])

# tissue

CD4CD8$ID <- paste0(CD4CD8$tissue,'_',CD4CD8$patient)

#calculate tissue composition table

comp_tissue <-CD4CD8@meta.data %>% group_by(tissue, ID, Cell_cluster) %>%  
  summarise(n = n() ) %>% mutate(n = case_when(tissue == 'luLN'~ -n,
                                               tissue == 'lung' ~ n)) %>% group_by(ID) %>% 
  mutate(percent = case_when(tissue == 'luLN'~ -(n/sum(n)*100),
                             tissue == 'lung' ~ n/sum(n)*100  )) %>% as.data.frame() 

comp_tissue


cl_comp_flow <- ggplot(comp_tissue, aes(y = n, x = Cell_cluster, fill = ID, color = ID,
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
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+
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
comp_CD %>% filter(Cell_cluster == 17) %>% as.data.frame()


comp_CD
library('tidyr')
library('broom')
comp_CD %>% group_by(Cell_cluster, patient) %>%   t_test(formula = abs(percent) ~ CD4CD8,
                                                         order = c("CD4", "CD8"),
                                                         alternative = "two-sided")
cl_comp_flow_CD <- ggplot(comp_CD, aes(y = percent, x = Cell_cluster, fill = ID_CD, color = ID_CD,
                                        stratum = ID_CD  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha = .6)+
  scale_fill_manual(values = ID_cl  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('% of cluster per donor')+
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+
  NULL
cl_comp_flow_CD


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
# DEGs --------------------------------------------------------------------
DefaultAssay(CD4CD8) <- 'RNA'


CD4CD8_DEGs <- 
  FindAllMarkers(CD4CD8, test.use = 'bimod', min.pct = 0.10,   only.pos = FALSE )%>%
        filter(p_val_adj < 0.05) %>%
        mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))


CD4CD8<- FindVariableFeatures(CD4CD8,assay = 'RNA',nfeatures = 2500, selection.method = 'vst')

# to reduce the site of suerat object, not every gene will be scaled, only the potential DEGs




gene_to_scale <- unique(c(CD4CD8_DEGs$gene, CD4CD8@assays$RNA@var.features))



CD4CD8 %<>% ScaleData( vars.to.regress = c('patient',
                                           'percent.hspa',
                                           "percent.mito",
                                           "S.Score",
                                           'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA',
                                           'nFeature_RNA' ), assay = 'RNA', feature = gene_to_scale)
CD4CD8_DEGs

ClusterCompare(CD4CD8,  '1', '8')

DoHeatmap(CD4CD8, c('CD4', 'CD8A', 'HLA-C','AL133415.1') ) 


                     

top10deg <- CD4CD8_DEGs %>% filter(!grepl('^RP|^MT', gene)) %>% 
                  arrange(cluster, desc(avg_log2FC))%>%
                  group_by(cluster) %>% top_n(10, avg_log2FC )


top10DEG_heat <- 
         DoHeatmap(subset(CD4CD8, downsample=500),features = top10deg$gene,raster = T, 
                   group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>% figsave('top10DEG.pdf', 190, 370)
top10DEG_heat










# Trm makrers 
TrmMarker <- FindAllMarkers(subset(CD4CD8, Cell_cluster %in% c(3,9,11,17)),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA')

TrmMarker


top15_TrmMarker <- TrmMarker %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(15, avg_log2FC ) %>% as.data.frame()

top10_TrmMarker

DoHeatmap(subset(CD4CD8, Cell_cluster %in% c(3,9,11,17)), features = top15_TrmMarker$gene) %>% heat_theme()

# Tnv markers
Tnvmarkers <-FindAllMarkers(subset(CD4CD8, Cell_cluster %in% c(7,8,10)),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA') 
top15_Tnvmarkers <- Tnvmarkers %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(15, avg_log2FC ) %>% as.data.frame()

DoHeatmap(subset(CD4CD8, Cell_cluster %in% c(7,8,10)), features = top15_Tnvmarkers$gene) %>% heat_theme()



# DEG CD4CD8 --------------------------------------------------------------

CD4CD8$CD_pheno <-   paste0(CD4CD8$CD4CD8, '_',CD4CD8$Cell_cluster)


# gene signatrue scors ----------------------------------------------------

CD4CD8@meta.data <- CD4CD8@meta.data[,1:75]
colnames(CD4CD8@meta.data)

sigtable <- read.xlsx('abd5778_Table_S3.xlsx',sheet = 1)

colnames(sigtable) %<>%   str_remove('.\\(.+\\)')

sigtable %<>% as.list() %>% map(~ na.exclude(.x) %>% as.vector)

names(sigtable)

sigtable$Tissue.resident


#caculate scores
for (i in names(sigtable)) {
  CD4CD8 <- AddModuleScore(CD4CD8, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}


colnames(CD4CD8@meta.data)[76:91] <- names(sigtable)


Feature_rast(CD4CD8,names(sigtable), color_grd = 'grd')

ViolinPlot(CD4CD8, names(sigtable), colors = umap.colors, box = T)

CD4CD8@assays$GM <- NULL
GMS  <- CD4CD8@meta.data[,c(names(sigtable))] %>% as.data.frame() %>% t()
scaleGM <- scale(t(CD4CD8@meta.data[,c(names(sigtable))]))
scaleGMassay <- CreateAssayObject(data = scaleGM)
CD4CD8@assays$GM <- scaleGMassay

DoHeatmap(subset(CD4CD8, downsample = 700), 
          raster =T, draw.lines = T, angle = 45,
          lines.width = 10,group.colors = umap.colors,
          
          assay = 'GM', features = names(sigtable), slot = 'data', size = gs(8)) +hmp2 + mytheme+
  theme(legend.position = 'bottom',
        legend.key.height = unit(2,'mm'))+
  guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))


# GSEA --------------------------------------------------------------------

library(clusterProfiler)

ent2_7 <- entrezlist_generator(CD4CD8, '2','7')


library(msigdbr)
Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)


c2_c7_Msigdb_GSEA <- GSEA(geneList = ent2_7, TERM2GENE=Mc7,  nPerm = 100000, 
                          minGSSize    = 15,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")

c2_c7_Msigdb_GSEA@result  <- c2_c7_Msigdb_GSEA@result %>% arrange(desc(NES)) 

view(c2_c7_Msigdb_GSEA@result)
Feature_rast(CD4CD8)
ClusterCompare(CD4CD8, 'C1', 'C2')


ClusterCompare(CD4CD8, 'C13', 'C17')




# CITEseq -----------------------------------------------------------------

Feature_rast(CD4CD8,  g='CD4CD8', 'patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T,slot = 'scale.data') +facet_wrap(~ patient,ncol = 2)+geom_hline(yintercept = 0.6)+geom_vline(xintercept = 0.8) 

Feature_rast(subset(CD4CD8, patient %in% c('p27', 'p71')),  g='ident', 'tissue',  d1='CD103.protein', d2='CD49a.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T, slot = 'scale.data', sz = 4)


ViolinPlot(CD4CD8, 'CD103.protein', assay = 'CITE')


Feature_rast(CD4CD8, c('CD103.protein', 'CD49a.protein', 'KLRG1.protein', 'CD69.protein'),
             assay = 'CITE', color_grd = 'grd')
Feature_rast(CD4CD8, c('CD161.protein', 'CD26.protein'),assay = 'CITE', color_grd = 'grd',slot = 'scale.data')


# RNA vs CITE
surfacemakers <- c('CD4','CD8A','CD8B', 'CXCR3','CCR6',  'ITGAE',   'CD69','CCR7','KLRB1','CD27','KLRG1', 'IL7R',
                   'DPP4', 'ITGA1', 'KLRD1', 'PTPRC', 'PDCD1')


surface_RNA <- Feature_rast(subset(CD4CD8, patient %in% c('p27', 'p71')), surfacemakers, sz = 0.3)

surface_CITE <- Feature_rast(subset(CD4CD8, patient %in% c('p27', 'p71')),sz = 0.3, CD4CD8@assays$CITE@data %>% rownames(), assay = 'CITE', color_grd = 'grd', slot = 'scale.data')

RNAvsCITE  <- PG(list(surface_RNA, surface_CITE), labels = c('RNA', 'CITE'), ncol = 1,vjust = -10) %T>%  figsave('RNAvsCITE.pdf',200, 290)

PG(list(surface_RNA, surface_CITE), labels = c('RNA', 'CITE'), ncol = 1,label_y = -2)

map(lungTcell, ~ .x@assays$CITE@counts %>% rownames )



Feature_rast(CD4CD8, c('AREG', 'CSF2', 'ITGAE', 'ITGA1', 'ZNF683', 'CD4CD8', 'ITGB1'))


# Cluster adjustment ------------------------------------------------------

Feature_rast(CD4CD8)


integrated_snn_res.1.4 <- unique(CD4CD8$integrated_snn_res.1.4) %>% as.numeric()  %>% sort() %>% as.character()
integrated_snn_res.1.4


Cell_pheno <-  c(
              'Tcm_CD4/8_L_LMNA.hi', #0
               'Trm_CD4/8_M_type1.hi',#1
                 'Tcm_CD4_L_CD69.hi',#2
              'Trm_CD8_T_LMNA.hi',#3
              'Tmem_CD8_M_MHCII.hi',#4
              'Trm_CD4/8_T_type3',#5
              'Th1_CD4/8_M',#6
              'Tnv_CD4/8_L_CD69.hi' ,#7
              'Tnv_CD4/8_L_CCR7.hi', #8
              'Trm_CD8_T_TXNIP.hi', #9
              'Tnv_CD4/8_L_HSPA.hi',#10
              'Trm_CD8_T_CD29.hi', #11
              'Th17_CD4_M',#12
              'Th2_CD4/8_M',#13
              'Treg_CD4_L',#14
              'Trm_CD4/8_T_IFN.induced',#15
              'Tfh_CD4_L',#16
              'Trm_CD8_T_KIRs.hi',#17
              'unclear_1',#18
              'unclear_1',#19
              'Trm_CD8_T_LMNA.hi',#20
              'unclear_2'#21
              
              
              
              
)

phenotable <- data.frame(integrated_snn_res.1.4, Cell_pheno) %>% arrange(Cell_pheno)

phenotable$Cell_cluster <- factor(paste0('C',c(1:17, 17, 18,19,19,20)), 
                                  levels = paste0('C', 1:20))

phenotable$Cluster_pheno <- paste0(phenotable$Cell_cluster,': ', phenotable$Cell_pheno) 
phenotable %<>% mutate(Cluster_pheno = factor(Cluster_pheno, levels = unique(Cluster_pheno)))



CD4CD8$Cell_cluster <- NULL

# CD4CD8@meta.data %<>%  left_join(phenotable, by = 'integrated_snn_res.1.4') %>% `rownames<-`(CD4CD8$bc_backup)
CD4CD8@meta.data %<>%  left_join(phenotable, by = 'integrated_snn_res.1.4', suffix = c('','')) %>% `rownames<-`(CD4CD8$bc_backup)

Feature_rast(CD4CD8,'Cluster_pheno', colorset = 'gg' )

Idents(CD4CD8) <- CD4CD8$Cell_cluster


UMAP_CD4CD8 <-(Feature_rast(CD4CD8, noaxis = F,sz = 0.3) +
                 ggtitle('T cells from lung tissue and lymph nodes')+
                 scale_color_manual(labels = levels(CD4CD8$Cluster_pheno), values = umap.colors)+
                 guides(color = guide_legend(ncol = 2, title = 'phenotypes', override.aes = list(size = 1.5)))) %T>% figsave('umap_CD4CD8_with_phenotypes.pdf',  180, 100  )



saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0716.rds')

# !!!!!!!!!!!since here we adjusted and renamed the clustering, some of scripts above from line 432 to 711 need to be run again



# # distribution between tissue and CD4CD9 --------------------------------


library(ggalluvial)
library(RColorBrewer)
RColorBrewer::display.brewer.all()
ID_cl <- c(brewer.pal(9,'Blues')[5:9], brewer.pal(9,'Reds')[5:9])
CD_cl <- c(brewer.pal(9,'RdPu')[4:8], brewer.pal(9,'Greens')[4:8])

# tissue

CD4CD8$ID <- paste0(CD4CD8$tissue,'_',CD4CD8$patient)
Feature_rast(CD4CD8, 'ID', colorset = CD_cl, do.label = F, noaxis = F)
Umap_donor_tissue <- Feature_rast(CD4CD8, 'ID', colorset = ID_cl, do.label = F, noaxis = F)+ ggtitle(NULL)+NoLegend()

#calculate tissue composition table

comp_tissue <-CD4CD8@meta.data %>% group_by(tissue, ID, Cell_cluster) %>%  
  summarise(n = n() ) %>% mutate(n = case_when(tissue == 'luLN'~ n,
                                               tissue == 'lung' ~ -n)) %>% group_by(ID) %>% 
  mutate(percent = case_when(tissue == 'luLN'~ (n/sum(n)*100),
                             tissue == 'lung' ~ -n/sum(n)*100  )) %>% as.data.frame() 

comp_tissue


cl_comp_flow <- ggplot(comp_tissue, aes(y = n, x = Cell_cluster, fill = ID, color = ID,
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
  guides(fill = guide_legend(nrow = 2, title = NULL, byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(15,'mm'),
        axis.line = element_blank())+
  NULL
cl_comp_flow

figtissue <- PG(list(Umap_donor_tissue, cl_comp_flow), rw = c(1, 1.8))

# CD4CD8

CD4CD8$ID_CD <- paste0(CD4CD8$CD4CD8,'_',CD4CD8$patient)

Umap_donor_CD <- Feature_rast(subset(CD4CD8, CD4CD8 %in% c('CD4', 'CD8')), 'ID_CD', colorset = CD_cl, do.label = F, noaxis = F)+ ggtitle(NULL)+NoLegend() 
Umap_donor_CD

comp_CD <-CD4CD8@meta.data %>% filter(CD4CD8 %in% c('CD4', 'CD8') ) %>% 
  group_by(CD4CD8, patient, Cell_cluster) %>%  
  summarise(n = n() ) %>%group_by(patient, Cell_cluster)  %>% 
  mutate(percent = case_when(CD4CD8 == 'CD8'~ -(n/sum(n)*100),
                             CD4CD8 == 'CD4' ~ n/sum(n)*100  ))  %>%
  mutate(n = case_when(CD4CD8 == 'CD8'~ -n,CD4CD8 == 'CD4' ~ n) ) %>% ungroup() %>%  tidyr::complete(Cell_cluster, patient,CD4CD8,fill = list(n = 0, percent = 0))%>%group_by(patient, Cell_cluster)  %>% mutate(ID_CD = paste0(CD4CD8, '_', patient), total = sum(abs(n))) %>% filter(total>0)


comp_CD
library('tidyr')
library('broom')
comp_CD %>% group_by(Cell_cluster, patient) %>%   t_test(formula = abs(percent) ~ CD4CD8,
                                                         order = c("CD4", "CD8"),
                                                         alternative = "two-sided")


meant_comp_CD<- comp_CD  %>% group_by(Cell_cluster, CD4CD8) %>% summarise( SD = sd(percent),percent = mean(percent))   %>% as.data.frame()
meant_comp_CD


comp_CD%>%ungroup %>%  group_by(Cell_cluster, patient)   %>%
  summarise_each(pv=funs(t.test(.[CD4CD8 == "CD4"], .[CD4CD8 == "CD8"])$p.value), vars=abs(percent))

library(ggpubr)

cl_comp_flow_CD<-ggplot(meant_comp_CD, aes(y = percent, x = Cell_cluster  ))+
  geom_bar(stat = 'identity', aes(fill = CD4CD8), color='black', size=0.3, width = 0.95)+
  # geom_errorbar(color = 'black', size = 0.2,aes(ymax = percent+SD, ymin=percent-SD))+
  theme_bw(base_line_size = 0)+
  geom_point(data=comp_CD, aes(color=ID_CD, y = percent, x = Cell_cluster, group = patient), position = position_dodge(width = 0.8) , size = 0.8 )+  
  scale_color_manual(values = c(CD_cl), labels = rep(sort(unique(CD4CD8$patient)),2)  )+
  scale_y_continuous(labels = abs)+
  scale_fill_manual(values = alpha(c('purple','green' ), 0.2))+
  theme(legend.position = 'bottom')+
  guides(  fill = guide_legend(nrow = 2, title = NULL),
    color = guide_legend(nrow = 2,  byrow = T, label.position = 'right', title = 'patient',override.aes = list(size = 2))
       )
cl_comp_flow_CD
figCD <- PG(list(Umap_donor_CD, cl_comp_flow_CD), rw = c(1, 1.8))


compositions <- PG(list(figtissue, figCD), ncol = 1) %T>% figsave('tissue_donor_CD4_8.pdf',200,220)

comp_CD2<- comp_CD %>% mutate(percent = abs(percent) )  %>% ungroup() %>% group_by(Cell_cluster)  %>% mutate(nrows = NROW(CD4CD8)) %>% filter(nrows >4 ) 

pvs <- comp_CD2 %>% ungroup() %>% group_by(Cell_cluster) %>% summarise(pv = t.test(abs(percent) ~ CD4CD8, paired = T)$p.value)

pvs$fdr <- p.adjust(pvs$pv, method = 'fdr', n = length(pvs$pv))

pvs$p.adjust <- p.adjust(pvs$pv, method = 'bonferroni', n = length(pvs$pv))
pvs

  




# DEGs --------------------------------------------------------------------
DefaultAssay(CD4CD8) <- 'RNA'


CD4CD8_DEGs <- 
  FindAllMarkers(CD4CD8, test.use = 'bimod', min.pct = 0.10,   only.pos = FALSE )%>%
  filter(p_val_adj < 0.05 | abs(log2FC) >0.5) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))
CD4CD8_DEGs %>% filter(avg_log2FC >0) %>% count(cluster)

CD4CD8_DEGs  %<>% filter(p_val_adj < 0.05 | abs(avg_log2FC) >0.5) 
write.xlsx(CD4CD8_DEGs, 'lungabT_DEGs.xls')

top10deg <- CD4CD8_DEGs %>%
  # filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC )


top10DEG_heat <- 
  DoHeatmap(subset(CD4CD8, downsample=500),features = top10deg$gene,raster = T, 
            group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>% figsave('top10DEG.pdf', 190, 370)
top10DEG_heat


# Trm makrers 
TrmMarker <- FindAllMarkers(subset(CD4CD8, Cell_cluster %in% paste0('C',13:18)),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA')

TrmMarker


top15_TrmMarker <- TrmMarker %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(15, avg_log2FC ) %>% as.data.frame()

 
DoHeatmap(subset(CD4CD8, Cell_cluster %in% paste0('C',13:18)), features = top15_TrmMarker$gene) %>% heat_theme()

# Tnv markers
Tnvmarkers <-FindAllMarkers(subset(CD4CD8, Cell_cluster %in% c(7,8,10)),min.pct = 0.1, logfc.threshold = 0.5,assay = 'RNA') 
top15_Tnvmarkers <- Tnvmarkers %>% filter(!grepl('^RP|^MT', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(15, avg_log2FC ) %>% as.data.frame()

DoHeatmap(subset(CD4CD8, Cell_cluster %in% c(7,8,10)), features = top15_Tnvmarkers$gene) %>% heat_theme()



# DEG CD4CD8 --------------------------------------------------------------

CD4CD8$CD_pheno <-   paste0(CD4CD8$CD4CD8, '_',CD4CD8$Cell_cluster)

Feature_rast(subset(CD4CD8, CD4CD8 %in% c('CD4', 'CD8')), 'CD_pheno', colorset = rainbow(40))


# gene signatrue scors ----------------------------------------------------


modulescore_Vln <-ViolinPlot(CD4CD8, names(sigtable), colors = umap.colors, box = T, x.angle = 45, ylabtext = ' score')

figsave(modulescore_Vln, 'modulescore_violin.pdf',180,270)


modulescore_Feature <-Feature_rast(CD4CD8, names(sigtable), color_grd = 'grd')



DoHeatmap(subset(CD4CD8, downsample = 700), 
          raster =T, draw.lines = T, angle = 45,
          lines.width = 10,group.colors = umap.colors,
          
          assay = 'GM', features = names(sigtable), slot = 'data', size = gs(8)) +hmp2 + mytheme+
  theme(legend.position = 'bottom',
        legend.key.height = unit(2,'mm'))+
  guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))



DEmodule <- FindAllMarkers(CD4CD8,assay = 'GM')
DEmodule %>% filter(avg_log2FC >0)

# GSEA --------------------------------------------------------------------

library(clusterProfiler)



library(msigdbr)
Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

ent13_17 <- entrezlist_generator(CD4CD8, 'C13','C17')
C1317_Msigdb_GSEA <- GSEA(geneList = ent13_17, TERM2GENE=Mc7,  nPerm = 100000, 
                          minGSSize    = 15,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
C1317_Msigdb_GSEA@result  <- C1317_Msigdb_GSEA@result %>% arrange(desc(NES)) 
view(C1317_Msigdb_GSEA@result)









# TCR data frome processing & integration  --------------------------------------------------------------------

# read raw TCR files 

colnames(CD4CD8) %>% str_extract('_.') %>% unique()
# [1] "_1" "_2" "_3" "_4" "_5" "_6"
proj
# [1] "p32"     "p45"     "p25_CD4" "p25_CD8" "p27"     "p71" 
pt


TCRdirs <- list.dirs(path = 'raw') %>% str_subset('VDJ.+outs$' )
# remove patient 31
TCRdirs <- TCRdirs[c(1:5,7)]

TCRs <- map2(TCRdirs, 1:6, ~  rbind(read.csv(paste0(.x, '/filtered_contig_annotations.csv' )),
                                    read.csv(paste0(.x, '/filtered_contig_annotations_2.csv' )))%>% 
               filter(productive == 'True'& is_cell == 'True') %>%
               dplyr::select(c(1, 5:10, 13,14))  %>%
               mutate( patient = pt[[.y]],  
                       bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% 
  set_names(proj) %>% reduce(.f = rbind)
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

saveRDS(TCRs_paired, 'abTCR_new.RDS')

# join the data 
CD4CD8$bc_backup <- rownames(CD4CD8@meta.data)


# in case there is duplicated data, remove old data at first

CD4CD8@meta.data  %<>%   select_at(.vars = vars(-contains(c('TRA','TRB','paired', '.x', '.y'))))

colnames(CD4CD8@meta.data)

CD4CD8@meta.data %<>% left_join(TCRs_paired, by =c('bc_backup', 'patient')) %>% `rownames<-`(CD4CD8$bc_backup)




Feature_rast(CD4CD8, 'cdr3_paired_perc',sz = 0.5)+
  scale_color_gradient2( na.value = alpha('lightgrey',0.5),
                          low = "#00ccff",midpoint = 1,
                         mid = 'purple', 
                         high = "#ff0066"   )
Feature_rast(CD4CD8, c('cdr3_TRB_perc', 'tissue'), color_grd = 'grd', facets = 'patient')



saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0924.rds')
dim(CD4CD8@assays$RNA@counts)

# TCR analysis -----------------------------------------------------------
# mapping rate
# to see how many T cells are paired with a TCR 


CD4CD8@meta.data %<>% mutate(TCR_summary = case_when(!is.na(paired) ~ 'paired TCR',
                                                    !is.na(chain_TRA)~'single TRA',
                                                    !is.na(chain_TRB)~ 'single TRB') )

Feature_rast(CD4CD8, 'TCR_summary', do.label = F)
table(CD4CD8$TCR_summary)

TCRmapping <- CD4CD8@meta.data %>% group_by(tissue, patient) %>%  count(TCR_summary) %>% 
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
saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0924.rds')


CD4CD8$cdr3_paired_freq

CD4CD8@meta.data %<>% mutate(clonal_expansion =case_when(cdr3_paired_freq == 1 ~ 'monoclonal',
                                                       nr(cdr3_paired_freq, 2,4)~'low (2~4)',
                                                       nr(cdr3_paired_freq, 5,9)~'moderate (5~9)',
                                                       
                                                       cdr3_paired_freq >9 ~ 'high (>9)' ),
                             clonal_expansion = factor(clonal_expansion, levels = c('high (>9)', 'moderate (5~9)',
                                                                                    'low (2~4)', 'monoclonal'))   )

clone_exp_umap <- Feature_rast(CD4CD8, 'clonal_expansion',c('patient', 'tissue'),do.label = F, colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'), facetcol = 4) %T>% figsave('clonalexpansion_tissue_patient.pdf',270,160) 

Feature_rast(CD4CD8, 'clonal_expansion',do.label = F, colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'))%T>% figsave('TCR_clonalexpansion_UMAP.pdf',120,100) 

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

Feature_rast(CD4CD8, 'cdr3_TRA_perc', color_grd = 'grd')
CD4CD8$cdr3_paired_perc


# TCR VDJ -----------------------------------------------------------------

CD4CD8@meta.data %<>% mutate(TRAVJ = case_when(!is.na(v_gene_TRA)  & !is.na(j_gene_TRA) ~ paste(v_gene_TRA, j_gene_TRA)  ) 
                            )

Feature_rast(subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20')) ,
            c( 'TRAVJ'), noaxis = F, sz = 1, facets = 'patient')

Feature_rast(subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20')) ,
             c( 'TRAVJ', 'v_gene_TRB',  'cdr3_TRA_perc'), noaxis = F, sz = 1)

CD4CD8$tissue


BG = Feature_rast(CD4CD8, 'tissue', colorset = c('grey', 'grey'), do.label = F)+NoLegend()




BG +(geom_point_rast(data =    subset(CD4CD8, v_gene_TRA == 'TRAV1-2' & j_gene_TRA %in% c('TRAJ12', 'TRAJ30', 'TRAJ20' ))  %>% FetchData(c('UMAP_1', 'UMAP_2', 'TRAVJ', 'c_gene_TRB', 'patient', 'tissue')
), aes(x = UMAP_1, y = UMAP_2, color = TRAVJ, size  = 1) +scale_color_manual(values = ggplotColours(12)) )                  
                     )+ggtitle('TRAV1-2 J12 TRBV6')+ facet_grid('patient')


DimPlot(CD4CD8)


Feature_rast(CD4CD8, 'patient')

Feature_rast(CD4CD8, c('tissue', 'KLRB1', 'DPP4', 'RORC' , 'CD4', 'CD8A'), color_grd = 'grd', ncol =3)


# TCRsharing --------------------------------------------------------------



# TCR sharing between CD4 CD8 &  lung and LN

C14TCR <- CD4CD8@meta.data %>% filter(!is.na(cdr3_paired) & CD4CD8 %in% c('CD4', 'CD8')) %>% select(cdr3_paired,CD4CD8) %>% group_split(CD4CD8)


intersect(C14TCR[[1]]$cdr3_paired,C14TCR[[2]]$cdr3_paired)



tissueTC <- CD4CD8@meta.data %>% filter(!is.na(cdr3_paired) &CD4CD8 %in% c('CD8')) %>% select(cdr3_paired,tissue) %>% group_split(tissue)


TCRpatient <- CD4CD8@meta.data %>% filter(!is.na(cdr3_paired) ) %>% select(cdr3_paired,patient) %>% group_split(patient) %>%set_names(patientID) %>% map(~ .x %>% pull(cdr3_paired) %>% unique)
  
(TCRpatient %>% unlist() %>% table() %>% sort(decreasing = T) )[1:100]

intersect(tissueTC[[1]]$cdr3_paired,tissueTC[[2]]$cdr3_paired)

TCRbyCD4CD8 <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1  & CD4CD8 %in% c('CD4', 'CD8')) %>% select(cdr3_paired,CD4CD8) %>%group_by(CD4CD8,cdr3_paired) %>%  summarise(pairedfreq = n()) %>% arrange(CD4CD8) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )
  arrange(pairedfreq)
nrow(TCRbyCD4CD8 %>% filter())


TCRBbyCD4CD8 <-CD4CD8@meta.data %>% filter(cdr3_TRB_freq >1  & CD4CD8 %in% c('CD4', 'CD8')) %>% select(cdr3_TRB,CD4CD8) %>%group_by(CD4CD8,cdr3_TRB) %>%  summarise(pairedfreq = n()) %>% arrange(CD4CD8) %>% mutate(cdr3_TRB = factor(cdr3_TRB, unique(cdr3_TRB)) )
cdr3_TRB

library(ggalluvial)
TCRsharingCD4CD8 <- (TCRbyCD4CD8 %>% 
  ggplot(
    aes( x = CD4CD8 , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    ggtitle("TCR sharing between CD4 and CD8")+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(884))+
    geom_stratum(size = 0.1, color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)
TCRsharingCD4CD8

(TCRBbyCD4CD8 %>% 
    ggplot(
      aes( x = CD4CD8 , y = pairedfreq, fill = cdr3_TRB,  stratum= cdr3_TRB, alluvium  = cdr3_TRB))+
    ggtitle("TCR sharing between CD4 and CD8")+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(1484))+
    geom_stratum(size = 0.1, color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)



TCRbytissue <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 &CD4CD8 %in% c('CD8', 'CD4')) %>% select(cdr3_paired,tissue,CD4CD8) %>%group_by(tissue,cdr3_paired,CD4CD8) %>%  summarise(pairedfreq = n()) %>% arrange(CD4CD8, tissue) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )

nrow(TCRbytissue %>% filter())

TCRsharingtissue <- (TCRbytissue %>% 
                       ggplot(
                         aes( x = tissue , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
                       ggtitle("TCR sharing between lung and LN")+
                       geom_flow(stat = "alluvium",
                                 color = "darkgray") +
                       # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
                       theme_minimal_hgrid()+
                       scale_fill_manual(values = rainbow(1014))+facet_wrap(~CD4CD8)+
                       geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
                       xlab(NULL) +ylab("TCRab frequencies")+
                       theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
                       guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300)


TCRsharingtissue

TCRsharingplot <- PG(list(TCRsharingCD4CD8, TCRsharingtissue), rw = c(1,1.8)) %T>% figsave('TCRsharingplot_CD4_8_tissue.pdf', 150,120)

  

Feature_rast(CD4CD8)
# how TCRs shared within clusters?

TCRbytissue_cluster <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 &CD4CD8 %in% c('CD8', 'CD4')) %>%
  select(cdr3_paired,tissue,CD4CD8, Cell_cluster) %>%group_by(tissue,cdr3_paired,CD4CD8,Cell_cluster) %>%  summarise(pairedfreq = n()) %>% ungroup() %>% arrange(CD4CD8, tissue) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )
TCRbytissue_cluster

clusters <- c('C1', 'C2', 'C4',  'C6', 'C7', 'C12', 'C14', 'C15', 'C17', 'C18')

clusters2 <- c('C1', 'C2', 'C4',  'C6', 'C7', 'C12', 'C14', 'C15', 'C17', 'C18')


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


TCRby_cluster <-CD4CD8@meta.data %>% filter(cdr3_paired_freq >1 &CD4CD8 %in% c('CD8', 'CD4')) %>%
  # select(cdr3_paired,tissue,CD4CD8, Cell_group) %>%
  group_by(cdr3_paired,CD4CD8,Cell_group, patient) %>%  summarise(pairedfreq = n()) %>% ungroup() %>% arrange(CD4CD8) %>% mutate(cdr3_paired = factor(cdr3_paired, unique(cdr3_paired)) )

CD4groups <- c('group 2: CD4 helpers', 'group 3: Tcm and Th1','group 4: Tem','group 6: Th1_17 Trm')
CD8groups <- c('group 6: Th1_17 Trm','group 3: Tcm and Th1','group 4: Tem','group 5: Trm effector')

(TCRby_cluster %>% 
    filter(Cell_group %in% CD4groups & CD4CD8 == 'CD4' ) %>% mutate(Cell_group = factor(Cell_group, CD4groups)) %>% 
    ggplot(
      aes( x = Cell_group , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(1100))+
    facet_wrap(~patient, scales = "free", ncol = 5)+
    geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 

TCRsharing_CD8_cluster <-(TCRby_cluster %>% 
    filter(Cell_group %in% CD8groups & CD4CD8 == 'CD8' ) %>% mutate(Cell_group = factor(Cell_group, CD8groups)) %>% 
    ggplot(
      aes( x = Cell_group , y = pairedfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
    geom_flow(stat = "alluvium",
              color = "darkgray") +
    theme_minimal_hgrid()+
    scale_fill_manual(values = rainbow(1100)[280:810])+
    facet_wrap(~patient, scales = "free", ncol = 5)+
    geom_stratum(size = 0.05,color = alpha('black', 0.5))+ 
    xlab(NULL) +ylab("TCRab frequencies")+ggtitle('CD8 TCRab')+
    theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
    guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 
TCRsharing_CD8_cluster
# public data satija ------------------------------------------------------

lung_Satija <- readRDS('/home/big/tanlikai/Lung/public/Azimuth.lung.Satija.rds')



CD8diffcopd <- ClusterCompare(lung_Satija %>% subset(annotation.l1 == 'CD8 T'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

Mfdiffcopd <- ClusterCompare(lung_Satija %>% subset(annotation.l1 == 'Macrophage'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

Basaldiffcopd <- ClusterCompare(lung_Satija %>% subset(annotation.l1 == 'Basal'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

CD14monodiffcopd <-
  ClusterCompare(lung_Satija %>% subset(annotation.l1 == 'CD14+ Monocyte'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', do.plot = F)

batcheff_copd_normal <- intersect(CD8diffcopd$table$gene, Basaldiffcopd$table$gene) %>% intersect(CD14monodiffcopd$table$gene) 
intersect(CD8diffcopd$table$gene, CD14monodiffcopd$table$gene)



intersect(CD8diffcopd$table$gene, Basaldiffcopd$table$gene) %>% intersect(CD14monodiffcopd$table$gene)  %>% sort()

Feature_rast(lung_Satija, c('tissue','TRDC', 'TRGC1', 'TRGC2', 'CD3D', 'CD3E') ,colorset = 'gg')
Feature_rast(lung_Satija, c('TRAC','TRDC', 'TRGC1', 'TRGC2', 'CD3D', 'CD3E') ,sz = 0.2)

lung_Satija$tissue %>% table()
tcells <-  lung_Satija$cell_type %>% unique() %>% str_subset('T cell|thy')

diseases  <-  c('chronic obstructive pulmonary disease', 'normal', 'COVID-19')

lungT_Satija <-  readRDS('/home/big/tanlikai/Lung/public/Azimuth.lung.Satija.rds') %>% 
  subset(cell_type %in% c(.$cell_type %>% unique() %>% str_subset('T cell|thy')) & disease %in% diseases)
dim(lungT_Satija)

colnames(lungT_Satija@meta.data)


Feature_rast(lungT_Satija, c( 'CD3D', 'CD3E'), colorset = 'gg')
Feature_rast(lungT_Satija, c( 'disease', 'health_status', 'tissue'), colorset = 'gg', sz = 0.1)
lungT_Satija$disease %>% table
table(lungT_Satija$disease, lungT_Satija$health_status)


Feature_rast(lungT_Satija, c('CD4', 'CD8A', 'RORC', 'CCR6', 'FOXP3', 'TOX', 'ITGAE' , 'ITGA1'))



load('/home/big/tanlikai/Lung/public/GSE162498_NSCLC_CD3_4tumors_4Juxta_2Juxta.Rds')
Feature_rast(eleven.tils.cd3.integrated, 'IL17A')
Feature_rast(eleven.tils.cd3.integrated)



three.tissues.cd3.integrated$tissue %>% unique()
dim(three.tissues.cd3.integrated)

Feature_rast(three.tissues.cd3.integrated, c('CD4', 'CD8A', 'CD8B', 'ITGAE', 'ITGA1'))


three.tissues.cd3.integrated@assays$RNA@counts %>% dim()


table(three.tissues.cd3.integrated$tissue)





# integration with satija T celss -----------------------------------------

# clean the data  

lungT_Satija  %<>% 
PercentageFeatureSet( '^MT', col.name =  'percent.mito') %>% 
  PercentageFeatureSet('^RP', col.name = 'percent.ribo')    %>% 
  PercentageFeatureSet('^HSPA', col.name =  'percent.hspa')



ViolinPlot(lungT_Satija, 'nCount_RNA', box = T) +ylim(0, 10000)
ViolinPlot(lungT_Satija, 'nFeature_RNA', box = T) +ylim(0, 5000)
ViolinPlot(lungT_Satija, 'percent.mito', box = T) 

Feature_rast(lungT_Satija, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', color_grd = 'grd',
             
             noaxis = F, axis.number = T)+
  geom_smooth(method = "lm")+
  
  scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
  scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
  geom_hline(yintercept = c(500,2500))+geom_vline(xintercept = c(1000,6500))



lungT_Satija  %<>% subset( nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA > 1000 & nCount_RNA < 6500& percent.mito < 16)

dim(lungT_Satija)
lungT_Satija$donor %>%  unique()
table(lungT_Satija$donor) %>% sort(decreasing = T)

table(lungT_Satija$donor, lungT_Satija$disease)


TCELLgenes <- c('PTPRC','CD3D', 'CD3G', 'CD3E', 'TRAC', 'TRBC1', 'TRBC2','TRDC', 'TRGC1', 'TRGC2', 'CD4', 'CD8A', 'CD8B')

Feature_rast(lungT_Satija, c('donor', 'disease', 'tissue'), colorset = 'gg', ncol =1)
Feature_rast(lungT_Satija, TCELLgenes, noaxis = F)
ViolinPlot(lungT_Satija, c('CD3D', 'CD3G', 'CD3E', 'TRAC', 'TRBC1', 'TRBC2','TRDC', 'TRGC1', 'TRGC2'),colors = umap.colors)


# select real T cells 
lungT_Satija  %<>%  AddModuleScore(features = list(c('CD3D', 'CD3G', 'CD3E')),name = 'CD3score')

ViolinPlot(lungT_Satija, c('CD3score1'),colors = umap.colors)+geom_hline(yintercept = c(-0.1,0.6))

Feature_rast(lungT_Satija, 'CD3score1', color_grd = 'grd', sz = 0.1)

lungT_Satija$bc_backup <- rownames(lungT_Satija@meta.data)  

lungT_Satija@meta.data  %<>% mutate(CD3exp = case_when(CD3score1 < 0 ~ 'CD3neg',
                                                       nr(CD3score1, 0, 0.6) ~ 'CD3lo',
                                                       CD3score1 > 0.6 ~ 'CD3hi')) %>% `rownames<-`(lungT_Satija$bc_backup)


Feature_rast(lungT_Satija, 'CD3exp',  sz = 0.5, colorset = 'gg', noaxis = F, axis.number = T) +
  geom_hline(yintercept = c(0,9))+
  geom_vline(xintercept = c(-11,0))

table(lungT_Satija$CD3exp)



lungT_Satija  %<>%  subset(CD3exp %in% c('CD3lo', 'CD3hi') & UMAP_1 > -11 & UMAP_1 < 0 & UMAP_2 > 0 & UMAP_2 < 9 )

dim(lungT_Satija)

Feature_rast(lungT_Satija, c('donor', 'tissue', 'disease'), colorset = 'gg', ncol =1)


lungT_Satija  %<>%  FindVariableFeatures(selection.method = 'vst', nfeatures = 5000)

lungT_Satija@assays$RNA@data


vf1 <- lungT_Satija@assays$RNA@var.features
vf2 <- lungT_Satija@assays$RNA@var.features
intersect(vf1, vf2)

lungT_Satija@assays$RNA@var.features %<>%str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T) 
VariableFeaturePlot(lungT_Satija) +ylim(-1, 10)
VariableFeaturePlot(CD4CD8)

lungT_Satija@assays$RNA@var.features
multicores(mem = 200)
Feature_rast(lungT_Satija, 'annotation.l1')
lungT_Satija %<>%  ScaleData(verbose = FALSE, assay = 'RNA',feature = lungT_Satija@assays$RNA@var.features,

                       vars.to.regress = c('percent.mito',  'nCount_RNA','nFeature_RNA', 'donor' ) )


table(lungT_Satija$dataset_origin, lungT_Satija$disease)


fwlungT_Satija  %<>%  RunPCA(npcs =  100) %>% RunUMAP(dims = 1:60, 
                                                    reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:60)) %>% 
  FindClusters(resulution = 0.6)

Feature_rast(lungT_Satija, c('ident', 'tissue', 'disease'), colorset = 'gg')

saveRDS(lungT_Satija, 'lungT_Satija_cleaned_processed.rds')

table(lungT_Satija$donor, lungT_Satija$disease)


ClusterCompare(lungT_Satija %>% subset(annotation.l1 == 'CD8 T'), 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease')

lungT_Satija@assays$RNA@scale.data

table(lungT_Satija$donor)

table(lungT_Satija$tissue, lungT_Satija$disease)

lungT_Satija$project <- 'satija'



# integration -------------------------------------------------------------

CD4CD8@assays
CD4CD8$donor <- CD4CD8$patient
CD4CD8$project <- 'Prinz&Falk'
CD4CD8$dataset_origin <- 'Prinz&Falk'


CD4CD8$pj_tissue <- paste(CD4CD8$project, CD4CD8$tissue)
lungT_Satija$pj_tissue <- paste(lungT_Satija$project, lungT_Satija$tissue)

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


lungT_Satija.l  <- SplitObject(lungT_Satija, split.by = "dataset_origin")

map(lungT_Satija.l, ~ dim(.x))

lungT_Satija.l$mayr_2020 <- NULL
lungT_Satija.l$lukassen_2020 <- NULL


lungT_Satija.l <-   lapply(X = lungT_Satija.l, FUN = function(x) {
  # x  <- AddMetaData(x, t(x[['CITE']]@data)) 

  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

for (x in names(lungT_Satija.l))      {
  
  lungT_Satija.l[[x]]@assays$RNA@var.features %<>% str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T)  %>% 
    setdiff(batcheff_copd_normal)
  
}


multicores(mem = 400)

anchors <- FindIntegrationAnchors(object.list = append(CD4CD8.l,lungT_Satija.l ), 
                                  dims = 1:50)
lung.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)


lung.integrated$project


DefaultAssay(lung.integrated) <- "integrated"

lung.integrated@meta.data  %<>% mutate(dataset_origin = case_when(project == 'Prinz&Falk' ~ 'Prinz&Falk',
                                                                  project == 'satija' ~ dataset_origin)) 
# scale data and run PCA
lung.integrated %<>% ScaleData( vars.to.regress = c('donor', 
                                                    'dataset_origin',
                                           "percent.mito",
                                           # 'project',
                                           "percent.ribo",
                                           'nCount_RNA',
                                           'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)



ClusterCompare(lung.integrated, 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease')

ClusterCompare(lungT_Satija, 'chronic obstructive pulmonary disease', 'normal', group.by = 'disease', rm = 'NO-NAME|RP|MT|HIST|HSPA', genetoshow = 100)

ElbowPlot(lung.integrated,ndims = 100)


Feature_rast(lungT_Satija, 'disease')


lung.integrated <- RunUMAP(object = lung.integrated, dims = 1:70, 
                  reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:70))
for (i in seq(0.6,1,0.1) %>% rev()) {
  lung.integrated <- FindClusters(lung.integrated, resolution = i)
}

ClusterCompare(lung.integrated, '3', '7', do.plot = F)


Feature_rast(lung.integrated, c('tissue', 'CD4CD8','pj_tissue', 'disease', 'project'), ncol = 2, sz  = 0.2)

Feature_rast(lung.integrated, 'ident', colorset = 'gg')

lung.integrated$donor %>% unique()


lung.integrated$tissue %>% unique()

lung.integrated$CD4CD8


Feature_rast(lung.integrated, c('CD4', 'CD8A', 'CD8B', 'KLRB1', 'KLRG1', 'ITGA1',  'TRGC1', 'TRDC', 'FOXP3',
                                'ITGAE', 'CD69' ,'RORC', 'CCR6','PDCD1', 'TOX') , assay = 'RNA')



ClusterCompare(lung.integrated, '3', '4')
saveRDS(lung.integrated, 'lungT.integrated.prinz.satija.rds')

# harmony -----------------------------------------------------------------
library(harmony)
lung.integrated@meta.data  %<>% mutate(dataset_origin = case_when(project == 'Prinz&Falk' ~ 'Prinz&Falk',
                                                                  project == 'satija' ~ dataset_origin)) %>% 
  `rownames<-`(lung.integrated$bc_backup)
DefaultAssay(lung.integrated) <- 'RNA'

lung.integrated  %<>%  FindVariableFeatures(selection.method = 'vst', nfeatures = 5000)



lung.integrated@assays$RNA@var.features %<>%str_subset('^NO-NAME|^MT|^IG|^TRAV|^TRBV|^HSP|^RP', negate = T)  %>% 
  setdiff(batcheff_copd_normal)




lung.integrated %<>% ScaleData( vars.to.regress = c(
  # 'donor', 
                                                    "percent.mito",
                                                    'dataset_origin',
                                                    # 'G2M.Score',
                                                    "percent.ribo",
                                                    'nCount_RNA',
                                                    'nFeature_RNA' ), assay = 'RNA', feature = lung.integrated@assays$RNA@var.features) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40) 


lung.integrated %<>%  RunHarmony(group.by.vars = c('dataset_origin', 'donor')) %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(preduction = "harmony", dims = 1:40) %>% FindClusters(resolution = 0.6)
ElbowPlot(lung.integrated, reduction = 'harmony',ndims = 100)

Feature_rast(lung.integrated, c('tissue', 'CD4CD8','pj_tissue', 'disease', 'project'), ncol = 2, sz  = 0.2)

  
Feature_rast(lung.integrated, c('disease', 'dataset_origin', 'Cell_pheno'))


ClusterCompare(lung.integrated, '1', '6')
