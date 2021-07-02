
# lung CD4 and CD8 project
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
library(openxlsx)
library(stringr)
library(magrittr)
# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')


file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

#paralell computation
options(future.fork.enable = TRUE)
options(future.globals.maxSize= 1024*1024*1024*100)
future::plan(strategy = "multicore", workers = 60)

# work dir ----------------------------------------------------------------
setwd('/home/big/tanlikai/Lung/abt')


# read_rawdata ---------------------------------------------------------------

dirs <- list.dirs(path = 'raw') %>% str_subset('surface.+outs$' )

dirs

proj <- c('p32', 'p45','p25_CD4', 'p25_CD8','p27','p31','p71')
pt <- c('p32', 'p45','p25', 'p25','p27','p31','p71')

# read raw data 
rawdata <- dirs %>% map(~ Read10X_h5(paste0(.x, '/filtered_feature_bc_matrix.h5')) ) %>% setNames(proj)
# The quality of p31 is very poor, so we need remove it.
rawdata$p31 <- NULL
proj <- c('p32', 'p45','p25_CD4', 'p25_CD8','p27','p71')
pt <- str_extract(proj, 'p\\d\\d')
lungTcell <- map2(rawdata,proj, ~
                   CreateSeuratObject(counts = .x$`Gene Expression`,project = .y,
                                       min.features = 20) %>% 
                    AddMetaData(str_extract(.y, 'p\\d\\d'),col.name = 'patient')
                 )




map(lungTcell, ~ .@meta.data %>% colnames )
map(lungTcell, ~ .@meta.data %>% ncol )

map(lungTcell, ~ dim(.))


lungTcell$p31 <- NULL

# demultiplexing ----------------------------------------------------------



# introduce citeseq and hashtaq

hash <- unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )) %>%  str_subset('^\\d_')
hash

cite <- setdiff(unlist(map(rawdata, ~ .$`Antibody Capture` %>% rownames )), hash)
cite

for (x in proj) {
  lungTcell[[x]][['HTO']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(hash, rownames(rawdata[[x]]$`Antibody Capture`)),])
  lungTcell[[x]][['CITE']] <- CreateAssayObject(rawdata[[x]]$`Antibody Capture`[intersect(cite, rownames(rawdata[[x]]$`Antibody Capture`)),])
  rownames(lungTcell[[x]][['HTO']]@counts) <- paste0('X',rownames(lungTcell[[x]][['HTO']]@counts) )
  rownames(lungTcell[[x]][['HTO']]@data) <- paste0('X',rownames(lungTcell[[x]][['HTO']]@data) )
  
  
}

for (x in proj) {
  rownames(lungTcell[[x]][['HTO']]@counts) %<>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+') %>% paste0('S',.)
  rownames(lungTcell[[x]][['HTO']]@data)%<>% str_extract('(?<=\\w-).+') %>% str_extract('(?<=\\w-).+')%>% paste0('S',.)
  
  
}

hash_sample <- lungTcell %>% map(~ .[['HTO']] %>% rownames)
hash_sample
lungTcell$p32@assays$HTO@counts


for (i in proj) {
  lungTcell[[i]] %<>%  NormalizeData( assay = "HTO", normalization.method = "CLR") %>% 
    ScaleData(assay='HTO', features = rownames( .[['HTO']])) %>% 
    # RunUMAP(assay= 'HTO', features = rownames( .[['HTO']])) %>%
    # RunPCA(assay= 'HTO', features = rownames( .[['HTO']])) %>% 
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

lungTcell$p32@assays$HTO@counts
Feature_rast(lungTcell$p32, 'hash.ID', d1 = hash_sample$p32[[1]], d2=hash_sample$p32[[2]])


lungTcell$p32$hash.ID
map2(lungTcell,hash_sample,   ~ Feature_rast(.x, d1 =.y[[1]], d2= .y[[2]],
                                             g = 'hash.ID',noaxis =F, axis.number=T)+
       xlim(0,5)+ylim(0,5)+
       geom_hline(yintercept = 0.5)+
       geom_vline(xintercept = 0.8)) %>% PG(ncol =6) %T>%
  figsave('hash.id.automatic.pdf', 400, 100)
Feature_rast(lungTcell$p32, 'hash.ID')

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
lungTcell <- map(lungTcell, ~ PercentageFeatureSet(.x, '^MT', col.name =  'percent.mito') %>% 
             PercentageFeatureSet('^RP', col.name = 'percent.ribo')   )
  
QCvio <- map(lungTcell, ~ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                colors =umap.colors  ,box = T, jitter = T , ncol = 3  )) %>% 
    PG(labels = proj, ncol = 1) %T>%
  figsave('beforeQC_violin.pdf', 150, 270)

  QCvio

  figsave(QCvio,'beforeQC_violin_CD4CD8.pdf', 200, 200)
  

QC_scatter <- map(lungTcell, ~ Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', facet = 'tissue',
                                      noaxis = F, axis.number = T)+grd+facet_wrap(~tissue)+
                    scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
                    scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
                    geom_hline(yintercept = c(400,2000))+geom_vline(xintercept = c(1200,7000)))  %>% 
  PG(labels = proj)%T>%
  figsave('beforeQC_scatter.pdf',600,400)
QC_scatter

figsave(QC_scatter,'beforeQC_scatter.pdf',600,300)


saveRDS(lungTcell, 'lungTcell_pre_inte.rds')
# trim
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
    str_subset('^TRAC|^TRBC|^TRAV|^TRBV|^MT|^HIST', negate = T)
}


map(lungTcell, ~ .x@assays$RNA@var.features %>% length() )

# integration -------------------------------------------------------------
anchors <- FindIntegrationAnchors(lungTcell, dims = 1:50)

CD4CD8 <- IntegrateData(anchors, dims = 1:50)

# mark CD4 and CD8 based on CITEseq ---------------------------------------



CD4CD8 %<>% ScaleData( assay = 'CITE',vars.to.regress = c('patient'))
rownames(CD4CD8@assays$CITE@scale.data) %<>% str_replace('-TotalSeqC', '.protein')

Feature_rast(CD4CD8, c('CD4.protein','CD8.protein'), assay = 'CITE', color_grd = 'grd')
Feature_rast(CD4CD8, c('CD4','CD8A'), color_grd = 'grd')

CD4CD8@assays$CITE


Feature_rast(CD4CD8, d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8, d1='CD4', d2='CD8B',  color_grd = 'grd', noaxis = F)
Feature_rast(CD4CD8,  g='patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T) +facet_grid(~ patient)+geom_hline(yintercept = 0.6)+geom_vline(xintercept = 0.8)

dim(CD4CD8)


CD4CD8 %<>% AddMetaData(FetchData(CD4CD8, c('CD4.protein', 'CD8.protein'), slot = 'data'))

CD4CD8@meta.data %<>% mutate(CD4CD8= case_when(orig.ident == 'p25_CD4'~ 'CD4', orig.ident == 'p25_CD8'~ 'CD8', cite_CD4.protein >= 0.8 & cite_CD8.protein >= 0.6 ~ 'DP',cite_CD8.protein >= 0.6 ~ 'CD8',cite_CD4.protein >= 0.8 ~'CD4',cite_CD4.protein < 0.8 & cite_CD8.protein < 0.6 ~ 'DN'))

Feature_rast(CD4CD8,  g='CD4CD8', 'patient', d1='CD4.protein', d2='CD8.protein', assay = 'CITE', color_grd = 'grd', noaxis = F,  axis.number = T) +facet_grid(~ patient)+geom_hline(yintercept = 0.6)+geom_vline(xintercept = 0.8)


Feature_rast(CD4CD8, 'CD4CD8','patient',  colorset = 'gg')+facet_grid(~patient)

table(CD4CD8$CD4CD8)

ViolinPlot(CD4CD8, 'nCount_RNA', group.by = 'CD4CD8')

saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0619.rds')

# DP cells are significantly larger than other cells, thus we remove it as doublets
CD4CD8  %<>% subset(CD4CD8 %in% c('CD4', 'CD8', 'DN'))



# dimensional reduction by PCA and UMAP -----------------------------------


DefaultAssay(CD4CD8) <- "integrated"

CD4CD8 %<>% ScaleData( vars.to.regress = c('patient', 
                                           "percent.mito",
                                                   "S.Score",
                                                   'G2M.Score',
                                                   "percent.ribo",
                                                   'nCount_RNA',
                                                   'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)
CD4CD8 %<>% ScaleData( vars.to.regress = c('patient',
                                           "percent.mito",
                                           "S.Score",
                                           'G2M.Score',
                                           "percent.ribo",
                                           'nCount_RNA',
                                           'nFeature_RNA' ), assay = 'RNA')




ElbowPlot(CD4CD8, ndims = 50)


CD4CD8 %<>% JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
CD4CD8 %>% JackStrawPlot( dims = 1:50 ) %T>%
  figsave('CD4CD8_lung_5p.jackstraw.pdf' , w = 400, h = 400)
saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0619.rds')

CD4CD8 <- RunUMAP(object = CD4CD8, dims = 1:41, 
                   reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:41))
for (i in seq(0.6,1.5,0.1) %>% rev()) {
  CD4CD8 <- FindClusters(CD4CD8, resolution = i)
}

CD4CD8$integrated_snn_res.0.6

Feature_rast(CD4CD8, paste0('integrated_snn_res.',seq(0.6,1.3,0.1)), ncol =4) %T>% figsave('CD4CD8_different_reso.pdf', 400,200) 


Feature_rast(CD4CD8,c('patient','ident','tissue', 'CD4CD8'), sz = 1)

Feature_rast(CD4CD8,'tissue',facet = 'patient')+facet_grid(~patient)

cite

CD4CD8@assays$CITE@counts



# # distribution between tissue and CD4CD9 --------------------------------



# DEGs --------------------------------------------------------------------
DefaultAssay(CD4CD8) <- 'RNA'


ClusterCompare(CD4CD8, '3', '5')

DoHeatmap(CD4CD8, c('CD4', 'CD8A', 'HLA-C','AL133415.1') ) 


CD4CD8_DEGs <- 
  FindAllMarkers(CD4CD8, test.use = 'bimod', min.pct = 0.10,   only.pos = FALSE )%>%
        filter(p_val_adj < 0.05) %>%
        mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))
                     

top10deg <- CD4CD8_DEGs %>% filter(!grepl('^RP|^MT', gene)) %>% 
                  arrange(cluster, desc(avg_log2FC))%>%
                  group_by(cluster) %>% top_n(10, avg_log2FC )


top10DEG_heat <- 
         DoHeatmap(CD4CD8,features = top10deg$gene,raster = T, 
                   group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>% figsave('top10DEG.pdf', 190, 270)
top10DEG_heat

surfacemakers <- c('TIGIT','PDCD1','ITGAE','ITGA1',  'CD69','CCR7','KLRB1','CD27','KLRG1','KLRD1', 'IL7R',
                  'CCR6',  'DPP4','RORC','NCR1','NCR3','FCGR3A', 'FOXP1')


sfab <-  map(ABT, ~Feature_rast(.x, c('ident','tissue',  surfacemakers), ncol = 4)
)
  
  
figsave(sfab$CD4, 'sfmarkers_CD4t.pdf', 200,200)
figsave(sfab$CD8, 'sfmarkers_CD8t.pdf', 200,200)

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

TCRs <- map2(TCRdirs, 1:6, ~ read.csv(paste0(.x, '/filtered_contig_annotations.csv' ))%>% 
               filter(productive == 'True'& is_cell == 'True') %>%
               dplyr::select(c(1, 5:10, 13,14))  %>% 
               mutate( patient = pt[[.y]],  
                       bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% 
  set_names(proj) %>% reduce(.f = rbind)
TCRs
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


TCRs_paired %<>% left_join(cdr3TRA_freq, by = c('cdr3_TRA','patient')) %>% 
  left_join(cdr3TRB_freq, by = c('cdr3_TRB','patient')) %>%
  left_join(cdr3Paired_freq, by = c('cdr3_paired','patient'))



saveRDS(TCRs_paired, 'abTCR.RDS')

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
Feature_rast(CD4CD8, 'cdr3_TRB_perc', color_grd = 'grd', facet = 'patient')+  facet_grid(~patient)

saveRDS(CD4CD8, 'CD4CD8_integrated_2021_0619.rds')
dim(CD4CD8@assays$RNA@counts)

# TCR analysis -----------------------------------------------------------
# mapping rate
# to see how many T cells are paired with a TCR 


CD4CD8@meta.data %<>% mutate(TCR_summary = case_when(!is.na(paired) ~ 'paired TCR',
                                                    !is.na(chain_TRA)~'single TRA',
                                                    !is.na(chain_TRB)~ 'single TRB') )

Feature_rast(CD4CD8, 'TCR_summary', colorset = 'gg')


TCRmapping <- CD4CD8@meta.data %>% group_by(tissue, patient) %>%  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)
 (ggplot(TCRmapping,aes(x = tissue,y= percent, group = TCR_summary, fill = TCR_summary))+geom_bar(stat = 'identity',position = position_stack(reverse = T))+facet_grid(~patient)+fill_m()+theme_bw()+mytheme) %T>% figsave('mappingrate_TCR.pdf', 150, 60)





