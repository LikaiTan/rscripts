
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
 
TG_Casar@meta.data$bc_backup <- rownames(TG_Casar@meta.data)

TG_Casar@meta.data %>%  group_by(orig.ident) %>% slice(1L) %>% pull(bc_backup)

colnames(TG_Casar@meta.data )

TG_Casar@meta.data %>% head

TG_Casar@meta.data  %<>%  mutate(TCR = case_when(
  !is.na(CTgene) ~ 'abTCR'
))
Feature_rast(TG_Casar, 'TCR')



TG_Casar@meta.data %>%  head()

# creat project -----------------------------------------------------------


dir_all <- list.files('cellranger_out/')
# TCR library
dir_TCR <- str_subset(dir_all, 'tcr')
# gdTCR library
dir_gdTCR<- str_subset(dir_all, 'gdTCR')
# gene expression library
dir_GEX <- str_subset(dir_all, 'tcr|gdTCR', negate = T)

samplesheet <- data.frame(sample = dir_GEX, donor  = str_extract(dir_GEX, 'sample0\\d\\d'), tp = str_extract(dir_GEX, '\\w\\w$'))


TG_Casar$bc_backup


# introduce gdTCR ---------------------------------------------------------

dir_gdTCR
raw_gdTCR <- map(dir_gdTCR, ~ read.csv(paste0('cellranger_out/',.x, '/outs/all_contig_annotations.csv' ))%>% 
                   filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene))  %>%
                   dplyr::select(c(1, 5:10, 13,14))  %>%
                   mutate(  
                     bc_backup = paste0(gsub('_gdTCR', '_', .x), gsub('-1', '', barcode)))%>% dplyr::select(-barcode)
                 
                 )%>% reduce(.f = rbind)

raw_gdTCR


TRGs <- raw_gdTCR %>% filter(grepl('GV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
  # mutate(v_gene = str_remove(v_gene,'DV\\d')) %>%
  rename_at(vars(-bc_backup), funs(sub('$','_TRG',.)))
TRGs

TRDs <- raw_gdTCR %>% filter(grepl('DV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
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
TCRs<- TCRs %>% 
  mutate(paired_sp =
           case_when(paired %in% top9paired ~ paired, 
                     !is.na(paired) ~ "Other paired")    ) %>%
  mutate(paired_sp = factor(paired_sp, levels = c(top9paired, "Other paired")))
# paired TRG and J 
TCRs<- TCRs %>% mutate(TRGJP = case_when(!is.na(chain_TRG) ~ paste(v_gene_TRG, j_gene_TRG)))%>%
  mutate(TRGJP = str_replace(TRGJP, ' TRG',' '))


TG_Casar@meta.data  %<>% left_join(TCRs, by = 'bc_backup', suffix = c('','')) %>% `rownames<-`(TG_Casar@meta.data$bc_backup )

TCRs


TG_Casar@meta.data  %<>%  mutate(TCR_ab_gd = case_when(
  !is.na(CTgene)& ( !is.na(chain_TRD ) | !is.na(chain_TRG )) ~ 'dual-TCR',
  !is.na(CTgene) ~ 'abTCR',
  !is.na(chain_TRD ) | !is.na(chain_TRG )  ~ 'gdTCR'
))

Feature_rast(TG_Casar, 'TCR_ab_gd', facets = 'patient')

Feature_rast(TG_Casar, 'TCR_ab_gd')

table(TG_Casar$TCR_ab_gd, TG_Casar$orig.ident)

TG_Casar@meta.data$nFeature_RNA

Feature_rast(TG_Casar, noaxis = F,axis.number = T)


ViolinPlot(TG_Casar, 'nCount_RNA', group.by = 'TCR_ab_gd', colors = umap.colors, box = T)

DefaultAssay(TG_Casar) <- 'RNA'


Feature_rast(TG_Casar, c('TRAC', 'TRBC2', 'TRDC', 'TRGC1'),  titlesize = 12)


# gene expression ---------------------------------------------------------

# read raw gene expression data  
RAW_GEX <- map(dir_GEX, ~ Read10X_h5(paste0('cellranger_out/', .x,'/outs/filtered_feature_bc_matrix.h5')) ) %>%  setNames(dir_GEX)


names(RAW_GEX)


ClusterCompare(TG_Casar, 'CD8_CL2', 'CD8_CL3')


TG_gdT@assays$RNA@counts %>%  head() %>%  view()



names(RAW_GEX$sample001_BL)
# [1] "Gene Expression"  "Antibody Capture"

RAW_GEX$sample001_BL$`Gene Expression` %>%  str()

map(RAW_GEX, ~ dim(.x$`Gene Expression`))


RAW_GEX$sample001_BL$`Antibody Capture` %>%  str()


rownames(RAW_GEX$sample001_BL$`Antibody Capture`)


bcraw_gdTCR <- map(dir_gdTCR, ~ read.csv(paste0('cellranger_out/',.x, '/outs/all_contig_annotations.csv' ))%>% 
                   filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene))  %>%
                   dplyr::select(c(1, 5:10, 13,14))  
                 
) %>%  setNames(dir_gdTCR)


intersect(bcraw_gdTCR$sample002_BL_gdTCR$barcode, colnames(RAW_GEX$sample002_BL$`Gene Expression` ))

bcoverlap <- list()


for (a in 1:8) {
  for (j in 1:8) {
    bcoverlap[[paste0(a,j)]]<-   intersect(bcraw_gdTCR[[a]]$barcode, colnames(RAW_GEX[[j]]$`Gene Expression` )) %>% length()
    
  }
}

bcoverlap %>% as.data.frame()


dir_gdTCR

names(RAW_GEX)


bcraw_gdTCR[[1]]['barcode']

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


# only gdT ----------------------------------------------------------------

TG_gdT <- subset(TG_Casar, TCR_ab_gd == 'gdTCR')
Feature_rast(TG_gdT, 'condition')




ViolinPlot(TG_gdT, 'nFeature_RNA')


TG_gdT@active.assay

TG_gdT$S.Score

TG_gdT %<>% ScaleData(vars.to.regress = c("percent.mito", "patient",'percent.ribo','S.Score', 'G2M.Score' ))

dim(TG_gdT)

TG_gdT  %<>%  RunPCA(npcs= 100, nfeatures.print = 40) %>%  JackStraw(dims = 80) %>%
  ScoreJackStraw(dims = 1:80) 
JackStrawPlot(TG_gdT, dims = 1:80 )%T>%
  figsave('Transgender_gdt.jackstraw.pdf' , w = 400, h = 400,device = 'pdf')




# UMAP and clustering -----------------------------------------------------




TG_gdT <- RunUMAP(TG_gdT, dims = 1:37, 
                    reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:28))
Feature_rast(TG_gdT, 'condition')
Feature_rast(TG_gdT, 'condition', facets = 'patient')

Feature_rast(TG_gdT, 'v_gene_TRD', facets = 'patient')

TG_gdT$orig.ident


TRDmapp <- TG_gdT@meta.data %>% count(orig.ident, v_gene_TRD) %>%  group_by(orig.ident)  %>%
  mutate(percent = n/sum(n)*100) 
TRGmapp <- TG_gdT@meta.data %>% count(orig.ident,v_gene_TRG) %>%  group_by(orig.ident)  %>%
  mutate(percent = n/sum(n)*100) 
mappedTCR <- map2(list(TRDmapp, TRGmapp), list('v_gene_TRD','v_gene_TRG'),~
                    ggplot(.x, aes_string('orig.ident', y = "n", color = .y, fill = .y))+
                    geom_bar(stat = 'identity',  width = 0.5,position = position_stack(reverse = T))+
                    color_m()+
                    fill_m()+
                    xlab(NULL)+
                    ylab('%')+
                    ggtitle(.y)+
                    guides(color = F, fill = guide_legend(nrow =4, title = NULL))+
                    theme_minimal()+
                    mytheme+gglp('r')
                  
)


# TCR delta sharing between two time points 

TCRD_freq <- TG_gdT@meta.data %>% group_by(patient) %>%  count(cdr3_TRD) %>% `colnames<-`(c('patient', 'cdr3_TRD','TRD_freq'))

TCRD_percent <-  TG_gdT@meta.data %>% group_by(orig.ident)  %>%  count(cdr3_TRD)  %>%  mutate(TRD_perc = n/sum(n)*100 ) %>% select(-n)


sum(TCRD_percent$TRD_perc)

TCRD_freq
TG_gdT@meta.data %<>% left_join(TCRD_freq, c('patient', 'cdr3_TRD')  , suffix = c('', '') ) %>% `rownames<-`(TG_gdT$bc_backup)

TG_gdT@meta.data %<>% left_join(TCRD_percent, c('orig.ident', 'cdr3_TRD')  , suffix = c('', '') ) %>% `rownames<-`(TG_gdT$bc_backup)


Feature_rast(TG_gdT, 'TRD_perc', color_grd = 'grd', facets = 'condition', titlesize = 12,sz = 1.5)


Feature_rast(TG_gdT, c('ident', 'v_gene_TRD'), noaxis = F)

Feature_rast(TG_gdT, 'cdr3_TRD',colorset = 'gg',do.label = F) +NoLegend()


TCRD_condition <- TG_gdT@meta.data %>%filter(!is.na(v_gene_TRD) &TRD_freq>1) %>%  
  group_by(condition, patient) %>%  count(cdr3_TRD) %>% mutate(percent = n/sum(n)*100)%>%
arrange(desc(percent))


ggplot(data = TCRD_condition,
aes( x = condition,
y = percent, fill = cdr3_TRD,  stratum= cdr3_TRD,
alluvium  = cdr3_TRD))+
# ggtitle(paste0('GV9DV2 TCR'))+
geom_flow(stat = "alluvium",size = 0.2,
color = "darkgray") +
scale_y_continuous(breaks  = seq(0,100,20), limits = c(0,100))+
theme_minimal_hgrid()+
# scale_fill_manual(values = set_sample(ggplotColours(45)))+
geom_stratum(size = 0.1)+
xlab(NULL) +ylab("% of repertoire")+mytheme+
theme(axis.text.x = element_text(angle = 45,size =8, vjust = 0.5))+
theme(legend.position = 'none')+
facet_grid(~patient)



# phenotypes --------------------------------------------------------------

for (i in seq(0.5,0.9,0.1) %>% rev()) {
  TG_gdT <- FindClusters(TG_gdT, resolution = i)
}

TG_gdT$integrated_snn_res.0.9
Feature_rast(TG_gdT,g= 'integrated_snn_res.0.9', facets = c( 'patient', 'condition'))

Feature_rast(TG_gdT,g= paste0('integrated_snn_res.0.',5:9))

Idents(TG_gdT) <- TG_gdT$integrated_snn_res.0.5



table(TG_gdT$patient, TG_gdT$condition)
table(TG_Casar$patient, TG_Casar$condition)


DefaultAssay(TG_gdT) <- 'RNA'

ClusterCompare(TG_gdT, '1', '4',min.pct = 0.1, log2fc = 0.25)

TG_gdT %<>% ScaleData( vars.to.regress = c(
  "percent.mito",

    "S.Score",
  'G2M.Score',
  "patient",
  'nCount_RNA',
  'nFeature_RNA' ))



ClusterCompare(TG_gdT, group.by = 'condition', id1 = 'BL', id2 = 'V2')




# gene modules ------------------------------------------------------------


library(openxlsx)


gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()
gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)
gene_cluster$GM_B


gcanno <- as.vector(c('GM_A: Naive or immature T cell',
                      "GM_B: Innate T cell differentiation",
                      'GM_C: Proliferating',
                      "GM_D: Type-3 immunity",
                      "GM_E: Interferon induced",
                      "GM_F: CTL response",
                      "GM_G: CTL response (VD1)",
                      
                      "GM_H: Acute activation"))


for (i in gc_name) {
  TG_gdT %<>%  AddModuleScore( features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(TG_gdT@meta.data) %<>% str_replace("(?<=\\w)1", '')
  
}
Feature_rast(TG_gdT, c(gc_name), color_grd = 'grd') 
ViolinPlot(TG_gdT, gc_name,colors = umap.colors, box = T)


c(4,6) %>% map(~
                       Feature_rast(TG_gdT, gc_name[.x], color_grd = 'grd', mythe = F, titlesize  = 12)+
                       ggtitle(gcanno[.x])
) %>% PG(ncol = 2)


Feature_rast(TG_gdT, c('ident', 'CCR7','CD27',  'KLRB1', 'KLRD1', 'RORC', 'TBX21', 'ZBTB16', 'GZMB'))

Feature_rast(TG_gdT, 'v_gene_TRD')



Feature_rast(TG_gdT,c('ident', 'v_gene_TRD'), noaxis = F, titlesize = 12)







Feature_rast(TG_gdT, c('CD196', 'CD197', 'CD127', 'CD161', 'CD279'), assay = 'ADT', color_grd = 'grd')
Feature_rast(TG_gdT, c('CD196', 'CD197', 'CD127', 'CD161', 'CD183'), assay = 'ADT', color_grd = 'grd', facets = 'condition', ncol = 3)




ClusterCompare(TG_gdT, 'BL', 'V2',group.by = 'condition', assay = 'ADT')


cl_condition <- TG_gdT@meta.data %>%
  group_by(condition, patient) %>%  count(seurat_clusters) %>% mutate(percent = n/sum(n)*100)%>%
  arrange(desc(percent))


ggplot(data = cl_condition,
       aes( x = condition,
            y = percent, fill = seurat_clusters,  stratum= seurat_clusters,
            alluvium  = seurat_clusters))+
  # ggtitle(paste0('GV9DV2 TCR'))+
  geom_flow(stat = "alluvium",size = 0.2,
            color = "darkgray") +
  scale_y_continuous(breaks  = seq(0,100,20), limits = c(0,100))+
  theme_minimal_hgrid()+fill_m()+
  geom_stratum(size = 0.1)+
  xlab(NULL) +ylab("% of repertoire")+mytheme+
  theme(axis.text.x = element_text(angle = 45,size =8, vjust = 0.5))+
  theme(legend.position = 'right')+
  facet_grid(~patient)



TG_gdT$condition_cl <- paste0(TG_gdT$condition, '_', TG_gdT$seurat_clusters)


ClusterCompare(TG_gdT, 'BL_3', 'V2_3', group.by = 'condition_cl', log2fc = 0.25, min.pct = 0.1)

Feature_rast(TG_gdT, 'CXCR3', facets = 'condition')


