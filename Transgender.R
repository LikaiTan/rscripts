
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

# creat project -----------------------------------------------------------


dir_all <- list.files('cellranger_out/')
# TCR library
dir_TCR <- str_subset(dir_all, 'tcr')
# gdTCR library
dir_gdTCR<- str_subset(dir_all, 'gdTCR')
# gene expression library
dir_GEX <- str_subset(dir_all, 'tcr|gdTCR', negate = T)

samplesheet <- data.frame(sample = dir_GEX, donor  = str_extract(dir_GEX, 'sample0\\d\\d'), tp = str_extract(dir_GEX, '\\w\\w$'))





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

Feature_rast(TG_Casar, 'TCR_ab_gd')


TG_Casar@meta.data$nFeature_RNA


ViolinPlot(TG_Casar, 'nFeature_RNA', group.by = 'TCR_ab_gd')

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



