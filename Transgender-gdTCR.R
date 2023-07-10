
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


TG_Casar@meta.data$bc_backup <- rownames(TG_Casar@meta.data)

TG

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



DimPlot(TG_Casar, group.by = 'TCR_ab_gd')



# only gdT ----------------------------------------------------------------

TG_gdT <- subset(TG_Casar, TCR_ab_gd == 'gdTCR')



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
DimPlot(TG_gdT,group.by =  'condition')

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

FeaturePlot(TG_gdT, 'TRD_perc',split.by = 'condition')

TCRD_condition <- TG_gdT@meta.data %>%filter(!is.na(v_gene_TRD) &TRD_freq>1) %>%  
  group_by(condition, patient) %>%  count(cdr3_TRD) %>% mutate(percent = n/sum(n)*100)%>%
arrange(desc(percent))
library(ggalluvial)

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




