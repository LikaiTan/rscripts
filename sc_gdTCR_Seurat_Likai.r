# gdTCR incorporation script for Seurat object
# auothor: Likai Tan

library(Seurat)
library(purrr)
library(dplyr)
library(stringr)
library(magrittr)



library(Nebulosa)


# read TCR files ----------------------------------------------------------


# here we have TCR data from two different datasets: 1 and 2_3
# and we will combine them as one
# for each sample  we have two sequence run.

# 1.read the dir where you put your TCR csv 
TCRdirs <- list.dirs(path = 'raw') %>% str_subset('VDJ.+outs$' )

# reasult of TCRdirs, there are four TCR libraries
# [1] "raw/VDJ_322_325_453_456_gd/outs" "raw/VDJ_Falk1_gd/outs"           "raw/VDJ_falk3_gd/outs"          
# [4] "raw/VDJ_falk4_73_77_gd/outs"    


# put the barcodes of Seurat object as a column of the meta.data
# in this case, the Seuratobject contains data from FOUR libraries, so the barcode end with _1, _2, _3, _4
#  you need to check the barcodes format in TCR data and Seurat data to see if they are identical

Seuratobject$bc_backup <- rownames(Seuratobject@meta.data)


# read raw TCR data and combine into one dataframe
# the c(2,1,3,4) is the order I wish to read the libraries so the TCR data and Seurat data will match


AllTCRs <- map2(TCRdirs, c(2,1,3,4), ~ list.files(path = .x,pattern =  c('all_contig_annotations.csv'), full.names = T) %>% 
                read.csv()) %>% 
  # filter out low quality seqs and non-GDTCR seqs
                  filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene) &high_confidence  == 'True') %>%
  # discard useless column, can be skipped
                  dplyr::select(c(1, 5:10, 13,14))  %>%
  # modify barcode format to match with Seurat barcodes
                  mutate(   bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% 
  # combine to one data frame
  reduce(.f = rbind)


# instead, if there is only ONE library

AllTCRs <-    read.csv('thepathtoyourTCRcsv')        %>% 
                  filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene) &high_confidence  == 'True') %>%
                  dplyr::select(c(1, 5:10, 13,14))  %>%
                  mutate(    bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)



# split the datafrom the TRG and TRD, discard duplicates.

TRGs <- AllTCRs %>% filter(grepl('GV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
  # mutate(v_gene = str_remove(v_gene,'DV\\d')) %>%
  rename_at(vars(-bc_backup), funs(sub('$','_TRG',.)))
TRGs


TRDs <- AllTCRs %>% filter(grepl('DV', v_gene))%>% distinct(bc_backup, .keep_all = T) %>% 
  rename_at(vars(-bc_backup), funs(sub('$','_TRD',.))) %>%
  # # for TRD, I removed the TRAV gene name so it looks simpler, can skip
  mutate(v_gene_TRD = str_remove(v_gene_TRD, 'AV\\d\\d'))
TRDs


# put them back as one datafrome.
TCRs <- full_join(TRGs,TRDs, by = c('bc_backup'))

# looking for paired TCR
TCRs  %<>%  mutate(paired = 
                     case_when(!is.na(v_gene_TRD) & !is.na(v_gene_TRG) ~ 
                                 paste0(str_remove(v_gene_TRG,'TR'),' ',str_remove(v_gene_TRD,'TR'))))
TCRs  %<>% mutate(cdr3_paired = 
                    case_when(!is.na(v_gene_TRD) & !is.na(v_gene_TRG) ~ 
                                paste(str_remove(v_gene_TRG,'TR'),cdr3_TRG,str_remove(v_gene_TRD,'TR'),cdr3_TRD))) 


# # joinTCR to seurat --------------------------------------------------------
# GDTlung is the name of seurat object

# highly suggest to backup your medadata before incorporate TCR data into Seurat, so if anything messed up you can get it back.
metabackup <- Seuratobject@meta.data 


#  incorporate the TCR
Seuratobject@meta.data  %<>% left_join(TCRs, by = 'bc_backup', suffix = c('','')) %>% `rownames<-`(Seuratobject@meta.data$bc_backup )

Seuratobject@meta.data

#   select(!ends_with('.x')) 
# mapped TCR
# to calculated how many cells is with a TCR, ID here is the patient ID in myobject.

TRDmapp <- Seuratobject@meta.data %>% count(ID, v_gene_TRD) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 
TRGmapp <- Seuratobject@meta.data %>% count(ID,v_gene_TRG) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 


pairedmapp <-Seuratobject@meta.data %>% count(ID,paired) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 



Seuratobject@meta.data %<>% mutate(TCR_summary = case_when(!is.na(paired) ~ 'paired TCR',
                                                        !is.na(chain_TRG)~'single TRG',
                                                        !is.na(chain_TRD)~ 'single TRD') )


mappedTCR <- map2(list(TRDmapp, TRGmapp, pairedmapp), list('v_gene_TRD','v_gene_TRG', 'paired'),~
                    ggplot(.x, aes_string('ID', y = "percent", color = .y, fill = .y))+
                    geom_bar(stat = 'identity',  width = 0.5,position = position_stack(reverse = T))+
                    xlab(NULL)+
                    ylab('%')+
                    ggtitle(.y)+
                    guides(color = F, fill = guide_legend(nrow =4, title = NULL))+
                    theme_minimal()
                  
) %>% PG(nrow =1, labels = 'AUTO')
mappedTCR

TCRmapping <- Seuratobject@meta.data %>% group_by(tissue, patient) %>%  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)


Seuratobject@meta.data  %<>%  mutate(Vg9Vd2 = case_when(paired_sp == 'GV9 DV2' ~ 'GV9 DV2',
                                                     !is.na(paired_sp)  ~ 'other γδTCR'
)) %>% `row.names<-`(Seuratobject$bc_backup)






