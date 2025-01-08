
# lung project 6 patients -------------------------------------------------





library(Seurat)
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
library(Nebulosa)

# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')


file.edit('//home/big/tanlikai/script/rscripts/funcs.r')

# future::plan(strategy = "sequential")
# ?parallelly::supportsMulticore
#single core
setwd("/home/big/tanlikai/Lung")

# GDTlung_s <- SeuratDisk::LoadH5Seurat('GDTlung.h5Suerat.h5seurat')

GDTlung_s  %<>% subset(patient != "p31")

GDTlung.rds <- 'GDTlung2023july_7p.rds'



GDTlung_s <- readRDS(GDTlung.rds)
saveRDS(GDTlung_s, GDTlung.rds)


GDTlung_s@meta.data  %<>%  mutate(
  across(everything(), ~replace(.x,.x == 'NA', NA))
  
)


saveRDS(GDTlung_s,GDTlung.rds)

GDTlung_s$Vg9Vd2

# before integration ------------------------------------------------------

file.edit('GDTlung_before_integration.r')


# QC ----------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes



GDTlung %<>% map(., ~ PercentageFeatureSet(.x, '^MT|^Mt', col.name =  'percent.mito') %>% 
                   PercentageFeatureSet('^RP', col.name = 'percent.ribo') #%>%
                 # CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
)

GDTlung$p1$percent.mito

GDTlung$p23$Phase
QCvio_T <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                           colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>% 
  PG(ncol = 2) %T>% figsave('beforeQC_count_MT.pdf', 150, 200) 

QCvio_T


QCvio_P <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'patient',
                                           colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>% PG(ncol = 3)

QCvio_P


figsave(QCvio_T,'beforeQC_violin_8patient.pdf', 200, 200)

QC_scatter <- GDTlung %>%  map( ., ~Feature_rast(.x, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', facet = 'tissue',
                                                 noaxis = F, axis.number = T)+grd+facet_wrap(~tissue)+
                                  scale_x_continuous(breaks = seq(0, 10000, 1000), limits = c(0,10000))+
                                  scale_y_continuous(breaks = seq(0, 5000, 500), limits = c(0,5000))+
                                  geom_hline(yintercept = c(400,2000))+geom_vline(xintercept = c(1000,7000))) %>%
  PG(ncol = 3)


QC_scatter
figsave(QC_scatter,'beforeQC_scatter.pdf',600,250)

# clean 
GDTlung %<>% map(.,~ subset(.x, subset =  nCount_RNA %in% 1100:7000 &
                              nFeature_RNA %in% 400:2000 &  percent.mito <25&
                              tissue %in% c('Lung', 'LLN')) 
)
dim(GDTlung$p1)
dim(GDTlung$p23)
# cell number
map(GDTlung, ~ dim(.x))
#mean counts 
map(GDTlung,~ .x$nCount_RNA %>% mean)


QCvio_T_clean <- GDTlung %>% map(.,~  ViolinPlot(.x,c("nFeature_RNA", "nCount_RNA",'percent.mito'),group.by = 'tissue',
                                                 colors =umap.colors  ,box = T, jitter = T , ncol = 1  ) ) %>%
  PG(ncol = 2) %T>% figsave('afterQC_violin_8patient.pdf',150,150)

QCvio_T_clean

# figsave(QCvio_T_clean,'afterQC_violin_6patient.pdf',150,150)

saveRDS(GDTlung,'GDTlung_8p_beforeintegration.rds')


GDTlung <- readRDS("GDTlung_8p_beforeintegration.rds") 

# integration by Seurat ---------------------------------------------------



# normalization -----------------------------------------------------------
GDTlung%<>% map(.,~ NormalizeData(.x, normalization.method = 'LogNormalize',
                                  scale.factor = 10000, assay = 'RNA') %>%
                  ScaleData( assay = 'RNA',
                             vars.to.regress = c("percent.mito",
                                                 # "S.Score", 'G2M.Score',
                                                 "percent.ribo",
                                                 'nCount_RNA','nFeature_RNA' )) %>%
                  FindVariableFeatures(assay = 'RNA',nfeatures = 2500, selection.method = 'vst')
)

for (i in c('p1','p23', 'p456', 'p78')) {
  GDTlung[[i]]@assays$RNA@var.features%<>%  
    str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T)
}





# data integration --------------------------------------------------------
anchors <- FindIntegrationAnchors(GDTlung, dims = 1:60)

GDTlung_s <- IntegrateData(anchors, dims = 1:60)

DefaultAssay(GDTlung_s) <- "integrated"


GDTlung_s %<>% ScaleData(vars.to.regress = c("percent.mito", 
                                             'patient',
                                             "Phase",
                                             # 'G2M.Score',
                                             # "percent.ribo",
                                             'nCount_RNA' ))

GDTlung_s  %<>%   FindVariableFeatures(assay = 'RNA',nfeatures =4000, selection.method = 'vst')
# GDTlung_s@assays$RNA@var.features%<>%  
#   str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T)



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


GDTlung_s %<>%  CellCycleScoring( s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

D
# PCA  --------------------------------------------------------------------
dim(GDTlung_s)

GDTlung_s %<>% RunPCA( npcs = 100, verbose =
                         T, nfeatures.print = 40) %>% 
  JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 

ElbowPlot(GDTlung_s, ndims = 50)

JackStrawPlot(GDTlung_s, dims = 1:50 ) %>%
  figsave('GDTlung_8p.jackstraw.pdf' , w = 400, h = 400)
# saveRDS(GDTlung,GDTlung.rds)



# dimensional reduction & clustering --------------------------------------

GDTlung_s  %<>%  RunUMAP( dims = 1:35, 
                          reduction = 'pca') %>%
  FindNeighbors(dims = c(1:35))

# GDTlung_s %<>% RunUMAP(dims = 1:37)

for (i in seq(0.8,1.5,0.1) ) {
  GDTlung_s  %<>% FindClusters( resolution = i)
}
GDTlung_s$RNA_snn_res.1.2
# Feature_rast(GDTlung_s, c('ident', 'pheno', 'tissue', 'Cell_cluster'))


# Feature_rast(GDTlung_s, 'Cell_cluster')

GDTlung_s$ID <- paste0(GDTlung_s$tissue, GDTlung_s$patient)

Feature_rast(GDTlung_s,c('ident','tissue','patient', 'ID'), ncol = 2, noaxis = F)
ViolinPlot(GDTlung_s,'nFeature_RNA', group.by = 'integrated_snn_res.0.9')
DefaultAssay(GDTlung_s) <- "RNA"


Feature_rast(GDTlung_s,c('TRDV2','TRDV1','TRDV3','TRDC','TRGC1','TRGV9', "CD3D",'CD3E', 'PTPRC'), ncol =4, assay = 'RNA')
Feature_rast(GDTlung_s, axis.number = T, noaxis = F)+geom_vline(xintercept = c(-5.3,7))
# trim non-T cells 
Feature_rast(GDTlung_s,c('ident','tissue','patient', 'ID'), ncol = 2, noaxis = F, axis.number = T)
GDTlung_s %<>% subset( UMAP_1 < 6 )


Feature_rast(GDTlung_s, 'nCount_RNA', color_grd = 'gdr', facets = 'orig.ident')

Feature_rast(GDTlung_s, 'tissue', facets = 'patient')

Feature_rast(GDTlung_s, 'ITGAE', color_grd = 'threecolor', l =  'green')+
  scale_color_gradient2(low = 'lightgrey', high = 'red', mid = 'purple', midpoint = 3)

# 2nd dim-reduc &clustering -----------------------------------------------

DefaultAssay(GDTlung_s) <- "integrated"


GDTlung_s@meta.data %>%  group_by(v_gene_TRG) %>%  count(patient) 


grep('^TRA', GDTlung_s@assays$integrated@var.features, value = T)

GDTlung_s %<>% ScaleData(vars.to.regress = c("percent.mito",
                                             'orig.ident',
                                             "Phase",
                                             # 'G2M.Score',
                                             # "percent.ribo",
                                             'nCount_RNA' ))

GDTlung_s %<>% NormalizeData(normalization.method = 'LogNormalize',
                             scale.factor = 10000, assay = 'RNA') %>%  ScaleData(assay = 'RNA', 
                         vars.to.regress = c("percent.mito", 
                                             'orig.ident',
                                             "Phase",
                                             # "S.Score", 'G2M.Score',
                                             # "percent.ribo",
                                             'nCount_RNA'))


GDTlung_s %<>%RunPCA( npcs = 100, verbose =
                        T, nfeatures.print = 40)


ElbowPlot(GDTlung_s, ndims = 100)
GDTlung_s %<>% JackStraw(num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 
JackStrawPlot(GDTlung_s, dims = 1:50 )
JackStrawPlot(GDTlung_s, dims = 1:50 ) %T>%
  figsave('GDTlung_8p.jackstraw.pdf' , w = 400, h = 400)

GDTlung_s  %<>%  RunUMAP( dims = 1:35, seed=52,
                          reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:35))

for (i in seq(0.8,1.5,0.1) ) {
  GDTlung_s <- FindClusters(GDTlung_s, resolution = i)
}
Feature_rast(GDTlung_s, "Cell_cluster", noaxis = F, axis.number = T)
GDTlung_s  %<>%  subset(UMAP_1 <5)
GDTlung_s %<>%RunPCA( npcs = 100, verbose =
                        T, nfeatures.print = 40)


ElbowPlot(GDTlung_s, ndims = 100)
GDTlung_s %<>% JackStraw(num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 

for (i in seq(0.8,1.6,0.1) ) {
  GDTlung_s <- FindClusters(GDTlung_s, resolution = i)
}

Feature_rast(GDTlung_s, c("Cell_cluster", "integrated_snn_res.1"))

allres <- Feature_rast(GDTlung_s, paste0("integrated_snn_res.", seq(0.9,1.6,0.1)), ncol =4) %T>% 
  figsave('lunggdt_8p_allres_2023June.pdf', 400,200)
allres


ViolinPlot(GDTlung_s, 'nFeature_RNA', group.by = 'integrated_snn_res.0.8')

table(GDTlung_s$seurat_clusters)


DefaultAssay(GDTlung_s) <- "RNA"

Feature_rast(GDTlung_s,'integrated_snn_res.0.8', facets = 'patient')

# GDTlung_s<- PrepSCTFindMarkers(GDTlung_s,assay = 'SCT')

Feature_rast(GDTlung_s, 'Phase')

Feature_rast(GDTlung_s,c('TRDV2','TRDV1','TRDV3','TRDC','TRGC1','TRGV9', "CD3D",'CD3E'), ncol =4, 
            color_grd = 'grd')
ClusterCompare(GDTlung_s,'10','12', group.by = 'integrated_snn_res.1' , log2fc = 0.5)

# ClusterCompare(GDTlung_s,'6','3', group.by = 'integrated_snn_res.0.9' , log2fc = 0.5 )
# ClusterCompare(GDTlung_s,'1','7', group.by = 'integrated_snn_res.0.7' , log2fc = 0.5)
# 

ClusterCompare(GDTlung_s,'1','8',group.by = 'integrated_snn_res.1' , log2fc = 0.5)


# GDTlung_s$Cell_cluster <- Idents(GDTlung_s) 

GDTlung_s$tissue  %<>% str_replace_all('LN', 'LLN') %>% str_replace_all('Pulm', 'Lung')


GDTlung_s$ID <- paste0(GDTlung_s$tissue,'_',GDTlung_s$patient)

Feature_rast(GDTlung_s,c('TRAC','TRBC1','TRBC2'))


# tissue distribution -----------------------------------------------------






cl_comp <- GDTlung_s@meta.data %>% group_by(tissue, ID, Cell_cluster) %>%  
  summarise(n = n() ) %>% group_by(ID) %>% 
  mutate(percent = n/sum(n)*100) %>% as.data.frame()
# cl_comp$ID <- factor(cl_comp$ID, levels = c(paste0('luLN#',1:),paste0('lung#',1:3) ))
# cl_comp$group <- factor(cl_comp$group, levels = c(, "NCB"))

cl_comp[cl_comp$tissue =='Pulm',]$n <-  -cl_comp[cl_comp$tissue =='Pulm',]$n 
cl_comp[cl_comp$tissue =='LN',]$percent <-  -cl_comp[cl_comp$tissue =='LN',]$percent



cl_comp$percent %>% sum
cl_comp
library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[(2:9)], brewer.pal(9,'Greens')[(2:9)])

Feature_rast(GDTlung_s, "ID", colorset = c(brewer.pal(9,'Blues')[(2:9)], brewer.pal(9,'Greens')[(2:9)]), noaxis = F, sz = 0.5, do.label = F)

cl_comp_flow <- ggplot(cl_comp, aes(y = n, x = Cell_cluster, fill = ID, 
                                    stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,31,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = 'LN\n\nPulm', byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(6,'mm'), axis.text.x = element_text(size = 8, angle = 90),
        legend.text = element_text(size=8),
        axis.line = element_blank())+
  NULL
cl_comp_flow

ggplot(cl_comp, aes(y = n, x = Cell_cluster, fill = ID, 
                    stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,31,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = 'LN\n\nPulm', byrow = T, label.position = 'bottom'),
         color = F )+
  # mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(10,'mm'), axis.text.x = element_text(size = 15, angle = 90),
        legend.text = element_text(size=15),
        axis.line = element_blank())+
  NULL



Umap_cluster<-  Feature_rast(GDTlung_s,c('Cell_cluster'), sz = 0.5, labelsize = 8,
                             noaxis = F)

Umap_cluster_patient<-  Feature_rast(GDTlung_s,c('Cell_cluster'), sz = 0.5, facets = 'patient', facetcol = 3,
                             noaxis = F) %T>% figsave('Umap_cluster_patient.pdf', 200, 270)


UMAP_ID<- Feature_rast(GDTlung_s, 'ID', colorset = alpha(ID_cl, 0.8), do.label = F, noaxis = F, sz = 0.5)+NoLegend()+ggtitle('patient and tissue')

PG(list(UMAP_ID, cl_comp_flow),nrow =1)

Umap_GDTlung<-PG(list(Umap_cluster, PG(list(UMAP_ID, cl_comp_flow),nrow =1)
        ), nrow =1, rw=c(1,1.8), labels = 'AUTO')


figsave(Umap_GDTlung, 'fig1_Umap_GDTlung.pdf', 200, 80)


genes <- c('TRDV1','TRDV2','TRDV3','TRGV9','TRBC1',
           'ZBTB16','TCF7', 'CCR7','TBX21', 'ITGAE','ITGB2',  
           'GZMA','FCGR3A','IFNG',"XCL1", "AREG",
           'KLRB1', 'KLRC1', 'KLRD1','KLRG1',
           "RORC", 'CCR6','DPP4','CXCR3','CXCR6')

gene_UMAP <- Feature_rast(GDTlung_s,genes,  sz = 0.4)
gene_UMAP




library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[(2:9)], brewer.pal(9,'Greens')[(2:9)])

cl_comp_flow <- ggplot(cl_comp, aes(y = n, x = Cell_cluster, fill = ID, 
                                    stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,31,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = 'LN\n\nPulm', byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(6,'mm'), axis.text.x = element_text(size = 8, angle = 90),
        legend.text = element_text(size=8),
        axis.line = element_blank())+
  NULL
cl_comp_flow




# figsave(gene_UMAP2, 'test_geneumap_gdt2020_2.pdf', 400, 400)

# figsave(fig1, 'fig1_lunggdt_6p_UMAP.pdf', 210,270)

saveRDS(GDTlung_s, GDTlung.rds)


# DEGs --------------------------------------------------------------------
DefaultAssay(GDTlung_s) <- "RNA"

hgdTmarkers <- FindAllMarkers(GDTlung_s, test.use = 'bimod', min.pct = 0.05,   only.pos = FALSE )%>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct.dff = pct.1 - pct.2) %>% arrange(cluster, desc(avg_log2FC))

hgdTmarkers


hgdTmarkers%>% arrange( desc(avg_log2FC))
# hgdTmarkers  %<>% arrange(cluster)
write.xlsx(hgdTmarkers, 'lung_gdT_DEG_2025jan.xlsx')

hgdTmarkers

top10_hgdT <- hgdTmarkers  %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(10, avg_log2FC)%>% filter(!grepl('^RP|^MT|^TRA|^TRD|^AC|^HSP|^TRB|-AS1|LINC', gene) & pct.1 > 0.075)

top10_hgdT

Feature_rast(GDTlung_s, "AREG", sz = 0.2, color_grd = "threecolor")

DoHeatmap(subset(GDTlung_s,downsample = 500),features = top10_hgdT$gene,raster = T, 
          group.colors = umap.colors,size = gs(8))%>%heat_theme()
top10_hgdT$gene


cols.use <- list(Cell_cluster = umap.colors, tissue = c("#003399", "#990000"))

DoMultiBarHeatmap(subset(GDTlung_s,downsample = 500),features = top10_hgdT$gene,  group.by='Cell_cluster', size = gs(7),additional.group.sort.by = "tissue", angle = 20,
                  additional.group.by = c("tissue", "patient"),
                  cols.use =cols.use) %>% heat_theme()




Heat_GDTlung <- DoHeatmap(subset(GDTlung_s,downsample = 500),features = top10_hgdT$gene,raster = T, 
                          group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>%
  figsave('heatmap_GDTlung_8p_2023june.pdf',150,300)
Heat_GDTlung


heatDOT_top10_hgdT <-( DotPlot(GDTlung_s,dot.scale = 3.5,
                              features = rev(unique(top10_hgdT$gene)))+mytheme+heattheme+
  
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
    size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom'))) %T>%   figsave('GDTlung_Top10_Dotplot_2023june.pdf',140,250) %T>% print()
 
heatDOT_top10_hgdT



top5_hgdT <- hgdTmarkers %>% filter(!grepl('^RP|^MT|^TR|^AC|^HSP', gene)) %>% 
  arrange(cluster, desc(avg_log2FC))%>%
  group_by(cluster) %>% top_n(5, avg_log2FC)

DoHeatmap(subset(GDTlung_s,downsample = 500),features = top5_hgdT$gene,raster = T, 
          group.colors = umap.colors,size = gs(8))%>%heat_theme() %T>%
  figsave('GDTlung_heatmap_top5_8p_2023.pdf',150,200)


# pseudobulk

pseudobuklGD <- AggregateExpression(GDTlung_s, assays = "RNA", return.seurat = T, group.by = c(  "Cell_cluster", "patient"))

colnames(pseudobuklGD)

pseudobuklGD$orig.ident

bulkDEG_GDT <- FindAllMarkers(pseudobuklGD, test.use = "DESeq2", logfc.threshold = 0.1, only.pos = F)

bulkDEGCD4CD8


bulkDEG_GDT %>%  filter(cluster == "Vg9Vd2-M")  %>%  arrange(desc(avg_log2FC))
  filter(gene %in% c("SELL", "KLRG1", "S1PR1", "KLF2", "RORC", "IL23R", "DPP4", "CCR6"))


top5degGDT_bulk <- bulkDEG_GDT %>%
  arrange(cluster, desc(avg_log2FC))%>%
  filter(!grepl('^RP|^MT|^TR', gene)) %>%
  group_by(cluster) %>% top_n(5, avg_log2FC )



heat_bulk_gdt <- (DoHeatmap(pseudobuklGD,features = top5degGDT_bulk$gene,raster = T,  
                            assay = "RNA", 
                            group.colors = umap.colors,size = gs(8) )  )%>% heat_theme() %T>% print()







Cytokines <- c( 'IFNG', 'GZMA', 'GZMB', 'PRF1',   'AREG','CSF2', 'TNF','CSF1'
               )


Cytokines %>% Feature_density(GDTlung_s, .)
Cytokines %>% Feature_rast(GDTlung_s, ., color_grd = 'grd')



# TFs

TFs <- readLines('HumanTFs_1639.txt') %>% as.vector()
TFs_gd <- intersect(TFs, rownames(GDTlung_s))


Feature_rast(GDTlung_s)

GDTlung_s$pheno

TRMTEMRA_TFs<- ClusterCompare(GDTlung_s, 'gd_TEMRA_LG3', 'gd_TRM_LG6', features = TFs_gd)
TRMTEMRA_TFs$plot

TRMTEMRA_TFs<- ClusterCompare(GDTlung_s, 'Lung_memory', 'Lung_TRM',
                              group.by = "pheno",
                              features = TFs_gd)
TRMTEMRA_Cyts<- ClusterCompare(GDTlung_s, 'gd_TEMRA_P4', 'gd_TRM_P6', features = Cytlist)



GDTlung_s$Cluster_tissue <-  paste0(GDTlung_s$Cell_cluster, "_", GDTlung_s$tissue)


ClusterCompare(GDTlung_s, "Vg9Vd2_M_Lung", "Vg9Vd2_M_LLN", group.by = "Cluster_tissue")

GDTlung_s$tissue

ClusterCompare(GDTlung_s, "gd_TRM_LG6_Lung", "gd_TRM_LG6_LLN", group.by = "Cluster_tissue")

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

Treglist <- c('FOXP3',  'IL2RA', 'IL2RB', 'IKZF2', 'CTLA4', 'TNFRSF18', 'TNFRSF4', 'CAPG', 'CHCHD10', 'GPR83' , 'LGALS3', 'LAG3')


Feature_rast(GDTlung_s, Treglist)

GDTlung_s %<>%  AddModuleScore( features = list(Treglist), name = 'Treg_module', assay = 'RNA')
colnames(GDTlung_s@meta.data) %<>% str_replace("(?<=\\w)1", '')
Feature_density(GDTlung_s,Treglist, ncol = 6)



ViolinPlot(GDTlung_s, 'Treg_module', box = T, colors = umap.colors,x.angle = 90)+

  stat_compare_means()
Feature_rast(GDTlung_s, 'Treg_module', color_grd = 'grd')


GDTlung_s$CD8.Cytotoxictiy
ViolinPlot(GDTlung_s, 'CD8.Cytotoxictiy', box = T, colors = umap.colors,x.angle = 90)+
  
  stat_compare_means()

colnames(GDTlung_s@meta.data)

Feature_rast(GDTlung_s, c(gc_name), color_grd = 'grd') 


plot_density(GDTlung_s, gc_name) 



ViolinPlot(GDTlung_s, gc_name, colors = umap.colors, box = T)
gc_name[9] <- 'Treg_module'
gc_name

GMS  <- GDTlung_s@meta.data[,c(gc_name)] %>% as.data.frame() %>% t()
scaleGM <- scale(t(GDTlung_s@meta.data[,c(gc_name)]))
scaleGMassay <- CreateAssayObject(data = scaleGM)
GDTlung_s@assays$GM <- scaleGMassay
gms <- paste0('GM-', LETTERS[1:8])
gcanno <- as.vector(c(' Naive or immature T cell',
                      "Innate T cell differentiation",
                      'Proliferating',
                      "Type-3 immunity",
                      "Interferon induced",
                      "CTL response (innate)",
                      "CTL response (adaptive)",
                      
                      "Acute activation",
                       'Treg core module'))


c(1,4,7,9) %>% map(~
                       Feature_rast(GDTlung_s, gc_name[.x], mythe = F)+
                       ggtitle(gcanno[.x]) 
) %>% PG(ncol = 4) %T>% figsave('GM_score.pdf', 200, 150)


c(4,7,9) %>% map(~
                   ViolinPlot(GDTlung_s, gc_name[.x],ylabtext = "score",
                                x.angle = 90, mythe = F, size = 10,
                              colors =  umap.colors, box = T)+
                   ggtitle(gcanno[.x])+theme(plot.title = element_text(size = 14))
) %>% PG(ncol = 3) 


ViolinPlot(GDTlung_s, 'GM_D', box = T, colors = umap.colors, mythe = F,x.angle = 90)+ggtitle('Type-3 module')


GM_score_heat <- DoHeatmap(subset(GDTlung_s, downsample = 1000), 
                           raster =T, draw.lines = T, angle = 45,
                           lines.width = 10,group.colors = umap.colors,
                           
                           assay = 'GM', features = gc_name, slot = 'data', size = gs(8)) +hmp2 + mytheme+
  theme(legend.position = 'bottom',
        legend.key.height = unit(2,'mm'))+
  guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))
GM_score_heat

sigtable <- read.xlsx('abt/abd5778_Table_S3.xlsx',sheet = 1)

colnames(sigtable) %<>%   str_remove('.\\(.+\\)')

sigtable %<>% as.list() %>% map(~ na.exclude(.x) %>% as.vector)


for (i in names(sigtable)) {
  GDTlung_s <- AddModuleScore(GDTlung_s, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}

TRM <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "CD69", "GPR25")
GDTlung_s <- AddModuleScore(GDTlung_s, features = list(TRM), name = "TRM", assay = 'RNA')



colnames(GDTlung_s@meta.data) %<>% str_replace("(?<=\\w)1$", '')
GDTlung_s@meta.data

GDTlung_s$Tissue.resident <- GDTlung_s$TRM1
ViolinPlot(GDTlung_s, names(sigtable)[c(1:9, 14,15)], colors = umap.colors, box = T,x.angle = 315, ncol = 4)

# CITESEQ analysis --------------------------------------------------------





GDTlung_s@assays$CITE

GDTlung_s$orig.ident %>% unique



GDTlung_cite <- subset(GDTlung_s, orig.ident %in% c("p71&p31&p27", "p73p77"))

GDTlung_s$orig.ident


Feature_rast(GDTlung_cite)
DefaultAssay(GDTlung_cite) <- 'CITE'

# antibodies <- GDTlung_cite@assays$CITE@counts 
rownames(GDTlung_cite@assays$CITE@counts ) 
abnames <- paste0(c('CD4', 'CD8','CD45RA', 'CD45RO', 'PD1', 'CXCR3', 'CCR6','CD103', 'CD69', 'CCR7',
                    'KLRB1', 'CD27', 'KLRG1', 'IL7R', 'CD26', 'CD49a', 'KLRD1'), '.protein')

Feature_rast(GDTlung_cite, abnames, assay = 'CITE', sz = 0.2, colorgrd = "grd2")







Feature_density(GDTlung_cite, abnames, assay = 'CITE')



trm <- paste0(c('CD103', 'CD49a', 'CD69'), '.protein')


rownames(GDTlung_cite@assays$CITE@counts ) <- abnames
rownames(GDTlung_cite@assays$CITE@data ) <- abnames

rownames(GDTlung_s@assays$CITE@counts ) <- abnames
rownames(GDTlung_s@assays$CITE@data ) <- abnames

saveRDS(GDTlung_s, GDTlung.rds)
saveRDS(GDTlung_cite, 'GDTlung_CITEseqvisual.rds')


DefaultAssay(GDTlung_cite) <- 'CITE'
GDTlung_cite %<>% NormalizeData(normalization.method ='CLR', margin = 2, assay = 'CITE')

GDTlung_cite %<>%ScaleData( assay = 'CITE')


CITEprotein <- Feature_density(GDTlung_cite, abnames, ncol =5, sz = 0.5,assay = 'CITE') %T>% figsave('CITEseq.pdf', 200, 150)

Feature_rast(GDTlung_cite, abnames, ncol =5, sz = 0.2,assay = 'CITE')
Feature_rast(GDTlung_cite, trm,color_grd = 'grd', sz = 0.4, ncol =6, assay = 'CITE')

Feature_rast(GDTlung_s, c('ITGAE', 'ITGA1', 'ZNF683', 'CXCR6', 'SPRY1', 'S1PR5', 'ITGB2', 'KLF2', 'KLF3'))

genenames <- c('CD4', 'CD8A', 'CD8B' , 'PTPRC', 'PDCD1', 'CXCR3', 'CCR6','ITGAE', 'CD69', 'CCR7',
               'KLRB1', 'CD27', 'KLRG1', 'IL7R', 'DPP4', 'ITGA1', 'KLRD1')

CITEgene <- Feature_density(GDTlung_s, genenames, sz = 0.5, ncol =5, assay = 'RNA') 



MemeEff <- paste0(c('CD45RA', 'CD45RO', 'CD27' , 'CCR7', 'KLRG1', 'KLRD1'), '.protein')

Meme <- Feature_rast(GDTlung_cite, MemeEff, sz = 0.5, ncol =3, assay = 'CITE',titlesize = 15) 
Meme


Feature_density(GDTlung_cite, c('CD45RA.protein', 'CD45RO.protein'), ncol = 1, sz =2, titlesize = 15,assay = 'CITE')


Feature_density(GDTlung_s,c('CD27' , 'CCR7', 'KLRG1', 'KLRD1')  )

CITEgene

PG(list(CITEgene, CITEprotein), nrow = 2 ) %T>% figsave('genevsprotein.pdf', 200, 300)

TRmarkers <- Feature_density(GDTlung_cite, c('CD103.protein', 'CD49a.protein'), joint = T, ncol = 3, sz = 0.5)
figsave(TRmarkers, 'TRmarkers.pdf' ,150, 40 )

TRmarkers_vio <- ViolinPlot(GDTlung_cite, c('CD49a.protein', 'CD103.protein'), colors = umap.colors, slot = 'scale.data')
ViolinPlot(GDTlung_cite, c('CD49a.protein', 'CD103.protein'), colors = umap.colors,  box = T)

ViolinPlot(GDTlung_cite, c('CD49a.protein', 'CD103.protein', 'CD26.protein', 'CCR6.protein'), colors = umap.colors)

Feature_rast(GDTlung_cite, d1 ='CD26.protein',d2= 'CCR6.protein', slot = 'scale.data', noaxis = F)

Feature_rast(GDTlung_cite  , d1 ='CD45RA.protein',d2= 'CD27.protein', noaxis = F, axis.number = T, g = "pheno")

Feature_rast(GDTlung_cite  , d1 ='CD45RO.protein',d2= 'CD27.protein', noaxis = F, axis.number = T, g = "pheno")

Feature_rast(GDTlung_cite , d1 ='CD45RA.protein',d2= 'CD27.protein', sz = 5,
             noaxis = F, axis.number = T, slot = 'scale.data')


Feature_rast(GDTlung_cite , g = 'pheno', sz = 0.8, do.label = F,
             d1 ='CD45RA.protein',d2= 'CD45RO.protein', mythe = F,
             noaxis = F, axis.number = T, slot = 'scale.data', othertheme = theme(axis.text = element_text(size = 15)))


Feature_rast(GDTlung_s, 'MKI67')

saveRDS(GDTlung_s, GDTlung.rds)

Feature_rast(GDTlung_cite , g = 'pheno', do.label = F,
             d1 ='CD45RA.protein',d2= 'CD27.protein', noaxis = F, axis.number = T, slot = 'scale.data')



Feature_rast(GDTlung_cite , 
             d1 ='CD45RA.protein',d2= 'CD27.protein', noaxis = F, axis.number = T, slot = 'scale.data')
# ggplot(FetchData(subset(GDTlung_s, downsample = 1000), c('UMAP_1', 'UMAP_2', 'nFeature_RNA')) ,
#         aes(x = UMAP_1, y = UMAP_2, fill =nFeature_RNA))+ 
#   geom_density()+
#   # coord_fixed(ratio = 1)+ 
#   scale_fill_continuous(type = "viridis")

# CD45RACD27 CD103 CD49a --------------------------------------------------------------


(Feature_rast(GDTlung_cite,
              g = 'T_pheno',
              #w facets = c('CD4CD8', 'tissue'),
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 6) %>% rev(),
              
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE', facets = "tissue",
              axis.number = T)+
    geom_hline(yintercept = 0.3, linewidth = 0.2)+
    geom_vline(xintercept = 2.5,linewidth = 0.2)
)


(Feature_rast(GDTlung_cite,
              g = 'tissue',
              #w facets = c('CD4CD8', 'tissue'),
              d1 ='CD45RA.protein', d2 =  'CD27.protein',
              colorset = ggplotColours(12, r = T), sz = 1.5,
              
              # slot = 'scale.data', 
              noaxis = F, assay = 'CITE', facets = "Cell_cluster",
              # do.label = F,
              axis.number = T)+
    geom_hline(yintercept = 0.3, linewidth = 0.2)+
    geom_vline(xintercept = 2.5,linewidth = 0.2)
)


Feature_rast(GDTlung_cite, facets = c("Cell_cluster", 'tissue'), do.label = F, g = "CD103.protein")
Feature_rast(GDTlung_cite, facets = c("Cell_cluster", 'tissue'), do.label = F, g = "CD49a.protein")


ClusterCompare(GDTlung_s %>%  subset(Cell_cluster == "gd_TRM_LG5"), group.by = 'tissue',
               id1 = "LLN", id2 = "Lung")


bc_naive <- colnames(   
  subset(GDTlung_cite,CD27.protein > 0.6 &  CD45RA.protein > 2.5 )
)

bc_cm <- colnames(   
  subset(GDTlung_cite, CD27.protein > 0.6 &  CD45RA.protein <= 2.5 )
)

bc_em <- colnames(   
  subset(GDTlung_cite,CD27.protein <= 0.6 &  CD45RA.protein <= 2.5 )
)

bc_tmra <- colnames(   
  subset(GDTlung_cite,CD27.protein <= 0.6 &  CD45RA.protein > 2.5 )
)


GDTlung_cite@meta.data  %<>%  mutate(T_pheno= case_when( 
  bc_backup %in% bc_naive ~ 'naive',
  bc_backup %in% bc_cm ~ 'Tcm',
  bc_backup %in% bc_em ~ 'Tem',
  bc_backup %in% bc_tmra ~ 'Temra'
  
)   )


Feature_rast(GDTlung_cite, c("T_pheno", "ident"))


Feature_rast(GDTlung_cite,
             g = 'T_pheno',
             facets = c('patient'),
             d1 ='CD45RA.protein', d2 =  'CD27.protein',
             colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 4) %>% rev(),
             
             # slot = 'scale.data',
             noaxis = F, assay = 'CITE',
             # slot = "scale.data",
             axis.number = T)


Feature_rast(GDTlung_cite, c("T_pheno"), facets = "tissue")




(Feature_rast(GDTlung_cite,
              # g = 'pheno',
              facets = c( 'tissue'),
              d1 ='CD49a.protein', d2 =  'CD103.protein',
              # colorset = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 6) %>% rev(),
              
              # slot = 'scale.data',
              noaxis = F, assay = 'CITE',
              axis.number = T)+
    geom_hline(yintercept = 0.5, linewidth = 0.2)+
    geom_vline(xintercept = 1,linewidth = 0.2)
)


bc_Trm <- colnames(   
  subset(GDTlung_cite,CD49a.protein > 2 &  CD103.protein > 1 )
)


GDTlung_cite@meta.data  %<>%  mutate(Trm= case_when( 
  bc_backup %in% bc_Trm ~ 'Trm',
  !(bc_backup %in% bc_Trm) ~ T_pheno
  
)   )




Feature_rast(GDTlung_cite, "Trm")

table(GDTlung_cite$Trm, GDTlung_cite$Cell_cluster)


I wa# adjust cluster  ----------------------------------------------------------
Idents(GDTlung_s) <- GDTlung_s$integrated_snn_res.1
Feature_rast(GDTlung_s)

pheno <- data.frame(integrated_snn_res.1 = as.character(0:15), 
                    Cell_cluster = c("gd_TEMRA_P1", "gd_TEM_L2", "gd_TRM_P5", "gd_TRM_P6",
                                     "Vg9Vd2_M","gd_TEMRA_P2", "gd_TEMRA_L4",  "gd_TEMRA_P4", 
                                     "gd_Naive_L1", "gd_TEMRA_P3", "gd_TRM_P7", "gd_TEM_L3",
                                     "gd_Naive_L1","gd_TEMRA_L4","gd_TRM_P5","gd_TEMRA_P1"
                                    )
                    
                    
                    )


clod <- c("gd_TEMRA_P1", "gd_TEM_L2", "gd_TRM_P5", "gd_TRM_P6",
          "Vg9Vd2_M","gd_TEMRA_P2", "gd_TEMRA_L4",  "gd_TEMRA_P4", 
          "gd_Naive_L1", "gd_TEMRA_P3", "gd_TRM_P7", "gd_TEM_L3",
          "gd_Naive_L1","gd_TEMRA_L4","gd_TRM_P5","gd_TEMRA_P1"
)%>% unique() %>% sort() %T>% print() 


Feature_rast(GDTlung_s)
GDTlung_s@meta.data 

GDTlung_s@meta.data %<>% left_join(pheno, by = 'integrated_snn_res.1', suffix = c('', '')) %>% 
 mutate( pheno = case_when(
   grepl("Vg9", Cell_cluster) ~ "Vg9Vg2",
   grepl("Naive", Cell_cluster) ~ "Naive",
   grepl("L", Cell_cluster) ~ "LN_memory",
   
   
   grepl("TRM", Cell_cluster) ~ "Lung_TRM",
   grepl("TEMRA_P", Cell_cluster) ~ "Lung_memory" )) %>% 
  
    `rownames<-`(GDTlung_s$bc_backup)


GDTlung_s$Cell_cluster  %<>% factor(levels = clod)


# Feature_rast(GDTlung_s, c('pheno','tissue'))
Feature_rast(GDTlung_s, c('pheno','Cell_cluster'))

Feature_rast(GDTlung_s, c('pheno','tissue'))


Feature_rast(GDTlung_s, c('Trm','Cell_cluster'))

# # change L to LN and P to Lu, pulm to Lung , LN to LLN
GDTlung_s$Cell_cluster  %<>% str_replace("_L", "_LN")
GDTlung_s$Cell_cluster  %<>% str_replace("_P", "_Lu")
GDTlung_s$Cell_cluster  %<>% str_replace("_Lu", "_LG")

GDTlung_s$Cell_cluster  %<>% str_replace("_LG3", "_LG4")

GDTlung_s$Cell_cluster  %<>% str_replace("_LG4", "_LG3")
GDTlung_s$Cell_cluster  %<>% str_replace("_L2G4", "_LG4")

GDTlung_s$Cell_cluster  %<>% str_replace("TEM_LN2", "TCM_LN2") %>% 
  str_replace("TEM_LN3", "TEMRA_LN3")


GDTlung_s$Cell_cluster_old <-  GDTlung_s$Cell_cluster

GDTlung_s$Cell_cluster  %<>% str_replace("LN3", "_Mix")%>% 
  str_replace("Vg9Vd2_M", "Vg9Vd2_Mix") %>% 
  str_replace("_LN4", "_LN1")%>% 
  str_replace("Naive_LN1", "Naive_LN3") 
GDTlung_s$Cell_cluster %<>% 
  factor(levels = c("Vg9Vd2_Mix","gd_TEMRA__Mix" , "gd_TEMRA_LN1","gd_TCM_LN2",
                    "gd_Naive_LN3","gd_TEMRA_LG1",
                    "gd_TEMRA_LG2",  "gd_TEMRA_LG3" , "gd_TEMRA_LG4", 
                    "gd_TRM_LG5" ,   "gd_TRM_LG6"   , "gd_TRM_LG7"    ))



Idents(GDTlung_s) <- GDTlung_s$Cell_cluster



Feature_rast(GDTlung_s)



saveRDS(GDTlung_s, GDTlung.rds)






# tissue distribution -----------------------------------------------------

ID_cl <- c(brewer.pal(9,'Blues')[3:9], brewer.pal(9,'Reds')[3:9])




cl_comp_gd <- GDTlung_s@meta.data %>% filter(patient != "p31") %>% 
  mutate(ID  = paste0(tissue, "_", patient)) %>% 
 group_by(tissue, ID, Cell_cluster) %>%  
  summarise(n = n() ) %>% group_by(ID) %>% 
  mutate(percent = n/sum(n)*100) %>% as.data.frame()


cl_comp_gd[cl_comp_gd$tissue =='Lung',]$n <-  -cl_comp_gd[cl_comp_gd$tissue =='Lung',]$n 
cl_comp_gd[cl_comp_gd$tissue =='LLN',]$percent <-  -cl_comp_gd[cl_comp_gd$tissue =='LLN',]$percent



cl_comp_gd$percent %>% sum
cl_comp_gd
library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[3:9], brewer.pal(9,'Reds')[3:9])

cl_comp_flow_gd <- ggplot(cl_comp_gd, aes(y = n, x = Cell_cluster, fill = ID, 
                                    stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = 'LN\n\nPulm', byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(6,'mm'), axis.text.x = element_text(size = 8, angle = 90),
        legend.text = element_text(size=8),
        axis.line = element_blank())+
  NULL
figsave(cl_comp_flow, "GDTlung_tissue_distribution.pdf", 70,70)






Umap_cluster<-  Feature_rast(GDTlung_s,c('Cell_cluster'), sz = 0.5, labelsize = 8,
                             noaxis = F)

Umap_cluster_patient<-  Feature_rast(GDTlung_s,c('Cell_cluster'), sz = 0.5, facets = 'patient', facetcol = 3,
                                     noaxis = F) %T>% figsave('Umap_cluster_patient.pdf', 200, 270)


   




PG(list(UMAP_ID, cl_comp_flow),nrow =1)

Umap_GDTlung<-PG(list(Umap_cluster, PG(list(UMAP_ID, cl_comp_flow),nrow =1)
), nrow =1, rw=c(1,1.8), labels = 'AUTO')


figsave(Umap_GDTlung, 'fig1_Umap_GDTlung.pdf', 200, 80)

zx
 genes <- c('TRDV1','TRDV2','TRDV3','TRGV9','TRBC1',
           'ZBTB16','TCF7', 'CCR7','TBX21', 'ITGAE','ITGB2',  
           'GZMA','FCGR3A','IFNG',"XCL1", "AREG",
           'KLRB1', 'KLRC1', 'KLRD1','KLRG1',
           "RORC", 'CCR6','DPP4','CXCR3','CXCR6')

gene_UMAP <- Feature_rast(GDTlung_s,genes,  sz = 0.4)
gene_UMAP







cl_comp <- GDTlung_s@meta.data %>% group_by(tissue, ID, integrated_snn_res.1.1) %>%  
  summarise(n = n() ) %>% group_by(ID) %>% 
  mutate(percent = n/sum(n)*100) %>% as.data.frame() %>% rename( Cell_cluster=integrated_snn_res.1.1)


cl_comp[cl_comp$tissue =='Pulm',]$n <-  -cl_comp[cl_comp$tissue =='Pulm',]$n 
cl_comp[cl_comp$tissue =='LN',]$percent <-  -cl_comp[cl_comp$tissue =='LN',]$percent



cl_comp$percent %>% sum
cl_comp
library(ggalluvial)
library(RColorBrewer)

ID_cl <- c(brewer.pal(9,'Blues')[(2:9)], brewer.pal(9,'Greens')[(2:9)])

cl_comp_flow <- ggplot(cl_comp, aes(y = n, x = Cell_cluster, fill = ID, 
                                    stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,31,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell number')+
  xlab('cluster')+ 
  xlab(NULL)+
  scale_y_continuous(labels = abs)+
  geom_hline(yintercept = 0,size = 0.5)+
  guides(fill = guide_legend(nrow = 2, title = 'LN\n\nPulm', byrow = T, label.position = 'bottom'),
         color = F )+
  mytheme+
  theme(legend.position = 'bottom', legend.key.height  = unit(2, 'mm'), legend.key.width   = unit(6,'mm'), axis.text.x = element_text(size = 8, angle = 90),
        legend.text = element_text(size=8),
        axis.line = element_blank())+
  NULL
cl_comp_flow




# figsave(gene_UMAP2, 'test_geneumap_gdt2020_2.pdf', 400, 400)

# figsave(fig1, 'fig1_lunggdt_6p_UMAP.pdf', 210,270)

saveRDS(GDTlung_s, GDTlung.rds)





# TCR  --------------------------------------------------------------------

# here we have TCR data from two different datasets: 1 and 2_3
# and we will combine them as one
# for each sample  we have two sequence run.

TCRdirs <- list.dirs(path = 'raw') %>% str_subset('VDJ.+outs$' )
TCRdirs
GDTlung_s$bc_backup <- rownames(GDTlung_s@meta.data)
GDTlung_s@meta.data %>% group_by(orig.ident) %>% slice(1) %>% select(orig.ident, bc_backup)





lungTCR <- map2(TCRdirs, c(2,1,3,4), ~ 
                  (
                    do.call(rbind, lapply( list.files(path = .x, full.names = T,
                                                      pattern = 'filtered_contig_annotations')  , read.csv))
                   
                 ) %>% 
                  filter(productive == 'True'& is_cell == 'True'& grepl('GV|DV', v_gene) & 
                           high_confidence  == 'True' & 
                           # filter out unproductive TRG genes 
                           !(v_gene %in% c('TRGV1' , 'TRGV5P', 'TRGV6', 'TRGV7', 'TRGV10', 'TRGV11'))) %>%
                  dplyr::select(c(1, 5:10, 13,14))  %>%
                  mutate(  
                    bc_backup = paste0(barcode, "_",.y)   ) %>% dplyr::select(-barcode)) %>% reduce(.f = rbind)


lungTCR$v_gene %>% unique





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
dim(TCRs)
dim(GDTlung_s)

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

GDTlung_s@meta.data  %<>%  select_at(.vars = vars(-contains(c('TRG','TRD','paired', '.x', '.y'))))
GDTlung_s

GDTlung_s@meta.data  %<>% left_join(TCRs, by = 'bc_backup', suffix = c('','')) %>% `rownames<-`(GDTlung_s@meta.data$bc_backup )

GDTlung_s@meta.data

Feature_rast(GDTlung_s, 'v_gene_TRG')
Feature_rast(GDTlung_s, 'v_gene_TRD')

# # this line to remove repeated join
# GDTlung_s@meta.data %<>% `colnames<-`(str_remove(colnames(GDTlung_s@meta.data), '.y')) %>%
#   select(!ends_with('.x')) 
# mapped TCR
TRDmapp <- GDTlung_s@meta.data %>% count(ID, v_gene_TRD) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 
TRGmapp <- GDTlung_s@meta.data %>% count(ID,v_gene_TRG) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 


pairedmapp <-GDTlung_s@meta.data %>% count(ID,paired_sp) %>%  group_by(ID)  %>%
  mutate(percent = n/sum(n)*100) 


GDTlung_s$paired_sp %>% table
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
                
                    mytheme+
                    theme(axis.text.x = element_text(angle = 270))+
                    gglp('b')
                  
) %>% PG(nrow =1, labels = 'AUTO')
mappedTCR
# ViolinPlot(subset(GDTlung, length_TRD >0 & patient == '#3'),'length_TRD')

saveRDS(GDTlung_s, GDTlung.rds)

TCRmapping <- GDTlung_s@meta.data %>% group_by(tissue, patient) %>%  count(TCR_summary) %>% 
  mutate(percent = n/sum(n)*100)


GDTlung_s@meta.data  %<>%  mutate(Vg9Vd2 = case_when(paired_sp == 'GV9 DV2' ~ 'GV9 DV2',
                                                     !is.na(paired_sp)  ~ 'other γδTCR'
)) %>% `row.names<-`(GDTlung_s$bc_backup)


Feature_rast(GDTlung_s, 'Vg9Vd2')




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



patient <- c(GDTlung_s$patient %>% unique()) %>% sort()

TCRD_FREQ <-map(patient,~
                  Feature_rast(subset(GDTlung_s, patient == .x),'cdr3_TRD_perc',color_grd = 'grd', navalue = 'transparent')+
                  ggtitle(paste('patient',.x,'\nTRD_clone (%)'))
) %>% PG(nrow =2)


map(patient,~
       Feature_rast(subset(GDTlung_s, patient == .x),'cdr3_paired_perc',color_grd = 'grd', navalue = 'transparent')+
       ggtitle(paste('patient',.x,'\nTRD_clone (%)'))
) %>% PG(nrow =2)


map(patient,~
      Feature_rast(subset(GDTlung_s, patient == .x),'cdr3_TRG_perc',color_grd = 'grd', navalue = 'transparent')+
      ggtitle(paste('patient',.x,'\nTRD_clone (%)'))
) %>% PG(nrow =2)
TCRD_FREQ


Feature_rast(GDTlung_s, 'cdr3_TRD_perc', color_grd = 'grd', mythe = F, noaxis = T, navalue = 'transparent')+
  ggtitle('TRD_clone (%)')


GDTlung_s@meta.data %<>% 
  mutate(clonal_expansion =case_when(cdr3_TRD_freq == 1 ~ 'monoclonal',
                                      nr(cdr3_TRD_freq, 2,9)~'low (2~9)',
                                         nr(cdr3_TRD_freq, 10,20)~'moderate (10~20)',
                                           cdr3_TRD_freq >20 ~ 'high (>20)' ),
     clonal_expansion = factor(clonal_expansion, levels = c('high (>20)', 'moderate (10~20)',
                                                                                    'low (2~9)', 'monoclonal'))   )


GDTlung_s@meta.data %<>% mutate(clonal_expansion_TCRGD =case_when(cdr3_paired_freq == 1 ~ 'monoclonal',
                                                            nr(cdr3_paired_freq, 2,4)~'low (2~4)',
                                                            nr(cdr3_paired_freq, 5,9)~'moderate (5~9)',
                                                            
                                                            cdr3_paired_freq >9 ~ 'high (>9)' ),
                                clonal_expansion = factor(clonal_expansion, levels = c('high (>9)', 'moderate (5~9)',
                                                                                       'low (2~4)', 'monoclonal'))   )

Feature_rast(GDTlung_s, "clonal_expansion", facets = "patient", do.label = F,
            colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'), navalue = "transparent")


Feature_rast(GDTlung_s, c('SELL', 'CD27', 'CD28', 'IL7R', 'TCF7'))

Feature_rast(GDTlung_s, "clonal_expansion",  do.label = F,
             colorset =  c('#DC143C','#9400D3', '#1E90FF', '#FAFAD2'), navalue = "transparent")



Feature_rast(GDTlung_s, "cdr3_TRD_perc")

# analysis and make figures -----------------------------------------------



# TRD on umap
TCRumap <- Feature_rast(GDTlung_s,c('Cell_cluster','v_gene_TRD','paired_sp'),do.label = F, sz =1, labels =LETTERS[4:6]   )
TCRumap
Feature_rast(GDTlung_s,'ident', c('ID'))+facet_grid(~ID)


GDTlung_s$length_TRD

ViolinPlot(GDTlung_s %>% subset(v_gene_TRD == 'TRDV1') , 'length_TRD', box = T)+ ylim(350,700)


Feature_rast(GDTlung, facet = 'patient')+facet_grid(~patient)

# top expanded clones -----------------------------------------------------


# top TRD clones
Top5TRD <- GDTlung_s@meta.data %>% filter( !is.na(cdr3_TRD) ) %>%
  group_by(patient,v_gene_TRD, cdr3_TRD) %>%
  dplyr::count(patient,cdr3_TRD) %>% arrange(desc(n)) %>%
  group_by(patient) %>%  dplyr::slice(1:5) %>% filter(n>=5) %>%
  dplyr::mutate(Top5TRD = paste0(patient, ' ',v_gene_TRD," ",cdr3_TRD,' ',  n)) %>% 
  arrange(patient,desc(n)) %>%
  select(patient, cdr3_TRD, Top5TRD) %>%
  ungroup()



Top5TRD_Pulm <- GDTlung_s@meta.data %>% filter( !is.na(cdr3_TRD) & tissue == 'Pulm') %>%
  group_by(patient,v_gene_TRD, cdr3_TRD) %>%
  dplyr::count(patient,cdr3_TRD) %>% arrange(desc(n)) %>%
  group_by(patient) %>%  dplyr::slice(1:5) %>% filter(n>=5) %>%
  dplyr::mutate(Top5TRD_Pulm = paste0(patient, ' ',v_gene_TRD," ",cdr3_TRD,' ',  n)) %>% 
  arrange(patient,desc(n)) %>%
  select(patient, cdr3_TRD, Top5TRD_Pulm) %>%
  ungroup()




Top5TRD_LN <- GDTlung_s@meta.data %>% filter( !is.na(cdr3_TRD) & tissue == 'LN') %>%
  group_by(patient,v_gene_TRD, cdr3_TRD) %>%
  dplyr::count(patient,cdr3_TRD) %>% arrange(desc(n)) %>%
  group_by(patient) %>%  dplyr::slice(1:5) %>% filter(n>=5) %>%
  dplyr::mutate(Top5TRD_LN = paste0(patient, ' ',v_gene_TRD," ",cdr3_TRD,' ',  n)) %>% 
  arrange(patient,desc(n)) %>%
  select(patient, cdr3_TRD, Top5TRD_LN) %>%
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

GDTlung_s$Top5TRD_Pulm <- NULL
GDTlung_s$Top5TRD_LN <- NULL
GDTlung_s@meta.data %<>% 
  left_join(Top5TRD_Pulm,by = c('cdr3_TRD','patient') )%>%
  left_join(Top5TRD_LN,by = c('cdr3_TRD','patient') )%>%
  left_join(Top5TRD,by = c('cdr3_TRD','patient') )%>%
  mutate(Top5TRD = case_when( !is.na(Top5TRD) ~ Top5TRD,
                                   is.na(Top5TRD) & !is.na(cdr3_TRD)  ~ 'mapped TCR')) %>%
  mutate(Top5TRD_Pulm = case_when( !is.na(Top5TRD_Pulm) ~ Top5TRD_Pulm,
                                   is.na(Top5TRD_Pulm) & !is.na(cdr3_TRD)  ~ 'mapped TCR')) %>%
  mutate(Top5TRD_LN = case_when( !is.na(Top5TRD_LN) ~ Top5TRD_LN,
                                   is.na(Top5TRD_LN) & !is.na(cdr3_TRD)  ~ 'mapped TCR')) %>%
  `rownames<-`(GDTlung_s$bc_backup)
GDTlung_s$Top5TRD_Pulm %<>% factor(levels = GDTlung_s$Top5TRD_Pulm %>%unique() %>%  str_sort())
GDTlung_s$Top5TRD_LN %<>% factor(levels = GDTlung_s$Top5TRD_LN %>%unique() %>%  str_sort())
GDTlung_s$Top5TRD %<>% factor(levels = GDTlung_s$Top5TRD %>%unique() %>%  str_sort())


Feature_rast(GDTlung_s)

# GDTlung@meta.data %<>% `colnames<-`(str_remove(colnames(GDTlung@meta.data), '.y')) %>% select(!ends_with('.x')) 

Feature_rast(GDTlung_s, c('tissue', 'Top5TRD_Pulm','Top5TRD_LN'), do.label = F, ncol = 1,
             colorset = ggplotColours(60))

top10trdcl <-
  c('#993333','#ff3300','#ff0066','#cc0099','#cc3300','#391313','#ff6666','#ffa64d','#ff8080','#660066',
    '#000066','#0000ff', '#336699','#006666','#6666ff', '#004d3d','#00b3b3','#003399','#39ac73','#1ab2ff',
    '#ffffcc'
  )

luLN_expanded<-Feature_rast(GDTlung_s, "Top5TRD_LN", do.label = F, facets ='patient')+
  scale_color_manual(values = c(  
    top10trdcl[1:(length(unique(GDTlung_s$Top5TRD_LN))-2)],
    'mapped TCR' = alpha('yellow',0.5)), 
    na.value =alpha('lightgrey',0.5) )+
  ggtitle("Most Expanded TRD in luLN")




map(unique(GDTlung_s$patient), ~ Feature_rast(subset(GDTlung_s, patient == .x),  'Top5TRD_LN', do.label = F, sz = 0.3 , navalue = 'transparent')+
      ggtitle(.x )+ 
      guides(color = guide_legend(ncol =1, override.aes = list(size = 1.5))) ) %>% PG(ncol = 1)


map(unique(GDTlung_s$patient), ~ Feature_rast(subset(GDTlung_s, patient == .x,),  'Top5TRD_Pulm', do.label = F, sz = 0.3 , navalue = 'transparent')+
      ggtitle(.x )+ 
      guides(color = guide_legend(ncol =1, override.aes = list(size = 1.5))) ) %>% PG(ncol = 1)


map(unique(GDTlung_s$patient), ~ Feature_rast(subset(GDTlung_s, patient == .x,),  'Top5TRD', do.label = F, sz = 1 , navalue = alpha('lightgrey',0.5) , colorset = c('lightyellow', ggplotColours(5))  )+
      ggtitle(.x )+ 
      guides(color = guide_legend(ncol =1, override.aes = list(size = 1.5))) ) %>% PG(ncol = 2)

Feature_rast(GDTlung_s, 'v_gene_TRD', do.label = F, mythe = F, othertheme = theme(legend.text = element_text(size = 15)))+
  guides(color= guide_legend(override.aes = list(size =3)))

saveRDS(GDTlung_s, GDTlung.rds)


# gini index ---------------------------------------------------------

library(reldist)

GDTTRDfreq <-  GDTlung_s@meta.data%>% 
  filter(!is.na(cdr3_TRD)) %>% 
  dplyr::count( patient, Cell_cluster, cdr3_TRD) %>%
  arrange(desc(n))

GDTTRDfreq



gini_TRD <-  GDTTRDfreq %>% 
  dplyr::group_by( patient, Cell_cluster) %>%
  summarise(Gini_Index = gini(n), sum = sum(n))  %>%
  mutate(Gini_Index = replace(Gini_Index, sum < 20, NA))%>% 
  ungroup() %>% 
  tidyr::complete(patient, Cell_cluster , 
                  fill = list(Gini_Index = NA, sum = NA))








giniindex_TRD <- ggplot(gini_TRD, aes(x = Cell_cluster , y = Gini_Index, 
                                           color = patient, group= Cell_cluster))+
  geom_boxplot()+geom_point()+
  color_m(color = set_sample(umap.colors, s = 22))+theme_minimal()  +
  theme(axis.text.x = element_text(angle = 90))






# venny nonVD2  TRD  ----------------------------------------------------

GDTlung_s$TRD
Feature_rast(GDTlung_s)
# GDTlung_s@meta.data %<>% mutate(big_cluster = case_when(Cell_cluster %in% c('c1','c2','c8', 'c9')~ 'nonVD2_lung_1',
#                                                         Cell_cluster %in% c('c5','c7' ,'12')~ 'nonVD2_lung_2',
#                                                         Cell_cluster %in% c('c3', 'c6', 'c10')~ 'nonVD2_luLN_1',
#                                                         Cell_cluster %in% c('c11')~ 'nonVD2_luLN_2',
#                                                         Cell_cluster %in% c('c4')~ 'VD2')
# )%>%
#   `rownames<-`(GDTlung_s$bc_backup)
# GDTlung_s$big_cluster

# UMAP_bigcluster <- Feature_rast(GDTlung_s, 'big_cluster')

# ClusterCompare(GDTlung_s,'nonVD2_lung_1','nonVD2_lung_2', group.by = 'big_cluster', log2fc = 1)

# patient
# 
# TRD_uniq <- c()
# for (i in patient) {
#   TRD_uniq[[i]] <-    map(c('P_Trm', 'P_Tcirc', 'L_Tcirc'), ~
#                             GDTlung_s@meta.data %>%
#                             dplyr::filter( pheno==.x  & !is.na(cdr3_TRD) &patient ==i &cdr3_TRD_freq>1 ) %>%
#                             pull(cdr3_TRD) %>% unique ) %>%
#     setNames(c('P_Trm', 'P_Tcirc', 'L_Tcirc'))
# }
# 
# 
# TRD_uniq$p25
# 
# 
# library(eulerr)
# library(ggplotify)
# 
# set.seed(8)
# # euler(TRD_uniq$`#3`, shape = 'ellipse')
# 
# set.seed(8)
# Venny_P_L <- map2(TRD_uniq, list(F,F,F,list(cex = .66),F,F,F,list(cex = .66)),~
#                     plot(euler(.x, shape = 'ellipse'),
#                             quantities = list(cex = 0.66),  prob = T,labels = NULL,legend =.y,
#                               fill = alpha(umap.colors,0.8),
# ) %>% ggplotify::as.ggplot() + ggtitle(' unique CDR3')+
#   theme(plot.title = element_text(size = 8, hjust = 0.5))) %>%
#   PG(nrow = 2, rw = c(1,1,1.8), labels = patient, label_fontface = 'plain',label_size = 8 )

# venny VD2  TRD  ----------------------------------------------------
GDTlung_s$tissue %>% unique()
patient <-  GDTlung_s$paired %>% unique()
TRD_VD2_uniq <- c()
for (i in patient) {
  TRD_VD2_uniq[[i]] <-    map(c('Lung','LLN'), function(x) GDTlung_s@meta.data %>% 
        dplyr::filter( Cell_cluster== 'Vg9Vd2_M' &
                                                v_gene_TRD == 'TRDV2' & cdr3_TRD_freq >=1 & patient ==i & tissue == x) %>% pull(cdr3_TRD) %>% unique ) %>% 
    setNames(c('Lung','LLN'))
}


TRD_VD2_uniq
set.seed(8)
map2(TRD_VD2_uniq, list(F,F,F,list(cex = .66),F,F,F,list(cex = .66)),~ plot(euler(.x, shape = 'ellipse'), 
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


TCRfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_paired_freq >1  & paired != 'GV9 DV2'    )  %>% group_by(pheno,tissue, cdr3_paired) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq) #%>% mutate(cdr3_paired=factor(cdr3_paired, levels = unique(cdr3_paired)))
TCRfreqtable

Feature_rast(GDTlung_s,"pheno")

gdTCRsharing <- TCRfreqtable %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 

  filter(pt %in% c('Lung_TRM_Lung', 'LN_memory_LLN', 'Lung_memory_Lung')) %>%

  mutate(pheno = factor(pheno, levels = c('Lung_TRM', 'LN_memory', 'Lung_memory'))) %>% 
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_paired, 
         stratum= cdr3_paired, alluvium  = cdr3_paired))+
  ggtitle("paired gd TCR sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = ggplotColours(211))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 315),
        legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  # facet_wrap(~patient)+
  NULL
gdTCRsharing



gdTCRsharing
gdTCRsharing %>% figsave("GD_pairedTCR_sharing.pdf", 100,100)
TRDfreqtable
TRDfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_TRD_freq >0 )  %>% group_by(pheno,tissue, cdr3_TRD) %>%  summarise(TCRfreq = n()) %>% ungroup() %>% 
  arrange(TCRfreq) %>% mutate(cdr3_TRD=factor(cdr3_TRD, levels = unique(cdr3_TRD)))

TRDsharing <- TRDfreqtable %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  
  filter(pt %in% c('Lung_TRM_Lung', 'LN_memory_LLN', 'Lung_memory_Lung')) %>%  
  
  mutate(pheno = factor(pheno, levels = c('Lung_TRM','Lung_memory',  'LN_memory'))) %>% 
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_TRD, 
         stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("TRD sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = ggplotColours(211))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  # facet_wrap(~patient)+
  NULL
TRDsharing
TRDsharing %>% figsave("GD_TRD_sharing.pdf", 100,100)


TRDfreqtable_NT <- GDTlung_s@meta.data %>% filter(cdr3_TRD_freq >0 )  %>% group_by(pheno,tissue, cdr3_nt_TRD) %>%  summarise(TCRfreq = n()) %>% ungroup() %>% 
  arrange(TCRfreq) %>% mutate(cdr3_nt_TRD=factor(cdr3_nt_TRD, levels = unique(cdr3_nt_TRD)))

TRDsharing_NT <- TRDfreqtable_NT %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  
  filter(pt %in% c('Lung_TRM_Lung', 'LN_memory_LLN', 'Lung_memory_Lung')) %>%  
  
  mutate(pheno = factor(pheno, levels = c('Lung_TRM','Lung_memory',  'LN_memory'))) %>% 
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_nt_TRD, 
         stratum= cdr3_nt_TRD, alluvium  = cdr3_nt_TRD))+
  ggtitle("TRD sharing")+
  geom_flow(stat = "alluvium",size = 0.2,
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = ggplotColours(211))+
  geom_stratum(size = 0.2)+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  # facet_wrap(~patient)+
  NULL
TRDsharing_NT
TRDsharing %>% figsave("GD_TRD_sharing.pdf", 100,100)




TCRfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_paired_freq >1 )  %>% group_by(pheno, cdr3_paired, patient, tissue) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq)


TCRfreqtable
TCRfreqtable %>%  mutate(pt = paste0(pheno, "_", tissue)) %>% 
  
  filter(pt %in% c('Lung_TRM_Pulm', 'LN_memory_LN', 'Lung_memory_Pulm')) %>% 
  mutate(pheno = factor(pheno, levels = c('Lung_TRM', 'LN_memory', 'Lung_memory')))%>%
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  ggtitle("TCR sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  scale_fill_manual(values = rainbow(300))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  facet_wrap(~patient, scales = 'free_y')+
  NULL


TCRfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_paired_freq >=1 & patient != 'p77')  %>% group_by(Cell_cluster,
                                                                                    cdr3_TRD) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq)



TCRfreqtable
library(tidyr)
crdf <- TCRfreqtable %>% filter(big_cluster %in% c( 'nonVD2_lung_1','nonVD2_lung_2', 'nonVD2_luLN_1')) %>%select(-patient) %>% 
  mutate(big_cluster = factor(big_cluster, levels = c('nonVD2_lung_1','nonVD2_luLN_1','nonVD2_lung_2'))) %>%
  spread(key = big_cluster, value = TCRfreq) %>% as.data.frame()

crdf[is.na(crdf)] <- 0
rownames(crdf) <- crdf$cdr3_paired
crdf$cdr3_paired<- NULL

chordDiagram(crdf, reduce = 0, gap.degree = 158)

TCRDfreqtable <- GDTlung_s@meta.data %>% filter(cdr3_TRD_freq >0 )  %>% group_by(big_cluster, cdr3_TRD, patient) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq)
TCRDfreqtable %>% filter(big_cluster %in% c( 'nonVD2_lung_1','nonVD2_lung_2', 'nonVD2_luLN_1')) %>%  
  mutate(big_cluster = factor(big_cluster, levels = c('nonVD2_lung_1','nonVD2_lung_2', 'nonVD2_luLN_1'))) %>%
  ggplot(
    aes( x = big_cluster, y = TCRfreq, fill = cdr3_TRD,  stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("TCRD sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  scale_fill_manual(values = rainbow(534))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle = 90))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  facet_wrap(~patient)+
  NULL



# VD2
GDTlung_s$paired
VD2freqtable <- GDTlung_s@meta.data %>%
  dplyr::filter(cdr3_paired_freq > 1 &   paired == 'GV9 DV2' )  %>%
  group_by(tissue, cdr3_paired) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
VD2freqtable %>% 
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  ggtitle("Vg9Vd2 TCR  sharing between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq
# none(nonVg9vd2)

VD2freqtable <- GDTlung_s@meta.data %>%
  dplyr::filter(cdr3_TRD_freq > 1 &   v_gene_TRD == 'TRDV2' & Cell_cluster == "Vg9Vd2_M" )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
VD2freqtable %>% 
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD  , alluvium  = cdr3_TRD            ))+
  ggtitle("Vg9Vd2 TCR  sharing between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq






NonVD2freqtable <- GDTlung_s@meta.data %>% dplyr::filter(cdr3_paired_freq >1 & 
                                                        paired != 'GV9 DV2' )  %>%
  group_by(tissue, cdr3_paired) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
NonVD2freqtable %>% 
  
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_paired,  stratum= cdr3_paired, alluvium  = cdr3_paired))+
  ggtitle("Non Vg9Vd2 TCR sharing in between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq



NonVD2freqtable_LN3 <- GDTlung_s@meta.data %>% dplyr::filter(cdr3_TRD_freq >1 & 
                                                           paired != 'GV9 DV2' & Cell_cluster == "gd_TRM_LG6" )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
NonVD2freqtable_LN3 %>% 
  
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("Non Vg9Vd2 TCR sharing in between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq




L1L2freqtable <- GDTlung_s@meta.data %>% dplyr::filter(cdr3_TRD_freq >1 & 
                                                         patient != 'p77'&
                                                        Cell_cluster %in% c('P1')  )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq)


ggplot(L1L2freqtable,
  aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("TCRD sharing in between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  scale_fill_manual(values = rainbow(310))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
GDTlung_s$cdr3_TRD_freq












# TCR similarity ----------------------------------------------------------
# library(scRepertoire)
library(immunarch)


GDTlung_s$Cluster.no = str_extract(GDTlung_s$Cell_cluster, "(?<=_)[^_]+$")

Total_list <- GDTlung_s@meta.data %>% mutate(aa = cdr3_TRD) %>% filter(!is.na(cdr3_TRD))   %>% 
  split(f = .$Cluster.no)




Total_list1.trimmed <- GDTlung_trimmed@meta.data %>% mutate(aa = cdr3_TRD) %>%
  filter(!is.na(cdr3_TRD) & cdr3_TRD_freq >= 1)   %>%  split(f = .$Cell_cluster)

Total_list$P1$cdr3_TRD
MG<- c()

# Morisita_GDTlung<- clonalOverlap(Total_list, cloneCall="cdr3_TRD", method="morisita")+
  # theme_gray()+xlim(0,13)

# clonalOverlap(Total_list, cloneCall="cdr3_TRD", method="morisita")+
  theme_gray()+scale_x_discrete('L1')


# Morisita_GDTlung_trimmed <- clonalOverlap(Total_list1.trimmed, cloneCall="cdr3_TRD", method="morisita")+
  # theme_gray()

Morisita_GDTlung_trimmed$data


# library(conflicted)

# conflicted::conflict_prefer_all('dplyr')
GDTlung_s$cdr3_nt_TRD
Total_list1.trimmed_tcR <- GDTlung_s@meta.data %>% 
  mutate(CDR3.aa = cdr3_TRD,
  CDR3.nt = cdr3_nt_TRD) %>% filter(!is.na(CDR3.nt) & Cluster.no != "M" ) %>% 
  group_by(Cluster.no) %>% 
count(CDR3.nt, name = 'Clones') %>% 
  mutate(Proportion = Clones/sum(Clones)) %>%  as.data.frame %>% 

split(f = .$Cluster.no)


repOverlap(immdata$data, .method = "morisita")


Mori_result_TRD <- repOverlap(Total_list1.trimmed_tcR, .method = "morisita", .verbose = F, .col = "nt")

vis(Mori_result_TRD)+ 
  # ggtitle('Morisita index')+
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  ggtitle('Morisita index of TCRD CDR3')+xlab('Cluster')+ylab('Cluster')+xlab('Cluster')+ylab('Cluster')


repOverlap(Total_list1.trimmed_tcR, .method = "morisita", .verbose = F)  %>%  vis()
repOverlapAnalysis(Mori_result, "mds") %>% vis()

glimpse(immdata$data)

glimpse(Total_list1.trimmed_tcR)
immdata$data$MS1$Clones %>%  range





Patient_list<- SplitObject(GDTlung_s, split.by = 'patient') %>%  map( ~ .x@meta.data %>% filter(!is.na(cdr3_TRD))    %>%  split(f = .$Cell_cluster) )

Morisita_patient <- map(Patient_list,~ 
                          clonalOverlap(.x, cloneCall="cdr3_TRD", method="morisita")  )


clonalOverlap(Patient_list$p27, cloneCall="cdr3_TRD", method="morisita") 


test <- subset(GDTlung_s, subset = is.na(cdr3_TRD) )



 
 test<- clonalOverlap(combined3, cloneCall="cdr3_nt_TRD", method="morisita")
 
 
 clonalOverlap(combined2, cloneCall="cdr3_TRD", method="overlap")

# TCR length --------------------------------------------------------------


GDTlung_s$length_TRD

library(stringr)

GDTlung_s@meta.data %<>% mutate(length_TRD = case_when(!is.na(v_gene_TRD) ~ str_length(cdr3_TRD) ), 
                                length_TRG = case_when(!is.na(v_gene_TRG) ~ str_length(cdr3_TRG) ) )


ViolinPlot(GDTlung_s %>% subset(v_gene_TRD %in% c('TRDV1')), 'length_TRD',
           facet = "patient",
           box = T)

GDTlung_s@meta.data %>%  filter(v_gene_TRD == 'TRDV1'  
                               ) %>%  
   group_by(Cell_cluster, patient,v_gene_TRD) %>% 
   
  summarise(TRDlength = median(length_TRD),n =  n() ) %>%  
  filter(n > 10) %>% 
  ggplot(aes(x = Cell_cluster, y = TRDlength, color = patient, group= Cell_cluster))+
  geom_boxplot()+

  geom_point()+
  geom_line(aes(group = patient))


GDTlung_s$pheno



# geom_boxplot()


GDTlung_s$length_TRD
GDTlung_trimmed$length_TRD 

pub <- GDTlung_s@meta.data %>%  filter(public_TRD %in% c( 'publicTRDV1')) %>%  group_by(patient) 

table(pub$patient, pub$cdr3_TRD)


GDTlung_s$detected.no %>% table()


grep(shared, alluniqueTRD)

# Gini --------------------------------------------------------------------

library(reldist)

GDTRDfreq <-  GDTlung_s@meta.data%>% 
  filter(!is.na(cdr3_TRD)) %>% 
  dplyr::count( patient, Cell_cluster, cdr3_TRD) %>%
  arrange(desc(n))

GDTRDfreq



gini_GDTRD <-  GDTRDfreq %>% 
  filter(patient != 'p31') %>%
  dplyr::group_by(patient, Cell_cluster, ) %>%
  summarise(Gini_Index = gini(n), sum = sum(n))  %>%
  mutate(Gini_Index = replace(Gini_Index, sum  <= 10, NA))%>%
  # dplyr::rename(CD4orCD8 = CD4CD8 ) %>%
  ungroup() %>% 
  tidyr::complete(patient, Cell_cluster, 
                  fill = list(Gini_Index = NA, sum = NA))


giniindex_GDT <- ggplot(gini_GDTRD, aes(x = Cell_cluster , y = Gini_Index, 
                                           color = patient,group= Cell_cluster
                                        ))+
  stat_summary(fun = mean, geom = "bar", fill = 'transparent', color = 'black') +
  geom_point()+
  # geom_line(aes(group = Cell_cluster))+
  color_m(color = set_sample(umap.colors, s = 22))+theme_minimal()  +mytheme


giniindex_GDT



# SCENIC result -----------------------------------------------------------

library(SCopeLoomR)
loom <- open_loom('GDTlung_pyscenic/lung_output_tracks_mingene10.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
# 
# dim(regulonAUC)
# AUCmat <- AUCell::getAUC(regulonAUC) 
# 
# glimpse(AUCmat)
# rownames(AUCmat)
# rownames(AUCmat)   %<>% gsub("\\(\\+\\)", "_REG", .)
# 
# names(regulons) %<>%  gsub("\\(\\+\\)", "_REG", .)
# 
# print(AUCmat)
# 
# colnames(AUCmat)
# 
# GDTlung_s[['AUC']] <- CreateAssayObject(data = AUCmat)
# 
# GDTlung_s <- ScaleData(GDTlung_s, assay = 'AUC')
# 
# DE_Reglon <- FindAllMarkers(GDTlung_s, only.pos = T, assay = "AUC",min.pct = 0.25, test.use =  "wilcox", slot = "scale.data")

nrow(D)
DE_Reglon %>%  filter(cluster == "gd_TRM_Lu6")

ClusterCompare(GDTlung_s, id1 = "gd_TRM_LG6", id2 = "gd_TEMRA_LG3", 
               assay = "AUC",  log2fc = 0.25)

ClusterCompare(GDTlung_s, id1 = "gd_TRM_LG7", id2 = "gd_TEMRA_LG3", 
               assay = "AUC",  log2fc = 0.25)

AUCcells %>% dim()




gdAUCcells <- vroom::vroom("GDTlung_pyscenic/lung_output_tracks_mingene10.csv")  %>% `colnames<-`(
  gsub("\\(\\+\\)", "_REG", colnames(.) )  
) %>%  column_to_rownames("Cell") %>% t()

dim(GDTlung_s[['AUC']] )

dim(gdAUCcells)
rownames(gdAUCcells)

gdAUCcells

Reggd <-  cbind(t(gdAUCcells), FetchData(GDTlung_s, c("UMAP_1", "UMAP_2", "Cell_cluster")))


Feature_rast(Reggd, "TBX21_REG", colorgrd = "grd2")

  
GDTlung_s[['AUC']] <- NULL
GDTlung_s[['AUC']] <- CreateAssayObject(data = gdAUCcells)

GDTlung_s <- ScaleData(GDTlung_s, assay = 'AUC')
ggplot(Reggd, aes(x = Cell_cluster, y = TBX21_REG, fill = Cell_cluster))+geom_violin()+geom_boxplot()


Feature_rast(GDTlung_s, assay = "AUC", 
             g = c("GATA3-REG", "EOMES-REG", "TBX21-REG", "BATF-REG", "RORA-REG"), sz = 0.3, slot = "scale.data", colorgrd = "grd2")

Feature_rast(GDTlung_s, c('BATF-REG', "EOMES-REG", "GATA3-REG"), assay = "AUC", colorgrd = "grd2", slot = "scale.data", sz = 0.3) %>% list(
Feature_rast(Reggd, c("BATF_REG", "EOMES_REG", "GATA3_REG"), colorgrd = "grd2", sz = 0.3) ) %>% 
  PG(ncol = 1)


Feature_rast(Reggd, "EOMES_REG", colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  scale_m("grd2", c(0, 0.10,0.20, 0.30,0.40), c(0,0.4)),
  coord_fixed())   )

Feature_rast(Reggd, "GATA3_REG", colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  scale_m("grd2", c(0, 0.10,0.2), c(0,0.18)),
  coord_fixed())   )



GDTlung_s <- ScaleData(GDTlung_s, assay = 'AUC')

DE_Reglon <- FindAllMarkers(GDTlung_s, only.pos = T, assay = "AUC",min.pct = 0.25, test.use =  "wilcox", slot = "scale.data")

DE_Reglon %>% filter(gene == 'RORC-REG')


top5_reg <-  DE_Reglon %>% group_by(cluster) %>% top_n(5, avg_diff) %>%  pull(gene)


DoHeatmap(GDTlung_s, features = top5_reg, assay = "AUC") %>% heat_theme()

ClusterCompare(GDTlung_s, id1 = "gd_TRM_LG6", id2 = "gd_TEMRA_LG3", 
               assay = "AUC",  log2fc = 0.25)


reglon_TRM_TEMRA <-  ClusterCompare(GDTlung_s, id1 = "gd_TRM_LG6", id2 = "gd_TEMRA_LG3", test = "wilcox",
               assay = "AUC",  log2fc = 0.25)


reglon_TRM_TEMRA$table %>%  filter(gene %in% c("GATA3-REG", "EOMES-REG", "BATF-REG") )

reglon_TRM_TEMRA$table$gene
Feature_rast(GDTlung_s, assay = "AUC", slot = "scale.data", g = "RORC-REG")

# AUCmat <- AUCell::getAUC(vroom::vroom("GDTlung_pyscenic/lung_output_tracks_no_parameter.csv"))

rownames(AUCmat)


colnames(gdAUCcells) 

dim(gdAUCcells)

gdAUCcells$bc_backup

gdAUCcells$bc_backup <-  gdAUCcells$Cell
gdAUCcells$Cell <- NULL

GDTlung_s@meta.data  %<>%  select_at(.vars = vars(-contains("_REG")))


GDTlung_s@meta.data$RORC_REG <- NULL

gdAUCcells

GDTlung_s@meta.data  %<>% left_join(gdAUCcells, by = "bc_backup") %>% 
  `rownames<-`(GDTlung_s$bc_backup)

Feature_rast(GDTlung_s, c( "RORC_REG"),sz = 0.1
             )

Feature_rast(GDTlung_s, c( "RORC-REG"),sz = 0.1, assay = "AUC", slot = "scale.data"
)


ViolinPlot(GDTlung_s, "RORC_REG")

UMAPmeta <-  FetchData(CD4CD8, c("CD4CD8", "UMAP_1", "UMAP_2", "Cell_pheno", "bc_backup"))


data_frame(UMAPmeta, AUCcells)


UMAPmeta  %<>% data_frame(AUCcells) %>% as.data.frame()
head(UMAPmeta)

class(UMAPmeta)

dim(UMAPmeta)
class(UMAPmeta %>% as.data.frame())

Feature_rast(UMAPmeta, "RUNX3_REG")

gdregs <- GDTlung_s@meta.data %>%   select_at(.vars = vars(contains(c("_REG", "Cell_cluster"))))

write_csv(gdregs, "gdreg.csv")

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

ClusterCompare(CD4CD8, "TRM_1_P", "TRM_3_P", assay = "REG", log2fc = 0.01)




# GSEA --------------------------------------------------------------------


library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
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

TRMvsTEMRA

TRMvsTEMRA <- Genelist_generator(GDTlung_s, 'gd_TRM_LG6', 'gd_TEMRA_LG3')

GSEA_TRMvsTEMRA_allref<-GSEA(geneList = TRMvsTEMRA, TERM2GENE=ALL_msigdb_G,  
                       pvalueCutoff = 0.05, pAdjustMethod = "BH") 



view(GSEA_TRMvsTEMRA_allref@result)

# GSEAp2p8GO<-GSEA(geneList = p2p8en, TERM2GENE=GO,  nPerm = 100000, 
#                  minGSSize    = 15,
#                  pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")








write.xlsx(GSEAp2p8_c7@GSEAp4p6_c7, 'GSEAp4p6_c7.xlsx')

c('GSE25087_TREG_VS_TCONV_ADULT_UP',
  'GSE25087_TREG_VS_TCONV_FETUS_UP',
  'GSE7852_TREG_VS_TCONV_FAT_UP',
  'GSE14415_INDUCED_TREG_VS_FAILED_INDUCED_TREG_UP') %>%
  map(~ GSEA_multipplot(GSEAp4p6_c7,description_to_show = .x, base_size = 12,
                        title = str_remove(.x, 'GSE\\d+_'), legend.position = 'no',
                        c1='P2', c2 = 'P8' )) %>% PG()
view(GSEAp4p6_c7@result)


GSEA_multipplot(GSEAp4p6_c7,description_to_show = c(
  'GSE25087_TREG_VS_TCONV_ADULT_UP',
  'GSE7852_LN_VS_FAT_TREG_DN'
  
), 
                
                base_size = 8,
                # title = str_remove(.x, 'GSE\\d+_'),
                legend.position = 'bottom',
                c1='TRM_P6', c2 = 'TEMRA_P4' )



GSEA_multipplot(GSEAp2p8_c7)

gseaplot2(GSEAp4p6_c7,'GSE25087_TREG_VS_TCONV_ADULT_UP',title = 'TREG_VS_TCONV_FAT_UP')


GSEA_multipplot(GSEAp4p6_c7,description_to_show = c('GSE25087_TREG_VS_TCONV_ADULT_UP'
                                                    ), base_size = 10, 
                title = 'TREG_VS_TCONV_FAT_UP',
            legend.position = 'no',
                c1='P2', c2 = 'P8' )


gseaplot2(GSEAp4p6_c7,'GSE25087_TREG_VS_TCONV_ADULT_UP',base_size = 15 )

c('KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY') %>%
  map(~ GSEA_multipplot(GSEAp4p6KEGG,description_to_show = .x, base_size = 12,
                        title = str_remove(.x, 'KEGG_'), legend.position = 'no',
                        c1='P6', c2 = 'P4' )) %>% PG(ncol = 1)
gseaplot2(GSEAp4p6KEGG,'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',title = 'NATURAL_KILLER_CELL_\nMEDIATED_CYTOTOXICITY')

GSEA_multipplot(GSEAp4p6KEGG,description_to_show = c('KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'), base_size = 8,legend.position = 'bottom',
                title = "KEGG GSEA TRM vs TEMRA",
                
                c1='P6', c2 = 'P4' )

GSEA_multipplot(GSEA_TRMvsTEMRA_allref,description_to_show = c('GSE25087_TREG_VS_TCONV_ADULT_UP', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'), base_size = 8,legend.position = 'bottom',
                title = "KEGG GSEA TRM vs TEMRA",
                
                c1='TRM_LG6', c2 = 'TEMRA_LG3' )



GSEAp2p8KEGG@result$ID


view(GSEA_TRMvsTEMRA_allref@result)



# proportion of gdTs and CD4CD8s in each patient --------------------------

library("Hmisc")
library(corrplot)

gdTp <- table(GDTlung_s$patient, GDTlung_s$Cell_cluster) %>% prop.table(margin = 1)*100


gdCl <- colnames(gdTp)

abTp <- table(CD4CD8$patient, CD4CD8$Cell_pheno)%>% prop.table(margin = 1)*100

abCl <- colnames(abTp)


Tp <- cbind(gdTp, abTp) %>% as.data.frame()

Tp$patient <- rownames(Tp)



Feature_rast(CD4CD8)

Feature_rast(Tp, g = "patient", d1 = "gd_TRM_LG6", d2 = "TRM_3_LG", noaxis = F, axis.number = T, sz = 3)

Feature_rast(Tp, g = "patient", d1 = "gd_TEMRA_LG3", d2 = "Temra_1_LG", noaxis = F, axis.number = T, sz = 3)

r36 <- cor.test(Tp$gd_TRM_LG6, Tp$TRM_3_LG, method = "spearman")

r36$p.value
r36$estimate





result2 <- rcorr(as.matrix(Tp %>% select(-"patient")),type = 'spearman')

result2

result2$r
result2$n

corrplot(result2$r, type="lower", col = rev(COL2('RdBu', 200)),
         p.mat = result2$P, sig.level = 0.05, insig = "blank")



corrplot(result2$P, type="lower", 
         p.mat = result2$P, sig.level = 0.1, insig = "blank", col.lim = c(0,1))

?corrplot


library(pwr)


pwr.r.test(
  n = 7,
  r = 0.5,
  sig.level = NULL,
  power = 0.8,
  alternative = 'two.sided'
)

install.packages("presize")



# if we only look at cells from lung?

gdl <- GDTlung_s@meta.data %>%  filter(tissue == "Lung"  & !grepl("LN", Cell_cluster) ) %>%   mutate(Cell_cluster = as.character(Cell_cluster))


abl <-  CD4CD8@meta.data %>%  filter(tissue == "Lung"& !grepl("LN|uni", Cell_pheno)) %>% 
  mutate(Cell_pheno = as.character(Cell_pheno))

gdTp <-  table(gdl$patient, gdl$Cell_cluster) %>% prop.table(margin = 1)*100
abTp <-  
  table(abl$patient, abl$Cell_pheno) %>% prop.table(margin = 1)*100  
Tp <- cbind(gdTp, abTp) %>% as.data.frame()
as.matrix(Tp)

result_lung <- rcorr(as.matrix(Tp),type = 'pearson')

result_lung

corL <-  corrplot(result_lung$r, type="lower", col = rev(COL2('RdBu', 200)),tl.col="black",
         p.mat = result_lung$P, 
         sig.level = 0.05, insig = "blank") 

corL

Tp
Feature_rast(Tp, g = "gd_TRM_LG7", d1 = "gd_TRM_LG6", d2 = "TRM_3_LG", noaxis = F, axis.number = T, sz = 3)
Tpd <- Tp %>% as.data.frame() %>%  rownames_to_column("patient")

Feature_rast(Tpd , g = "patient",do.label = F,
             d1 = "gd_TEMRA_LG3", d2 = "Temra_1_LG", noaxis = F, 
             axis.number = T, sz = 3)+
  geom_smooth(method = "lm", se = F)+
  stat_cor(method = "spearman")

Feature_rast(Tpd , g = "patient",do.label = F,
             d1 = "gd_TEMRA_LG3", d2 = "Temra_1_LG", noaxis = F, 
             axis.number = T, sz = 3)+
  geom_smooth(method = "lm", se = F)+
  stat_cor(method = "pearson")

Feature_rast(Tpd , g = "patient",do.label = F,
             d1 = "gd_TRM_LG5", d2 = "TCM_3_M", noaxis = F, 
             axis.number = T, sz = 3)+
  geom_smooth(method = "lm", se = F)+
  stat_cor(method = "spearman")

Feature_rast(Tpd , g = "patient",do.label = F,
             d1 = "gd_TEMRA_LG4", d2 = "TCM_2_M", noaxis = F, 
             axis.number = T, sz = 3)+
  geom_smooth(method = "lm", se = F)+
  stat_cor(method = "spearman")


Feature_rast(Tpd , g = "patient",do.label = F,
             d1 = "gd_Naive_LN1", d2 = "Naive_1_LN", noaxis = F, 
             axis.number = T, sz = 3)+
  geom_smooth(method = "lm", se = F)+
  stat_cor(method = "spearman")




gdln <- GDTlung_s@meta.data %>%  filter(tissue == "LLN") 
abln <-  CD4CD8@meta.data %>%  filter(tissue == "LLN")

gdTpln <-  table(gdln$patient, gdln$Cell_cluster) %>% prop.table(margin = 1)*100
abTpln <-  
  table(abln$patient, abln$Cell_pheno) %>% prop.table(margin = 1)*100  
Tpln <- cbind(gdTpln, abTpln) %>% as.data.frame()


result_ln <- rcorr(as.matrix(Tpln),type = 'spearman')

corrplot(result_ln$r, type="lower", col = rev(COL2('RdBu', 200)),
         p.mat = result_ln$P, sig.level = 0.05, insig = "blank")



gdcor <-  rcorr(cbind(gdTp, gdTpln) %>% as.matrix(),type = 'spearman') 
  corrplot(gdcor[['r']], type="lower", col = rev(COL2('RdBu', 200)),
           p.mat = gdcor[['P']], sig.level = 0.05, insig = "blank")


rcorr(cbind(gdTp, gdTpln) %>% as.matrix(),type = 'spearman') %>% print(.[["n"]]+ .[["P"]])



abcor <-  rcorr(cbind(abTp, abTpln) %>% as.matrix(),type = 'spearman') 
corrplot(abcor[['r']], type="lower", col = rev(COL2('RdBu', 200)),
         p.mat = abcor[['P']], sig.level = 0.01, insig = "blank")
                                                                                                                                                                                                                                             print(result_ln[["n"]])
# !!!!!!!!!!!stop!!!!!!!!!!!!!!!!!!!! -------------------------------------




# integration with gdT from adult PBMC ------------------------------------

Adult_GDT_2020AUG <- readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds') %>% 
  subset(group == 'Adult') 
Feature_rast(Adult_GDT_2020AUG)







Adult_GDT_2020AUG%<>%   ScaleData(assay = 'RNA', 
                                  vars.to.regress = c("percent.mito", 'orig.ident',
                                                      "S.Score", 'G2M.Score',
                                                      "percent.ribo",
                                                      'nCount_RNA','nFeature_RNA' ) )

DefaultAssay(Adult_GDT_2020AUG) <- 'RNA'

GDTlung_s$tissue
Adult_GDT_2020AUG$patient <- Adult_GDT_2020AUG$orig.ident
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


saveRDS(GDTlung_s, 'GDTlung_newTCR_2021.0914.rds')
GDTlung_inte  %<>% AddMetaData(GDTlung_inte@assays$CITE@data %>% t())
GDTlung_inte@meta.data

GDTlung_inte[['integrated']]<- NULL
GDTlung_inte[['HTO']]<- NULL
GDTlung_inte[['CITE']]<- NULL
GDTlung_inte[['GM']]<- NULL


# integtation by seurat integration
Adult_GDT_2020AUG
l.Adult_GDT_2020AUG <- SplitObject(Adult_GDT_2020AUG, split.by = "patient")
l.gdtlung <- SplitObject(GDTlung_inte, split.by = "patient")

l.all <- c(l.Adult_GDT_2020AUG, l.gdtlung)

l.all <-   lapply(X = l.all, FUN = function(x) {
  # x  <- AddMetaData(x, t(x[['CITE']]@data)) 
  
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

for (x in names(l.all))      {
  
  l.all[[x]]@assays$RNA@var.features %<>% str_subset('^HIST|^MT|^IG|^TRG|^TRD|^HSP|^RP', negate = T) 
  
}
anchors <- FindIntegrationAnchors(object.list = l.all, 
                                  dims = 1:50)
GDT_Lung_BL<- IntegrateData(anchorset = anchors, dims = 1:50)

DefaultAssay(GDT_Lung_BL) <- 'integrated'

GDT_Lung_BL %<>% ScaleData( vars.to.regress = c(
  "percent.mito",
  # 'project',
  "percent.ribo",
  'nCount_RNA',
  'nFeature_RNA' )) %>% 
  RunPCA(npcs = 100, verbose = T,nfeatures.print = 40)

ElbowPlot(GDT_Lung_BL,ndims = 100)

GDT_Lung_BL <- RunUMAP(object = GDT_Lung_BL, dims = 1:50, 
                       reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:70))
for (i in seq(0.6,1.6,0.1) %>% rev()) {
  GDT_Lung_BL <- FindClusters(GDT_Lung_BL, resolution = i)
}

Feature_rast(GDT_Lung_BL, c('ident','tissue','project', 'big_cluster', 'v_gene_TRD', 'patient'), ncol = 3) 
saveRDS(GDT_Lung_BL, 'GDT_Lung_BL_integrated_pb_and_lung.rds')

GDT_Lung_BL <- readRDS('GDT_Lung_BL_integrated_pb_and_lung.rds')
 

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


# Figs --------------------------------------------------------------------


Feature_density(GDTlung_s, 'IL1B1')




Feature_rast(GDTlung_s, 'paired_sp', colorset = c(umap.colors[1:9], 'lightblue'), sz = 1.2)


umap_idnet_lunggdt <- Feature_rast(GDTlung_s, c('ident'), noaxis = F, labelsize = 12)+ggtitle('Phenotypes of lung γδT cells')+
  theme(plot.title = element_text(size = 12, face = 'plain'))


barplot(101:123, col = umap.colors, names.arg = 1:23)


vg9d2 <- Feature_rast(GDTlung_s, 'Vg9Vd2',  do.label = F, noaxis = F, mythe = F)+ggtitle('γδTCR')+
  theme(plot.title = element_text(size = 12, face = 'plain'))
tissueUMAP <- Feature_rast(GDTlung_s, 'tissue', do.label = F, noaxis = F, mythe = F, sz = 0.5)+scale_color_manual(labels = c('LN', 'Lung'), values = umap.colors)+
  theme(plot.title = element_text(size = 12, face = 'plain'))




GM_D <- Feature_rast(GDTlung_s, 'GM_D', do.label = F, noaxis = F, color_grd = 'grd', mythe = F)+ggtitle('Gene module: type-3 immunity')+
  theme(plot.title = element_text(size = 12, face = 'plain'))


PG(list(vg9d2,tissueUMAP, GM_D), ncol = 3, rw = c(1, 0.85, 0.85)) %T>% figsave('Vg9Vd2areType3.v2.pdf', 220, 80) 



Feature_rast(GDTlung_s, 'GM_G', do.label = F, noaxis = T, color_grd = 'grd')


GMD_vio <- ViolinPlot(GDTlung_s, 'GM_D', colors = umap.colors, box = T)+ggtitle('Gene module: type-3 immunity')+
  theme(plot.title = element_text(size = 12, face = 'plain'))


PG(list(umap_idnet_lunggdt, vg9d2, GMD_vio), ncol = 3) %T>% figsave('Vg9Vd2areType3', 200, 80) 

GDTlung_s$public_TRD


Feature_rast(GDTlung_s, 'public_TRD', colorset = umap.colors[c(21, 5,11,6)], do.label = F)

GDT_2020AUG <- readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds')


Feature_rast(GDT_2020AUG, 'ZBTB16')+ggtitle('ZBTB16 (PLZF)')+
  theme(plot.title = element_text(size = 16))


# Trm vs Tcirc 

TrmTc <- c('ITGAE', 'ITGA1', 'ZNF683', 'ITGB2', 'KLF2', 'KLRG1')

Feature_density(GDTlung_s, TrmTc, ncol = 3, titlesize = 15)


TrmTcprotien <- c('CD103', 'CD49a', 'KLRG1') %>% paste0('.protein')

Feature_density(GDTlung_cite, TrmTcprotien, ncol = 3, titlesize = 15)

ctl <- c('CD8A', 'CD8B', 'EOMES', 'FCGR3A', 'NKG7', 'GZMA')

Feature_density(GDTlung_s, ctl, ncol = 6, titlesize = 15)

c('FOXP3',  'IL2RA',  'IKZF2', 'CTLA4', 'TNFRSF18',  'LGALS3')

Feature_density(GDTlung_s, c('FOXP3',  'IL2RA',  'GATA3', 'CTLA4', 'TNFRSF18',  'LGALS3',ctl), ncol = 6, titlesize = 15)
  
Feature_rast(GDTlung_s, c('FOXP3',  'IL2RA',  'IKZF2', 'CTLA4', 'TNFRSF18',  'LGALS3',ctl), ncol = 6, titlesize = 15)

Cytokines <- c( 'IL10', 'AREG','TGFB1', 'IL2','IFNG','TNF'
)


Cytokines %>% Feature_density(GDTlung_s, ., ncol =6, titlesize = 15)

Cytokines %>% Feature_rast(GDTlung_s, ., ncol =6, titlesize = 15)

Cytokines %>% ViolinPlot(GDTlung_s, ., ncol =3, colors = umap.colors)

Feature_density(GDTlung_s, c('AREG', 'CSF2', 'IFNG', 'PRF1'), ncol =4, titlesize = 15)

Feature_rast(GDTlung_s, 'Vg9Vd2', labelsize = 10, colorset = c('red', alpha('cyan',0.8)), navalue = alpha('transparent',0.5))


Feature_density(GDTlung_s, c('AREG', 'CSF2', 'IFNG', 'IL2'), ncol =6, titlesize = 15)



PPTf1 <- Feature_rast(GDTlung_s, c('ident'), noaxis = F, labelsize = 14, mythe = F)+NoLegend()+ggtitle('clustering')+
  theme(plot.title = element_text(size = 14))
PPTf2 <- Feature_rast(GDTlung_s, 'v_gene_TRD', do.label = F, mythe = F, noaxis = F)+
  theme(plot.title = element_text(size = 14))
PPTf3 <- Feature_rast(GDTlung_s, 'cdr3_TRD_perc', navalue = 'transparent', mythe = F, noaxis = F)+
  ggtitle('TRD CDR3 clonality')+
  theme(plot.title = element_text(size = 14))


PG(list(PPTf1, PPTf2,PPTf3), ncol = 3, rw = c(1, 1.15, 1.15))

# export to H5ad ----------------------------------------------------------
library(SeuratData)
library(SeuratDisk)
Feature_rast(GDTlung_s, 'Cell_cluster', facets = 'patient')

table(GDTlung_s$patient, GDTlung_s$pheno) %>%  prop.table(margin = 1)*100

GDTlung_trimmed <- subset(GDTlung_s, pheno %in% c('Lung_memory',   'Lung_TRM'  ,   'LN_memory') & 
                            patient %in% c('p25', 'p27','p32', 'p45', 'p71', 'p77'))


table(GDTlung_trimmed$patient, GDTlung_trimmed$pheno) %>%  prop.table(margin = 1)*100
Feature_rast(GDTlung_trimmed)

GDTlung_trimmed$Cell_cluster %>% unique()


GDTlung_trimmed$Cell_cluster <- as.vector(GDTlung_trimmed$Cell_cluster )
GDTlung_trimmed@meta.data  %<>%  mutate_if(is.factor, as.character) %>%  `rownames<-`(GDTlung_trimmed$bc_backup)

GDTlung_s_trim <- GDTlung_s %>% DietSeurat(dimreduces = "umap", assays = "RNA")
GDTlung_s_trim@assays$RNA$scale.data <-  NULL

GDTlung_s_trim@reductions$umap <-  GDTlung_s@reductions$umap

Feature_rast(GDTlung_s_trim)


GDTlung_s_trim@meta.data


GDTlung_s_trim@reductions$umap@key

GDTlung_s_trim@

SaveH5Seurat(GDTlung_s_trim, 'GDTlung_pyscenic/GDTlung.h5Suerat', overwrite = T)
# SaveH5Seurat(GDTlung_trimmed, 'GDTlung.trimmed')

Convert('GDTlung_pyscenic/GDTlung.h5Suerat.h5seurat', dest = 'h5ad', overwrite = T)

write.csv(GDTlung_s@meta.data, "GDTlung_pyscenic/GDTmeta.csv")
write.csv(FetchData(GDTlung_s, c("UMAP_1", "UMAP_2")), "GDTlung_pyscenic/GDTumap.csv")


# Convert('GDTlung.trimmed.h5seurat', dest = 'h5ad')
rm(GDTlung_s_trim)


# Public data fetal lung --------------------------------------------------

# library("anndata")
# renv::init()
# 
# renv::deactivate()
# 
# renv::activate()
# 
# Convert('Assembled10DomainsFiltered_WithRaw.h5ad', dest = "h5seurat", overwrite = TRUE)
# FetalLungT <- LoadH5Seurat(file = "Assembled10DomainsFiltered_WithRaw.h5seurat")
# 
# 
# Convert('He_et_al_fetallung.h5seurat', dest = 'h5ad')
# data <- read_h5ad("He_et_al_fetallung.h5ad")
# 
# lungH5ad <- read_h5ad('Assembled10DomainsFiltered_WithRaw.h5ad')
# str(lungH5ad)
# 
# glimpse(data) %>% View()
# data$var
# 
# data <- CreateSeuratObject(counts = t(lungH5ad$X), meta.data = lungH5ad$obs)
# 
# data$productive_summary
# data$obs$receptor_subtype %>%  table()
# 
# data$receptor_subtype %>%  unique()
# 
# 
# data$obs$productive_summary
# 
# rm(data)
# 
# data$varp
# 
# Feature_rast(data)
# 
# 
# GDTlung_s@assays$RNA@counts

# VIPER and PISCES --------------------------------------------------------
library(PISCES)

DefaultAssay(GDTlung_s) <- 'integrated'

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


Feature_rast(GDTlung_s, 'ident',  other = c('cdr3_TRD', 'patient'))+ geom_line(
  data = FetchData(subset(GDTlung_s, cdr3_TRD_freq >2), c('UMAP_1', 'UMAP_2', 'cdr3_TRD','patient')),
  aes(group = cdr3_TRD),color = alpha('black', 0.3))+
  facet_wrap(~patient)+
  NULL





Feature_rast(GDTlung_s, 'tissue', facets = 'patient')


GDTlung_s@meta.data$cdr3_paired

GDTlung_s@meta.data %>% top_n(1000,cdr3_TRD_freq) %>% pull(cdr3_TRG) %>% unique() %>% length()

GDTlung_s@meta.data %>% top_n(1000,cdr3_TRD_freq) %>% filter(!is.na(cdr3_paired)) %>% arrange(cdr3_TRD) %>% 
  select(cdr3_TRD, cdr3_TRD_freq, cdr3_TRG, cdr3_TRG_freq) %>% distinct()



DoHeatmap(GDTlung_s, c('CCR6', 'DPP4', 'KLRB1', 'RORC', 'IL23R', 'BLK', 'ZBTB16')) %>% heat_theme()



# barcodes 

allbc <- GDTlung_s@meta.data %>% group_split(orig.ident) %>% map(~ .x %>%pull(bc_backup) %>% str_remove(  '_\\d$') %>% as.vector ) %>% setNames(GDTlung_s$orig.ident %>% unique())

writeLines(allbc$p25, 'raw/surface_Falk1_gd/outs/filtered_feature_bc_matrix/filteredBC.tsv')

writeLines(allbc$`p32&p45`, 'raw/surf_322_325_453_456_gd/outs/filtered_feature_bc_matrix/filteredBC.tsv')

writeLines(allbc$`p71&p31&p27`, 'raw/surface_falk3_gd/outs/filtered_feature_bc_matrix/filteredBC.tsv')


writeLines(allbc$p73p77, 'raw/surface_falk77_73_gd/outs/filtered_feature_bc_matrix/filteredBC.tsv')
# public TRD --------------------------------------------------------------

alluniqueTRD <- readRDS('/home/big/tanlikai/Human_GDT_2019/alluniqueTRD.rds')

lung_uniqueTRD <- GDTlung_s@meta.data %>% filter(!is.na(v_gene_TRD)) %>% group_by(patient) %>% count(cdr3_TRD) %>%  pull(cdr3_TRD)

alluniqueTRD <- c(alluniqueTRD, lung_uniqueTRD)

publicTRD_ALL <- as_data_frame(alluniqueTRD) %>% count(value) %>%  filter(n >1) %>% arrange(-n) %>% 
  dplyr::rename('cdr3_TRD' = 'value' , 'detected.no' = 'n' )

GDTlung_s@meta.data %>% filter(!is.na(v_gene_TRD)) %>% mutate(cdr3_TRD = paste0(v_gene_TRD, cdr3_TRD))  %>% group_by(patient) %>% count(cdr3_TRD) %>%  pull(cdr3_TRD) %>% 
  
  as_data_frame() %>% count(value) %>%  filter(n >0) %>% arrange(-n) %>% 
  dplyr::rename('cdr3_TRD' = 'value' , 'detected.no' = 'n' )


grep('CALGELYGPLYWGPTPRTTDKLIF', alluniqueTRD)





publicTRD_ALL

GDTlung_s@meta.data  %<>%  left_join(publicTRD_ALL, by = 'cdr3_TRD') %>% 
  dplyr::mutate(public_TRD = case_when(
    cdr3_TRD %in% publicTRD_ALL$cdr3_TRD ~ paste0('public', v_gene_TRD),
    !is.na(cdr3_TRD) ~ 'private TRD'
  )) %>%  dplyr::mutate(
    public_TRD_freq = case_when(
      detected.no > 50 ~"most public (>50)",
      detected.no %in% 21:50 ~  "common public (21 ~ 50)",
      detected.no %in% 5:20 ~  "public (5 ~ 20)",
      detected.no %in% 2:4 ~ "rare public (2 ~ 4)",
      !is.na(cdr3_TRD) ~ 'private TRD'
    )  ) %>% mutate(public_TRD_freq = factor(public_TRD_freq, levels = c(
      "most public (>50)","common public (21 ~ 50)",
      "public (5 ~ 20)","rare public (2 ~ 4)",'private TRD'
    ))) 

GDTlung_s@meta.data%<>%  `rownames<-`(GDTlung_s$bc_backup)

Feature_rast(GDTlung_s, 'public_TRD', do.label = F, facets = 'patient')

GDTlung_s$detected.no %>%  table()


# FACS data  ----------------------------------------------------- --------

  
  GDTfacs <-  vroom::vroom("FACSdata/gdT_corrected_FACS.csv")
# TotalTfacs <- vroom::vroom("FACSdata/TotalT_corrected_FACS.csv")



# Vg9Vd2   
GDTfacs %<>% mutate(TCRtype = case_when(gdTtype == "NonVd2" ~ "NonVg9Vd2", TRUE ~ gdTtype))


Feature_rast(GDTfacs, g = "TCRtype", d1 = "Vd2", d2 = "Vg9", facets = "tissue", othertheme = list(coord_fixed()))+theme_bw()+mytheme+theme(panel.grid = element_line(size = 0) )+
  geom_hline(yintercept = 2.5, linetype  = 'dotted')+geom_vline(xintercept = 3.2, linetype = 'dotted')


vg9vd2st <- GDTfacs %>% group_by(batch, tissue, TCRtype) %>%  summarise(n = n()) %>% 
 ungroup() %>% group_by( tissue, batch) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
    Vg9Vd2 = sum(n[TCRtype == "Vg9Vd2"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(PercentVg9Vd2 = (Vg9Vd2 / Total) * 100)  %T>% print() %>% 
  group_by(tissue) %>% 
  summarise(mean = mean(PercentVg9Vd2) %>% round(2), sd =sd(PercentVg9Vd2) %>% round(2)  )%T>% print()

Feature_rast(TotalTfacs, g = "tissue", d1 = "TCRgd", d2 = "Vd2", noaxis = F, sz = 0.2)
# tissue batch Total Vg9Vd2 PercentVg9Vd2
# <chr>  <chr> <int>  <int>         <dbl>
#   1 LN     p25     983     80          8.14
# 2 LN     p31     468     35          7.48
# 3 LN     p32     199     13          6.53
# 4 LN     p77    1518    131          8.63
# 5 Lung   p25    3000     63          2.1 
# 6 Lung   p31     403    112         27.8 
# 7 Lung   p32    1260    101          8.02
# 8 Lung   p77    3000    108          3.6 

# tissue  mean     sd
# <chr>  <dbl>  <dbl>
#   1 LN      7.69  0.907
# 2 Lung   10.4  11.9 

Feature_rast(GDTfacs, g = "TCRtype", d1 = "Vd2", d2 = "Vg9", facets = "tissue", othertheme = list(coord_fixed()))+
  geom_hline(yintercept = 2.5, linetype  = 'dotted')+geom_vline(xintercept = 3.2, linetype = 'dotted')+
  geom_text(data = vg9vd2st, aes(label =paste0(mean,"%±",sd), x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1)+
  theme_bw()+mytheme+theme(panel.grid = element_line(size = 0) )



# CD26
GDTfacs  %<>% mutate(CD26pos = case_when(CD26 > 3 ~ "CD26pos", TRUE ~ "CD26neg"))
GDTfacs$CD26









stCD26 <- GDTfacs %>% group_by(batch, CD26pos, TCRtype) %>%  summarise(n = n()) %>% 
  ungroup() %>% group_by( TCRtype, batch) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
  CD26pos = sum(n[CD26pos == "CD26pos"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(PercentCD26 = (CD26pos / Total) * 100)  %T>% print() %>% 
  group_by(TCRtype) %>% 
  summarise(mean = mean(PercentCD26) %>% round(2), sd =sd(PercentCD26) %>% round(2)  )%T>% print()

stCD26

ggplot(GDTfacs , aes_string(x = "CD26", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + 
  geom_text(data = stCD26[1,], aes(label = paste0("CD16+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 0.7) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  geom_text(data = stCD26[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 0.5) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme



# KLRG1 
ggplot(GDTfacs , aes_string(x = "KLRG1", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0)+
  geom_vline(xintercept = 2.6, linetype = 'dotted')+fill_m()+mytheme


GDTfacs  %<>% mutate(KLRG1pos = case_when(KLRG1 > 2.6 ~ "KLRG1pos", TRUE ~ "KLRG1neg"))

stKLRG1 <- GDTfacs %>% group_by(batch, KLRG1pos, TCRtype) %>%  summarise(n = n()) %>% 
  ungroup() %>% group_by( TCRtype, batch) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
    KLRG1pos = sum(n[KLRG1pos == "KLRG1pos"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(PercentKLRG1= (KLRG1pos / Total) * 100)  %T>% print() %>% 
  group_by(TCRtype) %>% 
  summarise(mean = mean(PercentKLRG1) %>% round(2), sd =sd(PercentKLRG1) %>% round(2)  )%T>% print()
stKLRG1


ggplot(GDTfacs , aes_string(x = "KLRG1", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + 
  geom_text(data = stKLRG1[1,], aes(label = paste0("KLRG1+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 1.4) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  geom_text(data = stKLRG1[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 1.2) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme



# CD103

ggplot(GDTfacs , aes_string(x = "CD103", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0)


GDTfacs$TRM
stTRM <- GDTfacs %>% group_by(batch, TRM, TCRtype) %>%  summarise(n = n()) %>% 
  ungroup() %>% group_by( TCRtype, batch) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
    TRM = sum(n[TRM == "TRM"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(PercentTRM= (TRM / Total) * 100)  %T>% print() %>% 
  group_by(TCRtype) %>% 
  summarise(mean = mean(PercentTRM) %>% round(2), sd =sd(PercentTRM) %>% round(2)  )%T>% print()

ggplot(GDTfacs , aes_string(x = "CD103", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + 
  geom_text(data = stTRM[1,], aes(label = paste0("CD103+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 2.2) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  geom_text(data = stTRM[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 1.8) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme





 <- GDTfacs %>% group_by(batch, TRM, KLRG1pos, CD26pos, TCRtype) %>%  summarise(n = n()) %>% 
  ungroup() %>% group_by( TCRtype, batch) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
    CD103pos = sum(n[TRM == "TRM"]), 
    KLRG1pos = sum(n[KLRG1pos == "KLRG1pos"]), 
    CD26pos = sum(n[CD26pos == "CD26pos"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(
    PercentCD103= (CD103pos / Total) * 100,
    PercentKLRG1 = (KLRG1pos / Total) * 100,
    PercentCD26 = (CD26pos / Total) * 100)

c("PercentCD103", "PercentKLRG1", "PercentCD26")  %>% 
  map(~ ggplot(TRMKLRG1CD26_inVg9Vd2 , aes_string(x = "TCRtype", y = .x ))+geom_point(aes(color = TCRtype))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
        stat_compare_means(  paired = T, method = 't.test')
        
        ) %>% PG(ncol = 3)
+TRMKLRG1CD26_inVg9Vd2[c(1,5,2,6,3,7,4,8),]
compare_means(PercentCD103      ~ TCRtype   ,  data = TRMKLRG1CD26_inVg9Vd2[c(1,6,2,5,3,7,4,8),],
              paired = T,
              method = "t.test")



ggplot(TRMKLRG1CD26_inVg9Vd2, aes_string(x = "TCRtype", y = "PercentCD26"  ))+geom_point(aes(color = TCRtype))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
  stat_compare_means(  paired = T, method = "t.test")

GDTfacs$gdTtype

Feature_rast(GDTfacs %>%  filter(gdTtype == "NonVd2"), "tissue", d1 = "CD103", d2 = "CD49a", noaxis = F, axis.number = T , do.label = F, sz = 0.3, othertheme = list(coord_fixed()))  +theme_bw(base_line_size = 0)+
  ggtitle("Tissue residency")+mytheme


ggplot(GDTfacs %>%  filter(gdTtype == "NonVd2"), aes_string(x = "CD103", fill = "tissue"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0)+
  geom_vline(xintercept = 3, linetype = 'dotted')+fill_m()+mytheme
GDTfacs$Eomes_MFI




ggplot(TotalTfacs, aes_string(x = "GMCSF", fill = "TRM"))+geom_density(alpha = 0.5)+
  facet_wrap(~gdTcells+tissue)+geom_vline(xintercept = 3)


Feature_rast(TotalTfacs, "TRM", d1 = "TCRgd", d2 = "CD357", noaxis = F, axis.number = T , do.label = F, sz = 0.1)  + facet_wrap(~gdTcells+tissue)+
  geom_hline(yintercept = 3)



Feature_rast(TotalTfacs, "TRM", d1 = "TCRgd", d2 = "Eomes", noaxis = F, axis.number = T , do.label = F, sz = 0.1)  + facet_wrap(~gdTcells+tissue)+
  geom_hline(yintercept = 3)

Feature_rast(TotalTfacs, "TRM", d1 = "TCRgd", d2 = "GzmA", noaxis = F, axis.number = T , do.label = F, sz = 0.1)  + facet_wrap(~gdTcells+tissue)+
  geom_hline(yintercept = 2.5)



Feature_rast(GDTfacs, "TRM", d1 = "TCRgd", d2 = "CD357", noaxis = F, axis.number = T , do.label = F, sz = 1)  + facet_wrap(~batch+tissue)+
  geom_hline(yintercept = 3)




Feature_rast(GDTfacs, "TRM", d1 = "CD25", d2 = "CD127", noaxis = F, axis.number = T , do.label = F, sz = 1)  + facet_wrap(~batch+tissue)+
  geom_hline(yintercept = 3)






ggplot(GDTfacs, aes_string(x = "CD357", fill = "TRM"))+geom_density(alpha = 0.5)+
  # facet_wrap(~tissue)+
  geom_vline(xintercept = 3.2, linetype="dotted")+theme_classic()+mytheme


GDTfacs %>%  filter(TCRtype != "Vg9Vd2") %>% group_by(batch, TRM) %>%  summarise(medvalue= median(CD357_MFI), CD357_sd = sd(CD357_MFI)) %>%ungroup() %>%   compare_means(medvalue ~ TRM, paired = T, data = ., method = "wilcox.test")



GDTfacs %>% group_by(batch, TRM) %>% count()




results <- by(GDTfacs, GDTfacs$batch, function(subdata) {
  TRM <- subdata[subdata$TRM == "TRM", "CD357_MFI"]
  
  print(TRM)  # This will print the data for each patient separately
  # You can replace this with any statistical test or calculation
})

perm_test <- function(x, y, n.perm = 1000) {
  combined <- c(x, y)
  n.x <- length(x)
  test.stat <- mean(x) - mean(y)
  perm.stats <- replicate(n.perm, {
    permuted <- sample(combined)
    mean(permuted[1:n.x]) - mean(permuted[(n.x + 1):length(permuted)])
  })
  p.value <- mean(abs(perm.stats) >= abs(test.stat))
  list(test.statistic = test.stat, p.value = p.value)
}

GDTfacs %>%  filter(TCRtype != "Vg9Vd2") %>% group_by(batch, TRM) %>%  summarise(medvalue= median(GzmA), CD357_sd = sd(CD357_MFI)) %>%ungroup() %>%   
  compare_means(medvalue ~ TRM, paired = T,
  data = ., method = "t.test")

GDTfacs %>%  filter(TCRtype != "Vg9Vd2") %>% group_by(batch, TRM) %>%  summarise(medvalue= median(Eomes_MFI), CD357_sd = sd(CD357_MFI)) %>%ungroup() %>%   compare_means(medvalue ~ TRM, paired = T, data = ., method = "wilcox.test")



results <- by(GDTfacs, GDTfacs$batch, function(subdata) {
  TRM <- subdata[subdata$TRM == "TRM", "CD357_MFI"]
  Tcirc <- subdata[subdata$TRM == "Tcirc", "CD357_MFI"]
  perm_test(TRM, Tcirc) %T>% print()
})

print(results)


library(lme4)

# Assuming 'Expression' is the response, 'TissueType' as fixed effect, and 'PatientID' as a random effect
model <- lmer(CD357_MFI ~ TRM + (1|batch), data = GDTfacs, REML = FALSE)
summary(model)

par(mfrow=c(2,2))
plot(model)


lm2 <- lm(CD357_MFI ~ TRM + batch, data = GDTfacs)
summary(lm2)


library(glmmTMB)

glmmTMB(CD357_MFI ~ TRM + (1|batch), data = GDTfacs, REML = FALSE)



Feature_rast(GDTfacs, g = "tissue", d1 = "CD45RA", d2 = "CD127", noaxis = F)

GDTfacs  %<>% mutate(Gata3Gate = case_when(GATA3 > 2.5 ~ "GATA3+", TRUE ~ "GATA3-"),
                     CD357Gate = case_when(CD357 > 3 ~ "CD357+", TRUE ~ "CD357-"), 
                     GzmaGate =  case_when(GzmA > 2.5 ~ "GzmA+", TRUE ~ "GzmA-") )

GDTfacs_lung_Vd1 <- filter(GDTfacs, tissue == "Lung",Vg9Vd2 != "Vg9Vd2" )

map2(c("GATA3",  "CD357", "GzmA", "Eomes"),list(2.5,3,2.5,NULL) ,~
      ggplot(GDTfacs_lung_Vd1, aes_string(x = .x, fill = "TRM"))+geom_density(alpha=0.4)+
     theme_classic()+mytheme+fill_m()+NoLegend()+ geom_vline(xintercept = .y, linetype = 'dotted') ) %>% PG(ncol = 4, labels = 'Lung Vd1/3 T cells')


Feature_rast(GDTfacs_lung_Vd1,"Gata3Gate", d1 = "CD357", d2 = "GATA3", facets = c("batch", "TRM") )



Vd1Lungmarkers <- GDTfacs_lung_Vd1 %>% group_by(batch, TRM, Gata3Gate, CD357Gate,GzmaGate) %>%  summarise(n = n()) %>% 
  ungroup() %>% group_by( batch,TRM) %>%  summarise(
    Total = sum(n),                             # Total cells per batch and tissue
    Gata3pos = sum(n[Gata3Gate == "GATA3+"]), 
    CD357pos = sum(n[CD357Gate == "CD357+"]), 
   GzmApos = sum(n[GzmaGate == "GzmA+"]),       # Sum of Vg9Vd2 cells
    .groups = "drop"
  ) %>%
  mutate(
    Gata3Pos= (Gata3pos / Total) * 100,
    CD357Pos = (CD357pos / Total) * 100,
    GzmAPos = (GzmApos / Total) * 100)

Vd1Lungmarkers
Vd1Emmes <-   GDTfacs_lung_Vd1 %>% group_by(batch, TRM) %>% 
  summarise(EomesMFI_mean = mean(Eomes_MFI), EomesMFI_median = median(Eomes_MFI))
Vd1Emmes

Vd1Lungmarkers$EomesMFI_mean <-  Vd1Emmes$EomesMFI_mean
Vd1Lungmarkers
 c("Gata3Pos", "CD357Pos", "GzmAPos", "EomesMFI_mean")  %>% 
  map(~ ggplot(Vd1Lungmarkers , aes_string(x = "TRM", y = .x ))+geom_point(aes(color = TRM   ))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
        stat_compare_means(  paired = T, method = 't.test', )
      
  ) %>% PG(ncol = 4)


 
 
 
 
 
 

# TCR sharing with David Vermijlen data -----------------------------------
DVdata <- list.files('/home/big/tanlikai/Fetal_thymus_TCR_David_V/',full.names = T) %>% 
  map( ~ read.table(.x, header = T))



DVdata_D  <- list.files('/home/big/tanlikai/Fetal_thymus_TCR_David_V/') %>% 
  str_subset('_D.txt')%>% 
  map( ~ read.table(paste0('/home/big/tanlikai/Fetal_thymus_TCR_David_V/', .x), header = T) %>%  mutate(sample = .x)) %>%  purrr::reduce(.f = rbind) 


DV_TRD <- DVdata_D %>%  pull(cdr3aa)


lungTCRs <- GDTlung_s@meta.data %>% filter(!is.na(v_gene_TRD)) %>% pull(cdr3_TRD)


shared <- intersect(DV_TRD, lungTCRs)



GDTlung_s@meta.data %>% filter(cdr3_TRD == shared)

GDTlung_s@meta.data %>% filter(cdr3_TRD == 'CALGELGDDKLIF')


DVdata_D %>% filter(cdr3aa == shared)

saveRDS(DVdata_D, '/home/big/tanlikai/Fetal_thymus_TCR_David_V/delta.rds')




# FIgure 1, new , only gdT  -----------------------------------------------


# Fi1

F1A <-  NA

Sortting_s <- xlsx::read.xlsx(file = "LungTsorting .xlsx", sheetIndex = 1) %T>% print() 

Sortting <- xlsx::read.xlsx(file = "LungTsorting .xlsx", sheetIndex = 1) %>%
  rename("Total_T" = "Total_T_of_CD45") %>% 
  reshape2::melt(
    id.vars = c("donor", "Tissue"), 
    value.name = "Percent")  %>% 
  mutate(Percent = Percent * 100) %>% 
  rename("variable" = "Cell_Type")%T>% 
  print() 

Sortting$Percent %>% unique


pv_wilcox <- Sortting %>%   
  group_by(Cell_Type) %>% 
  wilcox_test(Percent ~ Tissue, paired = T,conf.level = 0.95, exact = TRUE, p.adjust.method = "fdr") %T>% print() 

Sortting %>%   
  group_by(Cell_Type) %>%
  wilcox_effsize(Percent ~ Tissue)
Sortting
new_labels <- c(CD4_of_CD45 = "CD4", CD8_of_CD45 = "CD8", gdT_of_CD45= "gdT", Total_T_of_CD45 = 'Total T')

F1B <- ggplot(Sortting %>%  filter(grepl("of_CD45", Cell_Type)), aes(x = Tissue, y = Percent, color = Tissue, group = donor)) +geom_point(size = 1)+
  geom_line(color = "grey",linewidth = 0.2)+ xlab(NULL)+ylab("Percentage of CD45")+
  facet_grid(~Cell_Type, shrink = T, labeller = labeller(Cell_Type = new_labels))+color_m(color = c("blue", "red"))+
  theme_minimal()+mytheme+NoLegend()+
  stat_compare_means(aes(group = Tissue), label = "p", method = "wilcox", paired = T, size = gs(8))
F1B
F1AB <- PG(list(NA, F1B), labels = "AUTO", rw = c(1.5,1)) %T>% print()

F1C <- (Feature_rast(GDTlung_s, sz = 0.3,
                     labelsize = 6,
                     noaxis = F, othertheme = list(theme(
                       legend.margin = margin(0,0,0,-10, "pt")), coord_fixed()) 
                     
)+ggtitle('gdT cells')) %T>% print() 
F1C


F1D15 <-  map( c("CD103.protein", "CD49a.protein", "CD45RA.protein", "CD27.protein" , "IL7R.protein" ), ~ 
                 Feature_rast(GDTlung_cite, ., colorgrd = "grd2",     sz = 0.2,othertheme = list(theme(
                   legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())   )             
)

F1D6 <-  Feature_rast(GDTlung_cite, "KLRG1.protein", colorgrd = "grd2",  ncol = 3, sz = 0.2,othertheme = list(theme(
  legend.margin = margin(0,0,0,-10, "pt")),
  scale_m("grd2", c(0,0.5, 1,1.5), c(0,1.5)),
  coord_fixed())   )

F1D15[[6]] <- F1D6

F1D <- PG(F1D15, ncol = 3)
F1D

F1CD <- PG(list(F1C, F1D), rw =c(1, 1.5),labels = c("C", "D") ) %T>% print()




F1E1 <- ( Feature_rast(GDTlung_s, "ID",  noaxis = F,
                      # legendcol = 7,
                     colorset = ID_cl, do.label = F, sz = 0.3, othertheme = list(theme(legend.title = element_text(size = 8),
                       legend.margin = margin(0,-10,0,0, "pt")), 
                       coord_fixed()) )+ggtitle('Sample source')+
           gglp(p = "b")+  
           scale_color_manual(values = ID_cl,
         labels = rep(paste0('p', c(25,27,32,45,71,73,77)),2) )+ 
           guides(color = guide_legend(nrow = 2, title = 'LN\n\nLung\n\npatient',  byrow = T, label.position = 'bottom',override.aes =list(size = 2)),  fill = F)
           
           )  %T>% print()


F1E2 <-   (ggplot(cl_comp_gd, aes(y = n, x = Cell_cluster, fill = ID, 
                                       color = ID,
                                       stratum = ID  )) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_stratum(alpha=0.8 ,size = 0.2)+
  scale_fill_manual(values = ID_cl,
                    labels = rep(paste0('p', c(25,27,32,45,71,73,77)),2)
  )+
  scale_color_manual(values = ID_cl)+
  theme_minimal() + 
  ylab('Cell no. of\ncluster per donor')+
  # xlab('cluster')+ 
  xlab(NULL)+
  ggtitle("gd T cell tissue distribution")+
  scale_y_continuous(labels = abs)+
    mytheme+
    theme(
          axis.text.x = element_text(size = 8, angle = 90),
          axis.line = element_blank())+
  geom_hline(yintercept = 0,size = 0.5)+NoLegend())%T>% print()

F1E <-  PG(list(F1E1, F1E2), rw = c(1, 1.1) ) %T>% print()




F1F  <-  Feature_rast(GDTlung_s, "v_gene_TRD",  noaxis = F, navalue = alpha("lightgrey", 0.5),
                      do.label = F, sz = 0.3, 
                      othertheme = list(theme(
                        legend.margin = margin(0,0,0,-15, "pt")), coord_fixed()) )  %T>% print()



?plot_grid

F1EF <-  PG(list(F1E, F1F), ncol = 2, rw = c(1.8, 1), labels = c("E", "F")) %T>% print()

RColorBrewer::brewer.pal(n = 7,name = "Accent")
cols.use <- list(Cell_cluster = umap.colors, tissue = c("blue", "red"), patient =RColorBrewer::brewer.pal(n = 7,name = "Accent")
  )

F1G <-  DoMultiBarHeatmap(subset(GDTlung_s,downsample = 500),features = top10_hgdT$gene,  group.by='Cell_cluster', size = gs(7),additional.group.sort.by = "tissue", angle = 20,
                  additional.group.by = c("tissue", "patient"),
                  cols.use =cols.use) %>% heat_theme(legend.position = "right") %T>% print()


F1_new <-  PG(list(F1AB, F1CD, F1EF,F1G), ncol = 1, rh = c(1, 1.1, 1.2, 3.3), labels = c(NA,NA, NA, "G")) %T>% figsave(path = figpath_ni, "Figure1_test_onlygdT.pdf", 195, 290 )


Feature_rast(GDTlung_s, "cdr3_TRD_perc",  noaxis = F, navalue = alpha("lightgrey", 0.5),
             do.label = F, sz = 0.3, 
             othertheme = list( ggtitle("TRD expansion"),
               theme( 
               legend.margin = margin(0,0,0,-15, "pt")), coord_fixed()) )  %T>% print()


patients <-  GDTlung_s$patient %>% unique()


map(patients,  ~ Feature_rast(GDTlung_s %>%  subset(patient == .x), "cdr3_TRD_perc",  noaxis = F, navalue = alpha("lightgrey", 0.5),
                              do.label = F, sz = 0.5, 
                              othertheme = list(
                               theme( 
                                legend.margin = margin(0,0,0,-15, "pt")), coord_fixed()) )  %T>% print() ) %>% PG()



# Fig4 --------------------------------------------------------------------

Feature_rast(GDTlung_s, "AREG", facets= "patient")

F4A <- (Feature_rast(GDTlung_s, sz = 0.5, "v_gene_TRD",
                  labelsize = 8,  do.label = F, othertheme =list(theme(
                    legend.margin = margin(0,0,0,-10, "pt")),coord_fixed() ),
                  noaxis = F
                  
)+ggtitle('TRDV genes')) %T>% print() 


VD2freqtable %<>% mutate(tissue = replace(tissue, tissue == "Pulm", "Lung"))


F4B <-  (VD2freqtable %>% 
           ggplot(
             aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD  , alluvium  = cdr3_TRD            ))+
           ggtitle("Vd2 TCRs sharing between lung & LN")+
           geom_flow(stat = "alluvium",
                     color = "darkgray") +
           # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
           theme_minimal_hgrid()+
           # scale_fill_manual(values = rainbow(270))+
           geom_stratum()+ 
           xlab(NULL) +ylab("TCRD frequencies")+
           theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
           guides(fill = guide_legend(ncol = 1, title = NULL)) +mytheme)%>% rasterise(dpi = 300) 






F4B

F4C <- Feature_rast(GDTlung_s, 
                    c("SELL",  "KLRG1", "S1PR1", "KLF2"),
                    sz = 0.3, ncol = 2,  othertheme =list(theme(
                      legend.margin = margin(0,0,0,-10, "pt")),coord_fixed() )  ) %T>% print() 

F4D  <- Feature_rast(GDTlung_cite, 
                    c("KLRB1.protein",  "CCR6.protein", "CD26.protein", "IL7R.protein"),
                    sz = 0.5, ncol = 2, colorgrd = "grd2",
                    othertheme =list(theme(
                      legend.margin = margin(0,0,0,-10, "pt")),coord_fixed() )  ) %T>% print() 
F4D

F4E  <- Feature_rast(GDTlung_s, 
                     c("RORC", "CCR6",  "KLRB1", "DPP4"),
                     sz = 0.5, ncol = 2,  othertheme =list(theme(
                       legend.margin = margin(0,0,0,-10, "pt")),coord_fixed() ) ) %T>% print() 
Feature_rast(GDTlung_s, 
             c("RORC", "CCR6",  "IL23R", "KLRG1"),
             sz = 0.5, ncol = 2,  othertheme =list(theme(
               legend.margin = margin(0,0,0,-10, "pt")),coord_fixed() ) ) 






F4F <-  (ViolinPlot(GDTlung_s, "GM_D", box= T, jitter = F,  x.angle = 90, 
                    colors = umap.colors, ylabtext = "score", othertheme = theme(plot.title = element_text(size = 8)))+
           ggtitle("Innate Type-3 module"))%T>% print() 

F4ABC <- PG(list(F4A, F4B,F4C), labels = "AUTO", ncol = 3, rw = c(1.2,0.8,1))
F4ABC

F4DEF <- PG(list(F4D, F4E,F4F), labels =  c("D", "E", "F"), ncol = 3, rw = c(1,1,1))


F4G <- (Feature_rast(GDTfacs, g = "TCRtype", d1 = "Vd2", d2 = "Vg9", facets = "tissue", sz = 0.3,
                     do.label = F, facetcol = 1,
                     othertheme = list(coord_fixed()))+
          geom_hline(yintercept = 2.5, linetype  = 'dotted')+geom_vline(xintercept = 3.2, linetype = 'dotted')+
          geom_text(data = vg9vd2st, aes(label =paste0(mean,"%±",sd), x = Inf, y = Inf), size = gs(8),
                    hjust = 1.1, vjust = 1.1)+
          theme_bw()+mytheme+theme(panel.grid = element_line(size = 0) )) %T>% print()

F4H1 <- ggplot(GDTfacs , aes_string(x = "CD103", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + ylab("density")+
  # geom_text(data = stTRM[1,], aes(label = paste0("CD103+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 2.2) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  # geom_text(data = stTRM[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 1.8) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme+NoLegend()+
  theme(panel.background = element_blank())

F4H2 <- ggplot(GDTfacs , aes_string(x = "KLRG1", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + ylab("density")+
  # geom_text(data = stKLRG1[1,], aes(label = paste0("KLRG1+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 1.4) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  # geom_text(data = stKLRG1[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 1.2) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme+NoLegend()+
  theme(panel.background = element_blank())

F4H3 <- ggplot(GDTfacs , aes_string(x = "CD26", fill = "TCRtype"))+geom_density(alpha = 0.5)+theme_bw(base_line_size = 0) + ylab("density")+
  # geom_text(data = stCD26[1,], aes(label = paste0("CD16+\nNonVg9Vg2:", mean,"%±",sd),x = Inf, y = 0.7) , hjust = 1.1, vjust = 1.1, size = gs(8))+
  # geom_text(data = stCD26[2,], aes(label = paste0("Vg9Vg2:", mean,"%±",sd),x = Inf, y = 0.5) , hjust = 1.1, vjust = 1.1, size = gs(8) )+
  geom_vline(xintercept = 3.2, linetype = 'dotted')+fill_m()+mytheme+NoLegend()+
  theme(panel.background = element_blank())




F4Hb <- c("PercentCD103", "PercentKLRG1", "PercentCD26")  %>% 
  map(~ ggplot(TRMKLRG1CD26_inVg9Vd2 , aes_string(x = "TCRtype", y = .x ))+geom_point(aes(color = TCRtype))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
        stat_compare_means(  paired = T, method = 't.test')
      
  ) %>% PG(ncol = 3)
+TRMKLRG1CD26_inVg9Vd2[c(1,5,2,6,3,7,4,8),]
compare_means(PercentCD103      ~ TCRtype   ,  data = TRMKLRG1CD26_inVg9Vd2[c(1,6,2,5,3,7,4,8),],
              paired = T,
              method = "t.test")


F4H <-  PG(list(F4H1, F4H2, F4H3), ncol = 3) %>%  list(F4Hb) %>%  PG(ncol = 1, rh = c(1,1)) %T>% print()



F4GH <- PG(list(F4G,F4H), labels = c("G", "H"), nrow = 1, rw= c(1.1,2)) %T>% print()







F4 <- (PG(list(F4ABC,F4DEF,F4GH), nrow = 3, greed = F, rh = c(1,1,1)) +
         draw_figure_label('Figure 4',
                           
                           position = 'top.right', size = 10, fontface = 'plain'))%T>% 
  print() %>% figsave( "Fig4_Vg9Vd2_2024_withFACS_stat.pdf", 
                       path = figpath_ni,
                       200, 210)


F4ABC
PG(list(F4ABC,F4DEF,F4G_J), nrow = 3, greed = F, rh = c(1,1,0.6))

# Figure 5 ----------------------------------------------------------------


# F5A morisita
write.csv(Mori_result_TRD, "Mori_result_TRD.csv")

Mori_result_TRD <- read.csv("Mori_result_TRD.csv")
data_long <- melt(Mori_result_TRD, id.vars = "Unnamed: 0")


data_long  %<>% mutate(Category.Var1 = case_when(
  Var1  == "LN1" ~ "Naive",
  Var1 %in% c("LN2", "LN3", "LN4") ~ "LN_TEM",
  Var1 %in% c("LG1", "LG2", "LG3", "LG4") ~ "Lung_TEMRA",
  Var1 %in% c("LG5", "LG6", "LG7") ~ "Lung_TRM"
),
Category.Var2 = case_when(
  Var2  == "LN1" ~ "Naive",
  Var2 %in% c("LN2", "LN3", "LN4") ~ "LN_TEM",
  Var2 %in% c("LG1", "LG2", "LG3", "LG4") ~ "Lung_TEMRA",
  Var2 %in% c("LG5", "LG6", "LG7") ~ "Lung_TRM"
)
)




F5C <-  ggplot(data_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + theme_nothing()+
  scale_fill_gradientn( na.value = 'white',
                        colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  mytheme+gglp("r")+
  theme(axis.text.x = element_text(size = 8 , angle  = 90))+
  xlab(NULL)+ylab(NULL)+theme_minimal()+
  # theme(text = element_text(size = 2))+ 
  ggtitle('Morisita index of TRD TCRs') +coord_fixed()+mytheme+heattheme

F5C
library(ComplexHeatmap)
library(circlize)

categories <- ifelse(colnames(Mori_result_TRD) %in% c("LN1", "LN2", "LN3", "LN4"), "LN_TEM",
                     ifelse(colnames(Mori_result_TRD) %in% c("LG1", "LG2", "LG3", "LG4"), "Lung_TEMRA",
                            "Lung_TRM"))
# col_map <- c("LN_TEM" = "red", "Lung_TEMRA" = "green", "Lung_TRM" = "blue")
col_map <- c("LN_TEM" = umap.colors[1], "Lung_TEMRA" = umap.colors[3], "Lung_TRM" = umap.colors[6])


# Create a HeatmapAnnotation object for the color bar
ha <- HeatmapAnnotation(df = data.frame(Category = factor(categories, levels = names(col_map))),
                        show_legend = F,
                        col = list(Category = col_map),
                        show_annotation_name = F)
ha
# Create the heatmap
Heatmap(Mori_result_TRD,
        name = "Morisita Index", 
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 12, name = "Spectral")))(100), 
        na_col = "white",
        show_row_names = T,cluster_columns = F, cluster_rows = F, 
        # border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize =8),
        heatmap_legend_param = list(
          heatmap_legend_direction = "horizontal",
          legend_direction = "horizontal",
          heatmap_legend_side="bottom",
          annotation_legend_side="right",
          legend_grouping = "original",
          legend_width = unit(2, "cm")), 
        column_names_gp = gpar(fontsize = 8),
        top_annotation = ha
)


F5A <-  Heatmap(Mori_result_TRD,
        name = "Morisita Index", 
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100), na_col = "white",
        show_row_names = T,cluster_columns = F, cluster_rows = F, 
        # border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize =8),
        heatmap_legend_param = list(
          heatmap_legend_direction = "horizontal",
          legend_direction = "horizontal",
          heatmap_legend_side="bottom",
          annotation_legend_side="right",
          legend_grouping = "original",
          legend_width = unit(2, "cm")), 
        column_names_gp = gpar(fontsize = 8),
        top_annotation = ha
            )   %>% ggplotify::as.ggplot()
F5A+coord_fixed()
F5A1 <- draw(F5A, heatmap_legend_side="bottom", annotation_legend_side="right",
     legend_grouping = "original") %>% grid.grabExpr() %>% ggplotify::as.ggplot()
F5A1
TCRfreqtable$TCRfreq



LULN_TCRgd <- TCRfreqtable %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  arrange(desc(TCRfreq)) %>% 
  # mutate(cdr3_paired = factor(cdr3_paired, levels = unique(cdr3_paired))) %>% 
  
  filter(pt %in% c('Lung_TRM_Pulm', 'LN_memory_LN', 'Lung_memory_Pulm'))



write.csv(LULN_TCRgd, "LULN_TCRgd.csv")

F5D <-   TRDfreqtable %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  
  filter(pt %in% c('Lung_TRM_Lung', 'LN_memory_LLN', 'Lung_memory_Lung')) %>%  
  
  mutate(pheno = factor(pheno, levels = rev(c('Lung_TRM','Lung_memory',  'LN_memory')))) %>% 
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_TRD, 
         stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("Non-VD2 TCR sharing")+
  geom_flow(stat = "alluvium", size =0.2,
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = ggplotColours(211))+
  geom_stratum(size = 0.1)+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"), 
        axis.text.x = element_text(angle = 90) )+
  guides(fill = guide_legend(ncol = 1, title = NULL))+mytheme+
  # facet_wrap(~patient)+
  NULL

F5D

library(circlize)
library(RColorBrewer)
library(igraph)
library(ggraph)


TCRgddata <- TCRfreqtable %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  arrange(desc(TCRfreq)) %>% 
  # mutate(cdr3_paired = factor(cdr3_paired, levels = unique(cdr3_paired))) %>% 
  
  filter(pt %in% c('Lung_TRM_Lung', 'LN_memory_LLN', 'Lung_memory_Lung'))


pheno_pairs <- unique(TCRgddata$pheno)
shared_matrix <- matrix(0, nrow = length(pheno_pairs), ncol = length(pheno_pairs),
                        dimnames = list(pheno_pairs, pheno_pairs))
ribbon_colors <- list()

cdr3_paired_colors <- setNames(rainbow(length(unique(TCRgddata$cdr3_paired))), unique(TCRgddata$cdr3_paired))


# Compute the sum of shared TCRfreq for each pair of phenos
for (i in seq_along(pheno_pairs)) {
  for (j in seq_along(pheno_pairs)) {
    if (i != j) {
      # Get the cdr3_paired and TCRfreq for the two phenos
      data_i <- TCRgddata[TCRgddata$pheno == pheno_pairs[i], c("cdr3_paired", "TCRfreq")]
      data_j <- TCRgddata[TCRgddata$pheno == pheno_pairs[j], c("cdr3_paired", "TCRfreq")]
      
      # Compute the sum of shared TCRfreq, and the colors for the ribbons
      shared_cdr3 <- intersect(data_i$cdr3_paired, data_j$cdr3_paired)
      shared_matrix[i, j] <- sum(data_i[data_i$cdr3_paired %in% shared_cdr3, "TCRfreq"]) +
        sum(data_j[data_j$cdr3_paired %in% shared_cdr3, "TCRfreq"])
      ribbon_colors[[paste(pheno_pairs[i], pheno_pairs[j], sep = "-")]] <- cdr3_paired_colors[shared_cdr3]
    }
  }
}

# Create the chord diagram
chordDiagram(shared_matrix, col = unlist(ribbon_colors))

edges <- TCRgddata[, c("cdr3_paired", "pheno")]
colnames(edges) <- c("from", "to")

# Create a graph from the edge list
g <- graph_from_data_frame(edges, directed=FALSE)
g
# Plot the graph
ggraph(g, layout = 'fr') + 
  geom_edge_link() + 
  geom_node_point(size = 5, aes(color = name )) + 
  # geom_node_text(aes(color = name), vjust = 1, hjust = 1) +
  theme_graph()+NoLegend()


F5E <- Feature_rast(GDTlung_s, 
                    c(  "ZNF683", "BATF","GATA3","EOMES"),
                    sz = 0.3, ncol = 2,  othertheme =list(theme(
                      legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())  ) %T>% print() 

F5AB <- NA


Feature_rast(GDTlung_s, 
                    c(  "ZNF683", "BATF","GATA3","EOMES"),
                    sz = 0.3, ncol = 2,  othertheme =list(theme(
                      legend.margin = margin(0,0,0,-10, "pt")), coord_fixed())  ) %T>% print() 


F5ABE <-  PG(list(F5AB, F5E), ncol = 1, rh = c(1.5,2))

F5ABE




# 
# 
# 
# F5F <- GSEA_multipplot(GSEA_TRMvsTEMRA_allref,description_to_show = c('GOBP_DEFENSE_RESPONSE_TO_BACTERIUM', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY'), base_size = 8,legend.position = 'bottom',col = umap.colors[c(9,4)],
#                        title = " Cytotoxictiy",
#                        
#                        c1='TRM_LG6', c2 = 'TEMRA_LG3' )
# F5D
# 
# 
# F5E <- GSEA_multipplot(GSEA_TRMvsTEMRA_allref,description_to_show = c(
#   'GOBP_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION',
#   'GSE7852_TREG_VS_TCONV_UP'
#   
# ), title = "GSEA regulatory",
#  
# base_size = 8, col = umap.colors[c(6,7)],
# legend.position = 'bottom',
# c1='TRM_LG6', c2 = 'TEMRA_LG3' )
Feature_rast(GDTlung_s, "BATF")



F5F <- NA

F5G <- Feature_rast(Reggd, c("BATF_REG", "EOMES_REG", "GATA3_REG"), colorgrd = "grd2", sz = 0.3, othertheme =list(theme(
  legend.margin = margin(0,0,0,-10, "pt")), coord_fixed()) ) %T>%  print()


ViolinPlot()



F5FG

Feature_rast(Reggd, c("BATF_REG", "EOMES_REG", "GATA3_REG"))


ggplot(Reggd, aes(x = Cell_cluster, y = GATA3_REG, fill = Cell_cluster))+geom_violin()+geom_boxplot()
GSEA_TRMvsTEMRA_allref@result %>% view

F5H <-  GSEA_multipplot(GSEA_TRMvsTEMRA_allref, c(
                                          'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY' ,
                                          "GOBP_PHAGOCYTOSIS",
                                          "GOBP_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION",
                                          'GSE7852_TREG_VS_TCONV_UP'),title = 'GSEA: gd_TRM_LG6 vs gd_TEMRA_LG3', legendpvalue = T, legend.position = "no",
                            plots = 1:2, rel_h = c(1,0.25),
                            
                            c1 = "gdTEMRA_LG3", c2 =  "gdTRM_LG6", base_size = 8 ) %T>% print()
F5H <-   PG(list(F5H, NA), ncol =2, rw = c(1, 1)) %T>% print()

figsave(F5H, "F5H.pdf", 150, 70)


F5FGH <- PG(list(F5F, F5G, F5H), ncol = 1, labels = c("F", 'G', "H"), rh = c(1.3,1,1.3))%T>% print()

F5FGH

F5I <-Feature_rast(GDTlung_s, 
                   c( "GZMA","CSF1", "CTLA4", "ENTPD1"),
                   sz = 0.3, ncol = 2,  othertheme =list(theme(
                     legend.margin = margin(0,0,0,-10, "pt")), coord_fixed()) ) %T>% print() 



Feature_rast(GDTlung_s, 
             c( "AREG","CSF1", "CTLA4", "ENTPD1"),
             sz = 0.3, ncol = 2,  othertheme =list(theme(
               legend.margin = margin(0,0,0,-10, "pt")), coord_fixed()) ) %T>% print() 



F5J <-  Feature_rast(GDTfacs %>%  filter(gdTtype == "NonVd2"), "TRM", d1 = "CD103", d2 = "CD49a", noaxis = F, axis.number = T , do.label = F, sz = 0.3,  facets = "tissue",
                     othertheme = list(coord_fixed()))  +theme_bw(base_line_size = 0)+
  ggtitle("Tissue residency")+mytheme

F5J


F5K1 <- map2(c("GATA3",  "CD357",  "GzmA","Eomes"), list(2.5,3,2.5,NULL) ,~
              ggplot(GDTfacs_lung_Vd1, aes_string(x = .x, fill = "TRM"))+geom_density(alpha=0.4)+
              theme_classic()+mytheme+fill_m()+NoLegend() +geom_vline(xintercept = .y, linetype = 'dotted')+
              theme(panel.background = element_blank()) ) %>% PG(ncol = 4, labels = 'Lung Vd1/3 T cells')


F5K2 <-  c("Gata3Pos", "CD357Pos", "GzmAPos", "EomesMFI_mean")  %>% 
  map(~ ggplot(Vd1Lungmarkers , aes_string(x = "TRM", y = .x ))+geom_point(aes(color = TRM   ))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
        stat_compare_means(  paired = T, method = 't.test',label = "p.signif" )
      
  ) %>% PG(ncol = 4)



c( "EomesMFI_mean","Gata3Pos")  %>% 
  map(~ ggplot(Vd1Lungmarkers , aes_string(x = "TRM", y = .x ))+geom_point(aes(color = TRM   ))+geom_line(linewidth = 0.3, aes(group = batch))+color_m()+theme_classic()+mytheme+NoLegend()+
        stat_compare_means(  paired = T, method = 't.test' ,label = "p")
      
  ) %>% PG(ncol = 2)

F5K <-  PG(list(F5K1, F5K2), ncol = 1) 

F5K


F5ABEIJ <-  PG(list(F5AB, F5E,F5I,F5J), ncol = 1, rh = c(1.5,2,2,1), labels = c("A", "E", "I", "J"))

F5ABEIJ


F5CD <- PG(list(F5C,F5D), labels = c("C", "D" ), ncol = 2, rw = c(1.1,0.8)) %T>%  print()


F5FGHK <- PG(list(F5F, F5G, F5H,F5K), ncol = 1, labels = c("F", 'G', "H", "K"), rh = c(1.3,1,1.3, 1.3))%T>% print()

F5C_H <- PG(list(F5CD,F5FGHK), ncol = 1, rh = c(1.3,5))%T>% print()

F5_test <-  PG(list(F5ABEIJ, F5C_H), ncol = 2, rw = c(1,2.5))%T>% figsave("Fig5_2024_test_FACS_Scenic.pdf", 200, 280, path = figpath_ni)



F5ABC <- PG(list(F5A1, F5B,F5C), labels = c("C", "D" ,"E"), ncol = 3, rw = c(1.1,0.8,0.9)) %T>%  print()


F5DEF <- PG(list(F5D, F5E,F5F), labels = c("F","G", "G"),
            greed = T, scale = 0.85,
            ncol = 3, rw = c(0.95,0.95,0.9)) %>% 
  list(NA) %>% PG(ncol = 1, rh = c(1, 0.05))%T>%  print()  

F5
F5 <- (PG(list(F5ABC,F5DEF), nrow = 2, greed = T, scale = 0.85) +
  draw_figure_label('Figure 5',  
                    
                    position = 'top.right', size = 10, fontface = 'plain'))%T>% 
  print() %T>% figsave("Fig5_2024.pdf", 200, 150, path = figpath_ni)


TCRfreqtable_TRD <- GDTlung_s@meta.data %>% filter(cdr3_TRD_freq >1  & paired != 'GV9 DV2'    )  %>% group_by(pheno,tissue, cdr3_TRD) %>%  summarise(TCRfreq = n()) %>% 
  arrange(TCRfreq) #%>% mutate(cdr3_paired=factor(cdr3_paired, levels = unique(cdr3_paired)))
TCRfreqtable_TRD


dTCRsharing <- TCRfreqtable_TRD %>% mutate(pt = paste0(pheno, "_", tissue)) %>% 
  
  filter(pt %in% c('Lung_TRM_Pulm', 'LN_memory_LN', 'Lung_memory_Pulm')) %>%
  
  mutate(pheno = factor(pheno, levels = c('Lung_TRM', 'LN_memory', 'Lung_memory'))) %>% 
  ggplot(
    aes( x = pheno, y = TCRfreq, fill = cdr3_TRD, 
         stratum= cdr3_TRD, alluvium  = cdr3_TRD))+
  ggtitle("paired gd TCR sharing")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = ggplotColours(211))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCR frequencies")+
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 315),
        legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))+
  # facet_wrap(~patient)+
  NULL
dTCRsharing


# color set test ----------------------------------------------------------



fC <-  Feature_rast(
  GDTlung_s, 
  c("RORC", "CCR6",  "IL23R", "KLF2"),             color_grd = "D",
             sz = 0.5, ncol = 2,  othertheme =theme(
               legend.margin = margin(0,0,0,-10, "pt"))  ) %T>% print() 

fA <-  Feature_rast(GDTlung_s, 
                    c("RORC", "CCR6",  "IL23R", "KLF2"),
                    color_grd = "A",
                    sz = 0.5, ncol = 2,  othertheme =theme(
                      legend.margin = margin(0,0,0,-10, "pt"))  ) %T>% print() 

fB <-  Feature_rast(GDTlung_s, 
                    c("RORC", "CCR6",  "IL23R", "KLF2"),
                    color_grd = "C",
                    sz = 0.5, ncol = 2,  othertheme =theme(
                      legend.margin = margin(0,0,0,-10, "pt"))  ) %T>% print() 


fD <-  Feature_rast(GDTlung_s, 
                    c("RORC", "CCR6",  "IL23R", "KLF2"),
                    color_grd = "threecolor",
                    sz = 0.5, ncol = 2,  othertheme =theme(
                      legend.margin = margin(0,0,0,-10, "pt"))  ) %T>% print() 

fE <-  Feature_rast(GDTlung_s, 
                    c("RORC", "CCR6",  "IL23R", "KLF2"),
                  
                    sz = 0.5, ncol = 2,  othertheme =  scale_colour_gradient2(
                      low = "blue",
                      mid = "white",
                      high ="red",
                      # midpoint = 1,
                      space = "Lab",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "colour"
                    ) ) %T>% print() 


PG(list(fA,fB,fC,fD,fE), labels = "AUTO", ncol = 2) %>% figsave("Gene_five_color_plan.pdf", 180, 270)



Feature_rast(GDTlung_s, 
             c("RORC", "CCR6",  "IL23R", "IL4I1"),
             sz = 0.5, ncol = 2,  color_grd =  "threecolor")  %T>% print() 





Feature_rast(GDTlung_s, 
             c("RORC", "CCR6",  "IL23R", "KLF2"),
             
             sz = 0.5, ncol = 2,  othertheme =  scale_colour_gradientn(
               colors = c("#eff4ff", "#ffcc66", "red", "#990000"),
               # midpoint = 1,
               space = "Lab",
               na.value = "grey50",
               guide = "colourbar",
               aesthetics = "colour"
             ) ) 



Feature_rast(GDTlung_s, facets = "patient") %T>%  figsave("GDT_seven_pt.pdf", 400, 400)


Feature_rast(CD4CD8, facets = "patient") %T>%  figsave("CD4CD8_seven_pt.pdf", 400, 400)


data.table::fwrite(GDTlung_s[["RNA"]]@data %>% as.data.frame(), "Normalized_count_gdtlung.csv", row.names = T)

 data.table::fwrite(data.frame(FetchData(GDTlung_s, c("UMAP_1", "UMAP_2")), GDTlung_s@meta.data), 
          "GDTlung_meta.csv", row.names = T)


allMetagdt <- data.frame(FetchData(GDTlung_s, c("UMAP_1", "UMAP_2")), GDTlung_s@meta.data) %>%  filter(!is.na(cdr3_TRD) & Cell_cluster != "Vg9Vd2_M")


allMetagdt %>% Feature_rast("v_gene_TRD", sz = 1.5)

allMetagdt$table()


allMetagdt$v_gene_TRD %>% table()



GDTlung_s$Vg9Vd2

Feature_rast(GDTlung_s, "Vg9Vd2")




colnames(GDTlung_s@meta.data)

Feature_rast(GDTlung_s, "GM_D",
             colorgrd =  c("#eff4ff", "#eff4ff", "purple", "#990000"),
             navalue ="transparent" 
             )



dim(AUCmat)





VD1freqtable <- GDTlung_s@meta.data %>%
  dplyr::filter(cdr3_TRD_freq > 1 &   v_gene_TRD == 'TRDV1' )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
Vd1s <-  VD1freqtable %>% 
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD  , alluvium  = cdr3_TRD            ))+
  ggtitle("Vd1 TCR  sharing between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))





VD2freqtable <- GDTlung_s@meta.data %>%
  dplyr::filter(cdr3_TRD_freq > 1 &   v_gene_TRD == 'TRDV2' &
                  Cell_cluster != "Vg9Vd2_M" )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
Vd2s <-  VD2freqtable %>% 
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD  , alluvium  = cdr3_TRD            ))+
  ggtitle("Vd2 TCR  sharing between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))
Vd2s


VD3freqtable <- GDTlung_s@meta.data %>%
  dplyr::filter(cdr3_TRD_freq > 1 &   v_gene_TRD == 'TRDV3' )  %>%
  group_by(tissue, cdr3_TRD) %>%  summarise(TCRDfreq = n()) %>% 
  arrange(TCRDfreq) 
Vd3s <-  VD3freqtable %>% 
  ggplot(
    aes( x = tissue, y = TCRDfreq, fill = cdr3_TRD,  stratum= cdr3_TRD  , alluvium  = cdr3_TRD            ))+
  ggtitle("Vd3 TCR  sharing between lung & LN")+
  geom_flow(stat = "alluvium",
            color = "darkgray") +
  # scale_y_continuous(limits = c(0, 40), breaks = c(0, 10,20,30,40))+
  theme_minimal_hgrid()+
  # scale_fill_manual(values = rainbow(270))+
  geom_stratum()+ 
  xlab(NULL) +ylab("TCRD frequencies")+
  theme(legend.position = 'none', legend.key.size = unit(2, "mm"))+
  guides(fill = guide_legend(ncol = 1, title = NULL))




PG(list(Vd1s,Vd2s,Vd3s), ncol = 3)





TCRcompo <-   GDTlung_s@meta.data  %>%  filter(!is.na(cdr3_paired)) %>% 
  group_by(Cell_cluster) %>%  count(paired) %>% 
  mutate(pct = n /sum(n)*100) 

ggplot(TCRcompo, aes(x = Cell_cluster, y = pct, fill = paired))+
  geom_col()


vg9vd2p <-  GDTlung_s@meta.data %>% filter(paired == "GV9 DV2") %>% count(Cell_cluster)  %>% 
  mutate(pct = n /sum(n)*100) 
levels(GDTlung_s$Cell_cluster)

Feature_rast(GDTlung_s, "Cell_cluster")

vg9vd2p$Cell_cluster <-  factor(vg9vd2p$Cell_cluster, levels = levels(GDTlung_s$Cell_cluster) )

ggplot(vg9vd2p, aes(x="", y=pct, fill=Cell_cluster)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+fill_m(color = umap.colors[c(1:8,10,11,9)])+ mytheme+ggtitle("Vg9Vd2 TCR distribution")+theme_minimal()+
  theme(axis.title = element_blank())


Feature_rast(GDTlung_cite, "KLRG1.protein", d1 = "KLRG1.protein", d2 = "CD103.protein", noaxis = F, axis.number = T)


Feature_rast(GDTlung_cite, "KLRG1.protein")
