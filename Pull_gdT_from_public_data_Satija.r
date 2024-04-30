# public lung satija 
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

source('/home/big/tanlikai/script/rscripts/funcs.r')
setwd('/home/big/tanlikai/Lung/')

# readdata
Lung_satija <- readRDS('public/Azimuth.lung.Satija.rds')

dim(Lung_satija)

Lung_satija@assays$RNA@key <- 'rna_'
glimpse(Lung_satija@meta.data)
Lung_satija$annotation.l1 %>% unique()

Lung_satija$dataset_origin %>% unique()



Feature_rast(Lung_satija,c('CD3D', 'TRDC', 'CD3G', 'CD3E'))

Feature_rast(Lung_satija, 'annotation.l1', colorset = 'gg')

Feature_rast(Lung_satija, d1 = 'TRDC', d2 = 'TRBC1')


Feature_rast(Lung_satija)
DimPlot(Lung_satija)
Lung_satija@reductions
# FetchData(Lung_satija, c('CD3D', 'CD4', 'PTPRC'), slot = 'RNA')

Lung_satija <- AddModuleScore(Lung_satija,list(c('CD3D', 'CD3G','CD3E')),assay = 'RNA', name = 'T_cell_score')



ViolinPlot(Lung_satija, c('CD3D', 'CD3G','CD3E', 'PTPRC','T_cell_score1'),colors = ggplotColours(20))

Feature_rast(Lung_satija,g ='PTPRC', d1 = 'CD3D', d2 = 'T_cell_score1', noaxis = F, axis.number = T)


Lung_satija$dataset_origin

ViolinPlot(Lung_satija, "T_cell_score1", group.by = "dataset_origin", color = umap.colors)

VlnPlot(Lung_satija, 'T_cell_score1', raster = T )
Feature_rast(Lung_satija, c('annotation.l1', 'T_cell_score1'), colorset = 'gg')


TCELL_adams <- subset(Lung_satija, subset = (T_cell_score1 >= 1 & PTPRC > 1.8 & dataset_origin == 'adams_2020') )
rm(Lung_satija)

dim(TCELL_adams)
Feature_rast(TCELL,c('CD3D', 'CD3E',  'CD3G', 'TRGC1', 'TRGC2', 'TRBC1', 'TRBC2', 'TRAC', 'CD4', 'CD8A', 'TRDC',"T_cell_score1"))

Feature_rast(TCELL_adams,c('ITGAE', 'ITGA1', 'CXCR6', 'KLRG1', 'KLF2', "CD4", "CD8A", "KLRB1"),sz = 0.2)
Feature_rast(TCELL_adams, d1 = 'TRDC', d2 = 'TRBC1')



# tests <- Feature_rast(TCELL,'annotation.l1', colorset = 'gg')


# ViolinPlot(TCELL, c('TRDC', 'TRGC1', 'TRGC2', 'TRAC', 'TRBC1', 'TRBC2'),colors = ggplotColours(20))


Feature_rast(TCELL_adams, d1 = 'TRDC', d2 = 'TRAC', 'annotation.l2', colorset = 'gg')
Feature_rast(TCELL_adams, d1 = 'TRDC', d2 = 'TRDV1', 'annotation.l2', colorset = 'gg')


TCELL$annotation.l1 %>%  unique()

# Feature_rast(TCELL, 'annotation.l1',colorset = 'gg')
# TCELL <- subset(TCELL, subset = annotation.l1 %in% c('CD4 T', 'CD8 T', 'Proliferating NK/T', 'Natural Killer',             'Natural Killer T' ))

# TRGD score and TRAB score

TRGD <- grep('^TRDV|^TRDC', rownames(TCELL),value = T)
TRAB <- grep('^TRAV|^TRBV|^TRAC|^TRBC', rownames(TCELL),value = T)



TCELL_adams@meta.data  %<>%  select(-starts_with('TCRGD'))
TCELL_adams <- AddModuleScore(TCELL_adams,list(TRGD),assay = 'RNA', name = 'TCRGD')
TCELL_adams <- AddModuleScore(TCELL_adams,list(TRAB),assay = 'RNA', name = 'TCRAB')
# TCELL_adams(TCELL, d2 = 'TCRGD1', d1 = 'TCRAB1', 'annotation.l2', colorset = 'gg', noaxis = F, do.label = F)
# Feature_rast(TCELL, d2 = 'TCRGD1', d1 = 'TCRAB1', c('TRGC1', 'TRGC2', TRGD,'TRBC1', 'TRBC2','TRAC'), noaxis = F, do.label = F, ncol = 3, sz = 1)
Feature_rast(TCELL_adams, d2 = 'TRDC', d1 = 'TCRAB1', 'annotation.l2', colorset = 'gg', noaxis = F, do.label = F)



TCELL_adams  %<>% NormalizeData(normalization.method = 'LogNormalize',
                        scale.factor = 10000, assay = 'RNA') %>%  ScaleData( assay = 'RNA',
                           vars.to.regress = c('donor' , 'nCount_RNA')) %>%  FindVariableFeatures(assay = 'RNA',nfeatures = 3000, selection.method = 'vst')
TCELL$dataset_origin %>%  table


TCELL_adams@assays$RNA@var.features <- TCELL_adams@assays$RNA@var.features%>%  
  str_subset('^RP|^MT|TRAV|TRBV|NO-NAME|IGL|IGH|LINC|MIR', negate = T)


TCELL_adams <- RunPCA(TCELL_adams, npcs = 100, verbose =
                    T, nfeatures.print = 40) %>%  JackStraw(dims = 50)

ElbowPlot(TCELL_adams,ndims = 100)
TCELL$dataset_origin


TCELL_adams <- RunUMAP(TCELL_adams, dims = 1:50, 
                   reduction = 'pca') %>%
  FindNeighbors(dims = c(1:50))

TCELL_adams <- FindClusters(TCELL_adams, resolution = 0.6)
Feature_rast(TCELL_adams, c("ident", 'ITGAE', "ITGA1", "CD4",
                            "CD8A", 'ZNF683', "KLRB1",
                            'DPP4', 'GZMA', 'GZMB', 'NKG7'),sz = 0.2,
             color_grd = "threecolor"
             )

Feature_rast(TCELL_adams, c("ident", "TRDC", "CD4", "CD8A", "CD8B", "ITGAE", "ITGA1",'ZNF683', 'CXCR6', 'CCR6', 'GZMA', 'GZMB', 'NKG7', 'PRF1', 'AREG', "DPP4", "KLRB1"
                            ))


Feature_rast(TCELL, c('ITGA1', 'ITGAE', 'KLRG1', 'KLF2','ITGB2', 'ZNF683', 'CXCR6', 'CD69', 'GZMA', 'GZMB', 'NKG7', 'PRF1', 'AREG'))

# select cells based on TCRgd and TCRab scores ----------------------------

library(ggiraph)
library(shiny)
library(plotly)

TCELL@assays$RNA@data
TCRabgd <- FetchData(TCELL, c('TRGC1', 'TRGC2', 'TRDC','TRBC1', 'TRBC2','TRAC','TCRGD1', 'TCRAB1'), slot = 'data')
TCRabgd$rowname <- rownames(TCRabgd)
TCRabgd$TCRs <- NA
TCRabgd$pointNumber          <- 1:nrow(TCRabgd)
rownames(TCRabgd) <- rownames(TCRabgd)
TCRabgd$rowname

TCRabgd_inte <- ggplot(TCRabgd, aes(y = TCRGD1, x = TCRAB1, color = TRAC, key = rowname,
                                tooltip = rowname, data_id = rowname))+geom_point(size = 0.3)+
  scale_color_gradientn( 
                         colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  theme_minimal()



# TCRab <-   Seurat::CellSelector(TCRabgd_inte, object =TCELL )
library(plotly)


ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brush")
)


# TCRab
server <- function(input, output, session) {
  frame2<<- NULL
  
  frame2 <<-data.frame()
  # frame2$number <<- row.names(TCRabgd)
  # nms <- TCRabgd$rowname
  
  output$plot <- renderPlotly({
    ggplotly(TCRabgd_inte) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRabgd[frame2$key, 'TCRs'] <<- 'TCRgd+'
  })
  
}

# frame2
shinyApp(ui, server)




TCRabgd$TCRs

TCRabgd[frame2$pointNumber, ]

# TCRgd
server <- function(input, output, session) {
  frame2 <<-data.frame()
  nms <- row.names(TCRabgd)
  
  output$plot <- renderPlotly({
    ggplotly(TCRabgd_inte) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRabgd[frame2$key, 'TCRs'] <<- 'TCRgd+'
  })
  
}
shinyApp(ui, server)
Feature_rast(TCRabgd, d1 = 'TCRGD1', d2 = 'TCRAB1', 'TCRs', colorset = 'gg', noaxis = F, do.label = F)


# double positive 
server <- function(input, output, session) {
  frame2<<- NULL
  frame2 <<-data.frame()
  nms <- row.names(TCRabgd)
  
  output$plot <- renderPlotly({
    ggplotly(TCRabgd_inte) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRabgd[frame2$key, 'TCRs'] <<- 'DP_TCR'
    
  })
  
}
shinyApp(ui, server)


# double negative 
server <- function(input, output, session) {
  frame2 <<-data.frame()
  nms <- row.names(TCRabgd)
  
  output$plot <- renderPlotly({
    ggplotly(TCRabgd_inte) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRabgd[frame2$key, 'TCRs'] <<- 'TCR_DN'
  })
  
}
shinyApp(ui, server)



ClusterCompare(TCELL, id1 = 'DP_TCR', id2 = 'TCRgd+', group.by = 'TCRs', features = c(TRGD),log2fc = 0,p_cutoff = 5,min.pct = 0.001)

TCRabgd$TCRs %>% table()

Feature_rast(TCRabgd, d1 = 'TCRGD1', d2 = 'TCRAB1', g = 'TCRs')
TCRabgd

TCELL$rowname <- colnames(TCELL)

TCELL@meta.data  %<>%  left_join(select(TCRabgd, TCRs, rowname), by = 'rowname', suffix = c('', ''))
TCELL@meta.data   %<>% `rownames<-`(TCELL$rowname)


Feature_rast(TCELL,c('annotation.l1', 'TCRs'), sz = 0.5, colorset = 'gg')

table(TCELL$TCRs, TCELL$annotation.l2 )

Feature_rast(TCELL)


ViolinPlot(subset(TCELL, subset = (dataset_origin == 'adams_2020')) , 'nCount_RNA', group.by = 'TCRs', box = T)+ylim(0,10000)






GDTcell <- subset(TCELL, subset = (TCRs ==  'TCRgd+'& dataset_origin == 'adams_2020'))

dim(GDTcell)
colnames(TCELL)

GDTcell$disease
Feature_rast(GDTcell_adams, c("ident", "disease"),
             do.label = F,mythe = F,
             sz = 2, ncol = 1)

Feature_rast(GDTcell_adams, c("ident"),
             do.label = F,mythe = T, colorset = set_sample(umap.colors,s = 1),
             sz = 2, noaxis = F)

Feature_rast(GDTcell_adams, c("disease"),
             do.label = F,mythe = T, colorset = set_sample(umap.colors,s = 3),
             sz = 2, noaxis = F)

GDTcell_adams$dataset_origin





Feature_rast(subset(TCELL, subset = (TCRs == 'TCRgd+' & dataset_origin == 'adams_2020')) , 'annotation.l1')

Feature_rast(GDTcell, c('ITGAE', 'KLF2'))


GDTcell$tissue %>% unique()


glimpse(GDTcell@meta.data)
GDTcell@assays$RNA@data

GDTcell %<>%  PercentageFeatureSet( '^MT', col.name =  'percent.mito') 

VlnPlot(GDTcell, c('nCount_RNA', 'percent.mito'))


Feature_rast(GDTcell, g = 'percent.mito', d1 ="nCount_RNA",d2 ='nFeature_RNA', 
             noaxis = F, axis.number = T)+grd



GDTcell<- NormalizeData(GDTcell, normalization.method = 'LogNormalize',
              scale.factor = 10000, assay = 'RNA')



GDTcell  %<>%   ScaleData( assay = 'RNA',
                           vars.to.regress = c('donor' , 'nCount_RNA', 'percent.mito')) 

GDTcell %<>%  FindVariableFeatures(assay = 'RNA',nfeatures = 3000, selection.method = 'vst')
GDTcell@assays$RNA@var.features
GDTcell@assays$RNA@var.features <- GDTcell@assays$RNA@var.features%>%  
  str_subset('^RP|^MT|TRAV|TRBV|NO-NAME|IGL|IGH|LINC|MIR', negate = T)


GDTcell <- RunPCA(GDTcell, npcs = 100, verbose =
                     T, nfeatures.print = 40)



ElbowPlot(GDTcell, ndims = 100)

GDTcell  <- JackStraw(GDTcell, num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 

JackStrawPlot(GDTcell, dims = 1:50 )


GDTcell <- RunUMAP(GDTcell, dims = 1:16, 
                    reduction = 'pca') %>%
  FindNeighbors(dims = c(1:16))

GDTcell <- FindClusters(GDTcell, resolution = 0.5)
Feature_rast(GDTcell, sz = 1)

saveRDS(GDTcell, 'public/GDTcell_Adams2020.rds')

dim(GDTcell)
# 
# GDTcell  %<>%  RunUMAP(dims = 1:16, 
#                        reduction = 'pca', min.dist = 1) 

Feature_rast(GDTcell, 'donor', colorset = 'gg')

Feature_rast(GDTcell,c( 'ITGA1', 'ITGAE', 'KLRG1', 'KLF2') ,sz = 1)

ClusterCompare(GDTcell, '0', '1')

ViolinPlot(GDTcell, 'nCount_RNA', group.by = 'donor')
ViolinPlot(GDTcell, 'nCount_RNA')


colnames(GDTcell@meta.data)

Feature_rast(GDTcell, c("ident", 'disease'))

GDTcell$donor %>% table() %>% sort(decreasing = T)
Feature_density(GDTcell,c('ITGA1', 'ITGAE', 'KLRG1', 'KLF2'))
Treglist <- c('FOXP3',  'IL2RA', 'IL2RB', 'IKZF2', 'CTLA4', 'TNFRSF18', 'TNFRSF4', 'CAPG', 'CHCHD10', 'GPR83' , 'LGALS3', 'LAG3')

GDTcell %<>%  AddModuleScore( features = list(Treglist), name = 'Treg_module', assay = 'RNA')


ViolinPlot(GDTcell, 'Treg_module1', box = T, colors = umap.colors)

VlnPlot(GDTcell, 'Treg_module1')

Feature_density(GDTcell, c('ITGA1', 'ITGAE', 'KLRG1', 'KLF2','ITGB2', 'ZNF683', 'CXCR6', 'CSF1', 'GZMA', 'GZMB', 'NKG7', 'PRF1'))

Feature_rast(GDTcell, c('ITGA1', 'ITGAE', 'KLRG1', 'KLF2', 'CSF1', 'CSF2', 'AREG','IL2RA','TNFRSF18', 'LGALS3'),sz = 1)

GDTcell@reductions
# 
# for (i in seq(0.5,0.9,0.1) %>% rev()) {
#   GDTcell <- FindClusters(GDTcell, resolution = i)
# }
ClusterCompare(GDTcell, '0', '1',min.pct = 0.1)



# GESA

CP12 <- entrezlist_generator(GDTcell, '0', '1')



library(clusterProfiler)
library(msigdbr)
Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

msigdbr::msigdbr_collections() %>% as.data.frame()

HALLMARK <-  msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

KEGG  <-   msigdbr::msigdbr(species = "Homo sapiens", category = "C2",subcategory = 'CP:KEGG') %>%
  dplyr::select(gs_name, entrez_gene)

GO<-  msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP') %>%
  dplyr::select(gs_name, entrez_gene)

ALL_msigdbr <- rbind(Mc7, HALLMARK,KEGG,GO)

GSEACP12 <- GSEA(geneList = CP12, TERM2GENE=ALL_msigdbr,
                    # minGSSize    = 10,
                    pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% 
setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")


GSEA_multipplot(GSEACP12, c("GSE25087_TREG_VS_TCONV_ADULT_UP",
                            'GOBP_CELL_KILLING'
                            ),title = 'Cluster 0 vs 1', c1 = 0, c2 = 1)
P

view(GSEACP12@result)


grep('TREG',GSEACP01@result$Description, value = T )

grep('CD8',GSEACP01@result$Description, value = T )


gseaplot2(GSEACP01, geneSetID = 'GSE25087_TREG_VS_TCONV_ADULT_UP')

GSEA_multipplot(GSEACP12, 'GSE25087_TREG_VS_TCONV_ADULT_UP',title = 'GSE25087_TREG_VS_TCONV_ADULT_UP', c1 = 0, c2 = 1)

GSEACP01_KEGG <- GSEA(geneList = KEGG, TERM2GENE=Mc7,
                 minGSSize    = 10,
                 pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")



GSEACP01_GO  <-gseGO(CP01,OrgDb = 'org.Hs.eg.db',ont = 'BP',pvalueCutoff = 0.05) %>%  simplify()%>% setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")

view(GSEACP01_GO@result)

GSEACP01_GO@result$Description %>% str_extract('killer')
GSEACP01_GO@result %>%  filter(grepl('killer', Description)) %>%  pull(ID)


GSEA_multipplot(GSEACP01_GO, 'GO:0042267',title = 'natural killer cell mediated cytotoxicity', c1 = 0, c2 = 1, base_size = 10)



Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

Hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'BP') %>%
  dplyr::select(gs_name, entrez_gene)

Mc2 <-  msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene)



Mc7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") %>%
  dplyr::select(gs_name, entrez_gene)

Hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'BP') %>%
  dplyr::select(gs_name, entrez_gene)

Mc2 <-  msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene)

Msig_refs <- list(Mc7, Mc2, Hallmarks, GO_BP) %>%  setNames(c('Mc7', 'Mc2', 'Hallmarks', 'GOBP'))
Msig_refs$GOBP
ent2_10 <- readRDS('ent2_10.rds') 


C2_C10_msigdb <- map2(Msig_refs,c('Mc7', 'Mc2', 'Hallmarks', 'GOBP'), function(x,y) {
  ge <-   GSEA(geneList = CP12, TERM2GENE=x,  
               pvalueCutoff = 0.05) %>%
    setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID") 
  ge@result %>% arrange(desc(NES))  %>%
    write.xlsx(paste0('CD4CD8c2_to_c10_Msigdb_',y,'_.xlsx') )
  print(nrow(ge@result))
  return(ge)
}) %>%  setNames(c('Mc7', 'Mc2', 'Hallmarks', 'GOBP'))

saveRDS(C2_C10_msigdb, 'CD4CD8_c2_vs_c10_GSEA.rds')


aGDTcell$donor %>% table

table(GDTcell$dataset_origin)[table(GDTcell$dataset_origin) > 1]


Feature_rast(GDTcell, 'disease', sz = 2)



GDTcell_adams <- subset(GDTcell, dataset_origin == 'adams_2020')



Feature_rast(GDTcell_adams, c('ITGA1', 'ITGAE', 'KLRG1', 'KLF2','ITGB2', 'ZNF683', 'CXCR6', 'CD69', 'GZMA', 'GZMB', 'NKG7', 'PRF1', 'AREG'))


GDTcell$Disease %>%  unique()




GDTcell@meta.data  %<>%  mutate(Disease = case_when(disease == 'chronic obstructive pulmonary disease' ~ 'COPD',
                                                 
                                                   disease == 'idiopathic pulmonary fibrosis' ~ 'IPF',
                                                   !is.na(disease) ~ 'normal'
                                                   )
                             )

Feature_rast(GDTcell, c('ident'), sz = 2, do.label = F, mythe = F)
Feature_rast(GDTcell, c('Disease'), sz = 2, do.label = F, mythe = F, colorset = 'gg')


ViolinPlot(GDTcell_adams, c('TRDC', 'TRGC1', 'TRGC2'), colors = umap.colors)





Feature_density(GDTcell,c('ITGA1', 'ITGAE', 'ZNF683' , 'KLRG1', 'KLF2','ITGB2'), ncol = 3, mythe = F, othertheme = NoLegend())


replace(c(1,2,3), 1, 3)



GDTcell$seurat_clusters

portions <- GDTcell@meta.data %>% count(seurat_clusters, disease, donor) %>%  group_by(disease, donor) %>%  mutate(per = n/sum(n)*100)

GDTcell$dataset_origin

Feature_rast(GDTcell, d1 ='TRAC', d2 = 'TRDC', assay = 'RNA')


ggplot(portions, aes(x = disease, y = per) ) +geom_jitter()+geom_boxplot() + facet_grid(~seurat_clusters)




# Figure 6 ----------------------------------------------------------------

Lung_satija$annotation.l1



GDTcell_adams@meta.data  %<>%  mutate(Disease = case_when(disease == 'chronic obstructive pulmonary disease' ~ 'COPD',
                                                    
                                                    disease == 'idiopathic pulmonary fibrosis' ~ 'IPF',
                                                    !is.na(disease) ~ 'normal'
)
)

Feature_rast(Lung_satija, c("CD3D", "TRBC1", "TRDC"))

F6A <- Feature_rast(GDTcell_adams, c("ident"),
             do.label = F,mythe = T, colorset = set_sample(umap.colors,s = 1),
             sz = 1, noaxis = F)

F6B <- Feature_rast(GDTcell_adams, c("Disease"),
             do.label = F,mythe = T, colorset = set_sample(umap.colors,s = 12),
             sz = 1, noaxis = F) %T>% print() 


F6AB <- PG(list(F6A, F6B), ncol = 1, align = "hv", axis  = "r", labels = "AUTO") %T>% print() 



GDTcell_adams$donor

numberofgdt <-   GDTcell_adams@meta.data %>% group_by(disease,donor) %>% count(seurat_clusters)


ggplot(numberofgdt , aes(x = seurat_clusters     , y = n , color =  disease )) +
  geom_jitter()+facet_wrap(~disease)


F6C <- Feature_rast(GDTcell_adams,c("ZNF683" ,'ITGA1', 'ITGAE', "KLRG1",  'ITGB2', 'KLF2') ,
                    sz = 1, ncol = 3,
                    othertheme =theme(
                      legend.margin = margin(0,0,0,-10, "pt"))) %>% list(NA) %>% PG(ncol = 1, rh = c(1, 0.1))
F6C

F6D <- GSEA_multipplot(GSEACP12, c("KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                                   'GOBP_CELL_KILLING'
),title = 'GSEA cytotoxicity', c1 = 0, c2 = 1, base_size = 8
)
  

F6D1 <- gseaplot2(GSEACP12, "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",base_size = 8, title = "NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY") %>% ggplotify::as.ggplot()
F6D2 <- gseaplot2(GSEACP12, "GOBP_CELL_KILLING",base_size = 8, title = "GOBP_CELL_KILLING")%>% ggplotify::as.ggplot()




F6D

F6E <-  GSEA_multipplot(GSEACP12, c("GSE25087_TREG_VS_TCONV_ADULT_UP",
                                    'GSE7852_TREG_VS_TCONV_FAT_UP'
),title = 'GSEA regulatory', c1 = 0, c2 = 1, base_size = 8
)


F6E1 <- gseaplot2(GSEACP12, "GSE25087_TREG_VS_TCONV_ADULT_UP",base_size = 8, title = "TREG_VS_TCONV_ADULT_UP")%>% ggplotify::as.ggplot()


F6E2 <- gseaplot2(GSEACP12, "GSE7852_TREG_VS_TCONV_FAT_UP",base_size = 8, title = "TREG_VS_TCONV_FAT_UP")%>% ggplotify::as.ggplot()


F6DE <- PG(list(F6D,NA, F6E,NA), labels = c("D",NA, "E", NA), ncol = 1, rh = c(1, 0.05, 1, 0.05))%T>% 
  print() %>% list(NA) %>% PG(ncol = 2, rw = c(1, 0.2))
F6DE <- PG(list(F6D, F6E), labels = c("D","E"), ncol = 2)%T>% 
  print() %>% list(NA) %>% PG(ncol = 1, rh = c(1, 0.05))



GSEACP12@result %>%  filter(ID == "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY")

F6DE


F6ABC <- PG(list(F6AB, F6C), labels = c(NA, "C"), ncol = 2, rw = c(1, 2)) %T>% print()

F6 <- (PG(list(F6ABC, F6DE), ncol = 1, rh = c(1, 1))+
         draw_figure_label('Figure 6',  
                           
                           position = 'top.right', size = 10, fontface = 'plain'))%T>% 
  print() %>% figsave("Fig6_publicdata.2.pdf", 160, 130)


Feature_density(GDTcell_adams,c("ZNF683" ,'ITGA1', 'ITGAE', "ENTPD1",
                                
                                "KLRG1",  'ITGB2', 'GZMB', "PRF1") ,
             sz = 1, ncol = 4) 

Feature_rast(GDTcell_adams,c("ZNF683" ,'ITGA1', 'ITGAE', "ENTPD1",
                                
                                "KLRG1",  'ITGB2', 'GZMB', "PRF1") ,
                sz = 1, ncol = 4) 
