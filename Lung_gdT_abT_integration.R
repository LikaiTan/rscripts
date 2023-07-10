
# lung project 6 patients -------------------------------------------------
# abT and gdT integration




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
# library(rstatix)
library(stringr)
library(magrittr)
library(Nebulosa)

# themes and functions ------------------------------------------------------------------
source('/home/big/tanlikai/script/rscripts/funcs.r')
file.edit('/home/big/tanlikai/script/rscripts/funcs.r')

# dir.create('/home/big/tanlikai/Lung/inte')
# dir.create('/home/big/tanlikai/Lung/inte/figs')

setwd("/home/big/tanlikai/Lung/inte")


# load projects -----------------------------------------------------------


GDTlung.rds <- '../GDTlung280622.rds'

GDTlung_s <- readRDS(GDTlung.rds)

Feature_rast(GDTlung_s)

CD4CD8 <- readRDS('../abt/CD4CD8_integrated_2022_8p.rds')

CD4CD8$TCR <- 'ABT'

GDTlung_s$TCR <- 'GDT'

Key(GDTlung_s)

CD4CD8$Cell_cluster <- paste0('C', Idents(CD4CD8))

Feature_rast(CD4CD8)



CD4CD8@assays$RNA@key <- 'rna_'

CD4CD8@assays$HTO@key <- 'hto_'

CD4CD8@assays$CITE@key <- 'cite_'

  CD4CD8@assays$integrated@key <- 'integrated_'

  CD4CD8@assays$GM@key <- 'GM_'
  
DefaultAssay(CD4CD8)
  

# split projects ----------------------------------------------------------
  
  
GDTlung_ss <-  GDTlung_s  %>% SplitObject(split.by ='orig.ident')


  CD4CD8s  <- CD4CD8%>% SplitObject(split.by ='orig.ident')

  anchors <- FindIntegrationAnchors(c(GDTlung_ss, CD4CD8s), dims = 1:60)
  
Lung_integration <- IntegrateData(anchors, dims = 1:60)

DefaultAssay(Lung_integration) <- "integrated"


Lung_integration %<>% ScaleData(vars.to.regress = c("percent.mito", 
                                             'patient',
                                             # "S.Score", 'G2M.Score',
                                             # "percent.ribo",
                                             'nCount_RNA',
                                             'nFeature_RNA' ))


Lung_integration %<>% RunPCA( npcs = 100, verbose =
                         T, nfeatures.print = 40) %>% 
  JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 

JackStrawPlot(Lung_integration, dims = 1:50 ) %T>%
  figsave('Lung_integration.jackstraw.pdf' , w = 400, h = 400)


Lung_integration  %<>%  RunUMAP( dims = 1:42, 
                          reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:42)) %>% FindClusters( resolution = 1.2)


Feature_rast(Lung_integration)

Feature_rast(Lung_integration, c('ident', 'pheno', 'tissue', 'TCR'), ncol = 2)

Feature_rast(Lung_integration, 'ITGAE', assay = 'RNA')

# 2nd strategy ------------------------------------------------------------
anchors2 <- FindIntegrationAnchors(c(GDTlung_s, CD4CD8), dims = 1:60)

multicores()

Lung_integration2 <- IntegrateData(anchors2, dims = 1:60)

DefaultAssay(Lung_integration2) <- "integrated"


Lung_integration2 %<>% ScaleData(vars.to.regress = c("percent.mito", 
                                                    'patient',
                                                    # "S.Score", 'G2M.Score',
                                                    # "percent.ribo",
                                                    'nCount_RNA',
                                                    'nFeature_RNA' ))

Lung_integration2 %<>% ScaleData(vars.to.regress = c("percent.mito", 
                                                     'patient',
                                                     # "S.Score", 'G2M.Score',
                                                     # "percent.ribo",
                                                     'nCount_RNA',
                                                     'nFeature_RNA' ), assay = 'RNA')



Lung_integration2 %<>% RunPCA( npcs = 100, verbose =
                                T, nfeatures.print = 40) %>% 
  JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 

JackStrawPlot(Lung_integration2, dims = 1:80) %T>%
  figsave('Lung_integration.jackstraw2.pdf' , w = 400, h = 400)


Lung_integration2  %<>%  RunUMAP( dims = 1:51, 
                                 reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:51)) %>% FindClusters( resolution = 1)


Feature_rast(Lung_integration2, colorset = 'gg')

Feature_rast(Lung_integration2, c('pheno', 'tissue', 'TCR','CD4CD8'), ncol = 2)

Feature_rast(Lung_integration2, 'ITGAE', assay = 'RNA')

Feature_rast(Lung_integration2, 'Cell_cluster', 'TCR', colorset = set_sample(ggplotColours(30))  )

ClusterCompare(Lung_integration2, 'P2', 'C2', group.by = 'Cell_cluster',rm = '^TRD|^TRG|^TRA|^TRB|^HSP|^RP|^MT' )

ClusterCompare(Lung_integration2, 'P2', 'C10', group.by = 'Cell_cluster', rm = '^TRD|^TRG|^TRA|^TRB|^HSP|^RP|^MT' )


ClusterCompare(Lung_integration2, 'P8', 'C7', group.by = 'Cell_cluster' )
ClusterCompare(Lung_integration2, 'M1', 'C4', group.by = 'Cell_cluster' )

rm(Lung_integration2)


# integration with human PBMC gdT ----------------------------------------------
humanGDT <- readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/GDT_2020AUG_woCOV.rds')
humanGDT@assays$GM <- 'GM_'
humanGDT$patient <- humanGDT$orig.ident
humanGDT$tissue <- 'blood'

humanGDT$Cell_cluster


humanGDT@meta.data  %<>%  mutate(pheno = case_when(
  Cell_cluster %in% c('c1', 'c2', 'c3', 'c4') ~ 'Naive/CM',
  Cell_cluster == 'c5' ~ 'type2-Vg9Vd2',
  Cell_cluster == 'c6' ~ 'type3-Vg9Vd2',
  Cell_cluster %in% c('c7', 'c8', 'c9') ~ 'type1-Vg9Vd2',
  Cell_cluster %in% c('c10', 'c11') ~ 'nonVd2-CTL',
  Cell_cluster == 'c12' ~ 'activated nonVd2'

  
  
))


Feature_rast(humanGDT, 'pheno', facets = 'group')


DefaultAssay(GDTlung_s)

humanGDTs <-  humanGDT  %>% SplitObject(split.by ='orig.ident')

humanGDTs$CB_donor1 <- NULL

humanGDTs$CB_donor2 <- NULL

GDTlung_ss%<>% map(.,~ NormalizeData(.x, normalization.method = 'LogNormalize',
                           scale.factor = 10000, assay = 'RNA') %>% 
           FindVariableFeatures(assay = 'RNA',nfeatures = 2500, selection.method = 'vst')
)


humanGDTs%<>% map(.,~ NormalizeData(.x, normalization.method = 'LogNormalize',
                                     scale.factor = 10000, assay = 'RNA') %>% 
                   FindVariableFeatures(assay = 'RNA',nfeatures = 2500, selection.method = 'vst')
)

for (i in names(GDTlung_ss)) {
  GDTlung_ss[[i]]@assays$RNA@var.features%<>%  
    str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T)
}

for (i in names(humanGDTs)) {
  humanGDTs[[i]]@assays$RNA@var.features%<>%  
    str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T)
}

# Feature_rast(humanGDT, 'Cell_cluster')
# Feature_rast(GDTlung_s, 'Cell_cluster')



anchors <- FindIntegrationAnchors(c(GDTlung_ss, humanGDTs), dims = 1:60)
GDT_integration<- IntegrateData(anchors, dims = 1:60)

DefaultAssay(GDT_integration) <- "integrated"



GDT_integration@assays$integrated@var.features %<>%  
  str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB|^HIST|^RP|^HSP', negate = T)





GDT_integration %<>% ScaleData(vars.to.regress = c(
                                                   'nFeature_RNA' )) %>%
  RunPCA( npcs = 100, verbose =T, nfeatures.print = 40)
ElbowPlot(GDT_integration, ndims = 100)
GDT_integration %<>% ScaleData(vars.to.regress = c(
  'nFeature_RNA' ), assay = 'RNA' ) 

GDT_integration %<>% 
  JackStraw( num.replicate = 100, dims = 80)%>%
  ScoreJackStraw(dims = 1:80) 
JackStrawPlot(GDT_integration, dims = 1:80)


 GDT_integration   %<>% RunUMAP( dims = 1:43, 
                                reduction = 'pca', min.dist = 0.2) %>%
  FindNeighbors(dims = c(1:43)) %>% FindClusters( resolution = 0.6)
 

Feature_rast(GDT_integration, c('integrated_snn_res.0.6', 'Cell_cluster', 'tissue'), colorset = 'gg' )
 ClusterCompare(GDT_integration, '8', '9',  assay = 'RNA')

GDT_integration$bc_backup <- colnames(GDT_integration)
GDT_integration$integrated_snn_res.0.6 %>%  unique()

Feature_rast(GDT_integration, c('ident', 'Cell_cluster'), facets = 'tissue',ncol = 1  ,colorset = 'gg' , sz = 0.2)





GDT_integration@meta.data %<>%  mutate(phenotype = case_when(integrated_snn_res.0.6 %in% c('10', '11') ~ 'type3_Vg9Vd2',
                                                            integrated_snn_res.0.6 %in% c('0','1') ~ 'type1_Vg9Vg2',
                                                            integrated_snn_res.0.6 == '8' ~ 'Naive/Tcm',
                                                            integrated_snn_res.0.6 %in% c('2', '4', '13') ~ 'cytotoxic Tem_nonVd2',
                                                            integrated_snn_res.0.6  %in%c('6', '7') ~ 'cytotoxic Teff_nonVd2',
                                                            integrated_snn_res.0.6 ==3 ~ 'lung Trm1',
                                                            integrated_snn_res.0.6 == 9 ~ 'lung Trm2',
                                                            integrated_snn_res.0.6== 5 ~ 'lung Tcric' )) %>% `rownames<-`(GDT_integration$bc_backup)

GDT_integration$GM_D
  

Feature_rast(GDT_integration, 'phenotype')


Feature_rast(GDT_integration, c('phenotype', 'Cell_cluster', 'tissue'), colorset = 'gg' )

saveRDS(GDT_integration, 'GDT_integration_BL_LNG.RDS')

GDT_integration <- readRDS('GDT_integration_BL_LNG.RDS')


Feature_rast(readRDS('GDT_integration_BL_LNG.RDS'), c('phenotype', 'Cell_cluster', 'integrated_snn_res.0.6'), colorset = 'gg')






# Label transfering with human blood data ---------------------------------



DefaultAssay(humanGDT) <- 'RNA'
DefaultAssay(GDTlung_s) <- 'RNA'



GDT.anchors <- FindTransferAnchors(reference = subset(humanGDT, group == 'Adult'),
                                   query = GDTlung_s,query.assay = 'integrated',
                                   reference.assay = 'integrated',
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = GDT.anchors, 
                            refdata = subset(humanGDT, group == 'Adult')$pheno,
                            dims = 1:30)
GDTlung_s <- AddMetaData(GDTlung_s, metadata = predictions)
Feature_rast(GDTlung_s, 'predicted.id')

predictions
DefaultAssay(GDTlung_s) <- 'RNA'


# # Harmony -----------------------------------------------------------------
# library(harmony)
# HGDT_integration<- c(GDTlung_ss, humanGDTs) %>% reduce(.f = merge, project = 'lungGDT')
# 
#
# 
# 
# 
# 
# HGDT_integration%<>% NormalizeData(verbose = FALSE, assay = 'RNA') %>% 
#   FindVariableFeatures(selection.method = 'vst', nfeatures = 3500) 
# 
# 
# HGDT_integration@assays$RNA@var.features %<>%str_subset('^TRG|^TRD|^MT|^IG|^TRA|^TRB^|^HIST|^MIR|^HSP', negate = T) #%T>% print(length())
# 
# #%>% 
# HGDT_integration%<>%   ScaleData(verbose = FALSE, assay = 'RNA',    
#                         vars.to.regress = c('percent.mito',  'nFeature_RNA' )
# ) %>% 
#   RunPCA(npcs = 100, verbose=T)
# HGDT_integration%<>% RunHarmony(c("orig.ident"), plot_convergence = TRUE)
# 
# ElbowPlot(HGDT_integration, ndims = 50,reduction = 'harmony')
# 
# saveRDS(HGDT_integration, 'HGDT_integration.rds')
# 
# Feature_rast(HGDT_integration, c('phenotype', 'Cell_cluster', 'tissue'), colorset = 'gg' )
# 


# Gene modules ------------------------------------------------------------
library(openxlsx)


gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_c_list_ent$GM
gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)
gene_cluster$GM_B

GDT_integration@meta.data  %<>% select(-starts_with('GM_'))


for (i in gc_name) {
  GDT_integration %<>%  AddModuleScore( features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(GDT_integration@meta.data) %<>% str_replace("(?<=\\w)1", '')
  
}
colnames(GDT_integration@meta.data)

Feature_rast(GDT_integration, c(gc_name), color_grd = 'grd') 

Feature_rast(GDT_integration, c('GM_D'), color_grd = 'grd', facets = 'tissue') 



ViolinPlot(GDT_integration, gc_name, colors = umap.colors)

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
                      "GM_F: CTL response (innate)",
                      "GM_G: CTL response (adaptive)",
                      
                      "GM_H: Acute activation"))


c(1,2,4,6,7) %>% map(~
                       Feature_rast(GDT_integration, gc_name[.x], color_grd = 'grd', mythe = F)+
                       ggtitle(gcanno[.x])
) %>% PG(ncol = 3) %T>% figsave('GM_score.pdf', 200, 150)

# 
# GM_score_heat <- DoHeatmap(subset(GDTlung_s, downsample = 1000), 
#                            raster =T, draw.lines = T, angle = 45,
#                            lines.width = 10,group.colors = umap.colors,
#                            
#                            assay = 'GM', features = gms, slot = 'data', size = gs(8)) +hmp2 + mytheme+
#   theme(legend.position = 'bottom',
#         legend.key.height = unit(2,'mm'))+
#   guides(color = FALSE, fill = guide_colourbar(title = 'Scaled modula score', title.position = 'top'))
# GM_score_heat


sigtable <- read.xlsx('/home/big/tanlikai/Lung/abt/abd5778_Table_S3.xlsx',sheet = 1)

colnames(sigtable) %<>%   str_remove('.\\(.+\\)')


sigtable %<>% as.list() %>% map(~ na.exclude(.x) %>% as.vector)

Feature_rast(GDTlung_s, sigtable$Tissue.resident)

sigtable$Tissue.resident

for (i in names(sigtable)) {
  GDT_integration <- AddModuleScore(GDT_integration, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}
colnames(GDT_integration@meta.data) %<>% str_replace("(?<=\\w)1$", '')


GDT_integration %<>%  AddModuleScore( features = list(sigtable$Tissue.resident[c(1:4,7,8)]), name = 'Tissue.resident', assay = 'RNA')
GDT_integration$Tissue.resident <- GDT_integration$Tissue.resident1



Feature_rast(GDT_integration,c(names(sigtable)[c(1:9, 14,15)]), color_grd = 'grd')


Feature_rast(GDT_integration, 'Tissue.resident')-Feature_rast(GDT_integration, 'GM_A')
Feature_rast(GDT_integration, 'tissue')

GDT_integration@meta.data  %<>%  mutate(Cell_cluster2 =  case_when( 
  tissue != 'blood' ~ Cell_cluster,tissue == 'blood' ~ 'blood gdT'
                                                                   ))

GDT_integration$Cell_cluster2  %<>%  factor(levels = c('blood gdT', paste0('P', 1:9), paste0('L', 1:4), 'M1'))


Feature_rast(GDT_integration, 'Cell_cluster2', colorset = c('black', umap.colors),sz = 0.2, do.label = F)



Feature_rast(GDTlung_s)

