# CIPD data 


# load packages -----------------------------------------------------------



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
file.edit('/home/big/tanlikai/script/rscripts/funcs.r')


setwd("/home/big/tanlikai/Lung/")


COPD_pub <- readRDS("public/pdd_copd_public_2023-09-18_revised.RDS")





saveRDS(COPD_pub, "public/pdd_copd_public_2023-09-18_revised.RDS")

# COPD_pub <- UpdateSeuratObject(COPD_pub)

# colnames(COPD_pub) <- COPD_pub@reductions$umap@cell.embeddings %>%  rownames()

rownames(COPD_pub@meta.data) <- COPD_pub@reductions$umap@cell.embeddings %>%  rownames()




figpath_ni <-  "/home/big/googledrive/Lungfigures/"

DefaultAssay(COPD_pub) <-  "RNA"

Idents(COPD_pub)

COPD_pub$Level_4


Feature_rast(COPD_pub, colorset = 'gg', g = c("TRDC", "v_gene", "TRAC", "CD3E",  "AREG", "ITGAE", "ITGA1", "GATA3"), 
             do.label = F, othertheme = NoLegend())


COPD_pub@meta.data$doublet_scrublet


Feature_rast(COPD_pub, c("qc.nCount.Outlier", "qc.mit.Outlier", "doublet_scrublet"))




ViolinPlot(COPD_pub, "doublet_scrublet", group.by = "orig.ident")

ViolinPlot(COPD_pub, "IL32", group.by = "disease")



# # add module score  -----------------------------------------------------



CD3list <- list(c("CD3D", "CD3E", "CD3G"))
CD3list

COPD_pub  %<>% AddModuleScore(features = CD3list, name = "CD3_score")


COPD_pub$CD3_score1



TRDlist <- list(c("TRDV1", "TRDV2", "TRDV3", 'TRDC'))

TRABlist <- list(
  grep("^TRAV|^TRBV|TRBC1|TRBC2|^TRAC", rownames(COPD_pub), value = T)
  
)

TRABlist



COPD_pub  %<>% AddModuleScore(features = TRDlist, name = "TRD_score")
COPD_pub  %<>% AddModuleScore(features = TRABlist, name = "TRAB_score")



Feature_rast(COPD_pub, c("CD3_score1", "TRD_score1", "TRAB_score1") )



# select T cells  ---------------------------------------------------------

COPD_pub@meta.data %>% filter(qc.nFeat.Outlier == "TRUE" | qc.mit.Outlier == "TRUE")


COPD_pub  %<>% subset(qc.nFeat.Outlier != "TRUE" ) 

ViolinPlot(COPD_pub, group.by =   "Level_1", g= "CD3_score1", x.angle = 90)


ViolinPlot(COPD_pub, group.by =   "Level_3", g= "TRD_score1", x.angle = 90)


ViolinPlot(COPD_pub, group.by =   "Level_3", g= "TRAB_score1", x.angle = 90)


COPD_Tcells <- subset(COPD_pub, CD3_score1 > 0.3)



#

# select gdT cells  -------------------------------------------------------


Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", 
             noaxis = F, axis.number = T, do.label = F, g = "v_gene", colorset = "gg", navalue = "blue")+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0.018, linewidth = 0.2)


COPD_Tcells@meta.data   %<>% mutate(gdTcells = case_when( 
  TRAB_score1 <= 0.018 & TRD_score1 >=  0.23 & is.na(v_gene) ~ TRUE, TRUE ~ FALSE
))


Feature_rast(COPD_Tcells, c("gdTcells"), colorset = "gg", do.label = F)



Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", 
             noaxis = F, axis.number = T, do.label = F, g = "gdTcells", 
             colorset = c("lightgrey", "red"))+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0.018, linewidth = 0.2)





table(COPD_Tcells$gdTcells, COPD_Tcells$Level_4)




COPD_gdTcells <- subset(COPD_Tcells,   gdTcells == TRUE)



dim(COPD_gdTcells)

table(COPD_gdTcells$disease
)

Feature_rast(COPD_Tcells, c("gdTcells", "TCRgamma-or-delta-prot", "TCRVdelta2-prot" ,),  assay = "Protein",
             colorset = c(
  alpha("lightgrey", 0.3), "red"
), do.label = F)


rownames(
  COPD_Tcells@assays$Protein@meta.features)

PG(list(
ViolinPlot(COPD_Tcells,  g = "CD3-prot", group.by = "Level_4" , assay = "Protein",x.angle = 90),


ViolinPlot(COPD_Tcells,  g = "CD3_score1", group.by = "Level_4" , assay = "Protein",x.angle = 90) ),
ncol = 1)



PG(list(
  ViolinPlot(COPD_gdTcells,  g = "CD3-prot", group.by = "Level_4" , assay = "Protein",x.angle = 90),
  
  
  ViolinPlot(COPD_gdTcells,  g = "CD3_score1", group.by = "Level_4" , assay = "Protein",x.angle = 90) ),
  ncol = 1)



PG(list(
  ViolinPlot(COPD_gdTcells,  g = "CD3-prot", group.by = "Level_4" , assay = "Protein",x.angle = 90),
  
  
  ViolinPlot(COPD_gdTcells,  g = "CD3_score1", group.by = "Level_4" , assay = "Protein",x.angle = 90) ),
  ncol = 1)



ViolinPlot(COPD_gdTcells, c("TRD_score1", "TCRgamma-or-delta-prot", "TCRVdelta2-prot" , "TCR-alpha-or-beta-prot") , assay = "Protein", 
           x.angle = 90,group.by = "Level_4" , colors = ggplotColours(35), ncol = 1
           )






# reclustering  -----------------------------------------------------------
COPD_gdTcells$disease %>% table()
COPD_gdTcells$orig.ident %>% table()

table(COPD_gdTcells$orig.ident,COPD_gdTcells$disease)

COPD_gdTcells  %<>% subset(disease != "End-stage COPD")


COPD_gdTcells<- NormalizeData(COPD_gdTcells, normalization.method = 'LogNormalize',
                        scale.factor = 10000, assay = 'RNA')



COPD_gdTcells  %<>%    FindVariableFeatures(assay = 'RNA',nfeatures = 6000, selection.method = 'vst')
%>%  ScaleData( assay = 'RNA',
                           vars.to.regress = c('nCount_RNA', 'percent.mito')) 
COPD_gdTcells@assays$RNA@var.features
COPD_gdTcells@assays$RNA@var.features <- GDTcell@assays$RNA@var.features%>%  
  str_subset('^RP|^MT|TRAV|TRBV|NO-NAME|IGL|IGH|LINC|MIR', negate = T)


COPD_gdTcells <- RunPCA(COPD_gdTcells, npcs = 100, verbose =
                    T, nfeatures.print = 40)



ElbowPlot(COPD_gdTcells, ndims = 100)

COPD_gdTcells  <- JackStraw(COPD_gdTcells, num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 

JackStrawPlot(COPD_gdTcells, dims = 1:50 )


COPD_gdTcells <- RunUMAP(COPD_gdTcells, dims = 1:13, reduction.name = "UMAP",reduction.key = "UMAP_",
                   reduction = 'pca') %>%
  FindNeighbors(dims = c(1:13))

COPD_gdTcells <- FindClusters(COPD_gdTcells, resolution = 0.7)
Feature_rast(COPD_gdTcells, c('ident', "disease", "AREG"))

Feature_rast(COPD_gdTcells, c('ident', "disease", "ITGAE", "KLRG1", "AREG"))
COPD_gdTcells$orig.ident %>% table()

Feature_density(COPD_gdTcells, c("ITGAE", "KLRG1", "AREG", "TRDV1", "TRDV2", 
                                 "TRGV9", "CCR6", "DPP4", "ITGB2", "ITGA1", "IFNG", "CSF1"))

COPD_gdTcells@reductions$umap

table(COPD_gdTcells$disease, COPD_gdTcells$gd_cluster) %>%  prop.table(margin = 1)*100

ViolinPlot(COPD_gdTcells, "AREG", group.by = "disease")


Feature_rast(COPD_gdTcells, "orig.ident",  facets = 'disease', do.label = F)

Feature_rast(COPD_pub, "orig.ident",  facets = 'disease', do.label = F, sz = 0.1, colorset =  ggplotColours(n = 29) )



Feature_rast(COPD_gdTcells, c("ITGAE", "KLRG1", "AREG", "TRDV1", "TRDV2", 
                                 "TRGV9", "CCR6", "DPP4", "ITGB2", "ITGA1", "IFNG", "CSF1"))
COPD_gdTcells$orig.ident


Feature_rast(COPD_gdTcells, c("ITGAE", "KLRG1", "AREG", "TRDV1", "TRDV2", "GZMA", "GZMB", "KLRD1",  "CD8A", "CD8B",
                              "TRGV9", "CCR6", "DPP4", "ITGB2", "ITGA1", "IFNG", "CSF1", "ZNF683", "KLF2", 'GATA3'))

# rename  -----------------------------------------------------------------





COPD_gdTcells@meta.data  %<>% mutate(gd_cluster = case_when(seurat_clusters == "0" ~ "gdTRM_1",
                                                            seurat_clusters ==  "1" ~ 'gdTemra_1',
                                                            seurat_clusters == "2" ~ "gdTRM_2",
                                                            seurat_clusters == "3" ~ "gdTRM_3",
                                                            seurat_clusters == "4" ~"Type1_Vd2",
                                                            seurat_clusters == "5" ~ "Type3_Vd2",
                                                            seurat_clusters == "6"~ "gdTemra_2",
                                                            seurat_clusters == "7" ~ "unidentified"))

COPD_gdTcells  %<>% subset(gd_cluster != "unidentified")



COPD_gdTcells$gd_cluster %>% unique() %>% sort()



COPD_gdTcells@meta.data  %<>% mutate(gd_cluster = factor(gd_cluster, 
                                                         level = 
                                                           gd_cluster %>% unique() %>% sort()))

Idents(COPD_gdTcells) <-  'gd_cluster'


F7D <-  Feature_rast(COPD_gdTcells, c("ITGAE", "ITGA1", "ZNF683", "KLRG1", "ITGB2",
                                      "KLF2", 'TRDV1' , 'TRDV2', "TRDV3", "TRGV9"),  othertheme = coord_fixed()) %T>% print()



ViolinPlot(COPD_gdTcells,c( "AREG", "CSF1", "IL22"), group.by = "disease", colors = umap.colors)

proteinlist <-  COPD_gdTcells@assays$Protein@meta.features %>% rownames() %>% sort()


proteinlist

COPD_gdTcells$orig.ident %>% table()
Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", assay = "Protein",
             noaxis = F, axis.number = T, do.label = F, g = "TCRgamma-or-delta-prot", colorset = c("lightgrey", "red"))



saveRDS(COPD_gdTcells, "public/COPD_gdTcells.rds")


table(COPD_gdTcells$gd_cluster, COPD_gdTcells$disease) %>%  prop.table(margin = 2)*100



Feature_rast(COPD_Tcells, "Level_4", colorset = "gg")
Feature_rast(COPD_gdTcells, "Level_4", colorset = "gg")



# abundance ---------------------------------------------------------------



copdratio <-  (table(COPD_gdTcells$disease, COPD_gdTcells$gd_cluster) %>%  
                 prop.table(margin = 1)*100) %>% data.frame() %>% filter(!is.na(Freq)) %>% 
  `colnames<-`(c("disease", "gd_cluster","Percent"))



copdratio_donor <-  COPD_gdTcells@meta.data %>%  
  group_by(orig.ident, disease) %>% count(gd_cluster) %>%  mutate(percent = n/sum(n)*100, SUM = sum(n))

copdratio_donor$SUM %>%  unique()





ggplot(copdratio_donor %>%  filter(SUM > 20), aes(x = disease, y = percent)) + geom_point()+facet_wrap(~gd_cluster )+ theme(axis.text.x = element_text(angle = 90))






copdratio_donor_T <-  COPD_Tcells@meta.data %>%  
  group_by(orig.ident, disease) %>% count(Level_4) %>%  mutate(percent = n/sum(n)*100, SUM = sum(n))






ggplot(copdratio_donor_T %>%  filter(SUM > 50), aes(x = disease, y = percent)) + geom_point()+facet_wrap(~Level_4 )+ theme(axis.text.x = element_text(angle = 90))




library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

COPD_gdTcells@meta.data$disease %>% unique
COPD_gdTcells@meta.data  %<>% mutate(COPD_vs_Ctrl = case_when(disease %in% c("Emphysema","moking-related ILD",
                                                                             "GOLD I/II") ~ "COPD and Emphysema",  disease %in% c("Healthy", "Donor") ~ "Control" ))

Feature_rast(COPD_gdTcells, "COPD_vs_Ctrl")




copdratio_donor <-  COPD_gdTcells@meta.data %>%  filter(!is.na(COPD_vs_Ctrl)) %>% 
  group_by(orig.ident, COPD_vs_Ctrl) %>% dplyr::count(gd_cluster) %>%  mutate(percent = n/sum(n)*100, SUM = sum(n))

copdratio_donor$SUM %>%  unique()

ggplot(copdratio_donor %>%  filter(SUM > 20), aes(x = COPD_vs_Ctrl, y = percent)) + geom_point()+facet_wrap(~gd_cluster )+ theme(axis.text.x = element_text(angle = 90))


COPD_gdTcells_ab <- subset(COPD_gdTcells, COPD_vs_Ctrl %in% c('COPD and Emphysema', "Control") )

Feature_rast(COPD_gdTcells, "COPD_vs_Ctrl", colorset =c( "cyan", "blue"))

COPD_gdTcells_ab <- as.SingleCellExperiment(COPD_gdTcells_ab)
COPD_gdTcells_ab_milo <- Milo(COPD_gdTcells_ab)

data.frame(colData(COPD_gdTcells_ab_milo))


COPD_gdTcells_ab_milo  %<>% buildGraph(k = 10, d = 30) %>% 
  makeNhoods( prop = 0.1, k = 10, d=30, refined = TRUE) %>% 
  countCells(meta.data = data.frame(colData(COPD_gdTcells_ab_milo)), samples="orig.ident")


COPD_design <- data.frame(colData(COPD_gdTcells_ab_milo))[,c("orig.ident", "COPD_vs_Ctrl")] %>% 
  distinct() 

COPD_design
rownames(COPD_design) <- COPD_design$orig.ident
## Reorder rownames to match columns of nhoodCounts(milo)
COPD_design <- COPD_design[colnames(nhoodCounts(COPD_design)), , drop=FALSE]

COPD_design

COPD_gdTcells_ab_milo <- calcNhoodDistance(COPD_gdTcells_ab_milo, d=30)

da_results <- testNhoods(COPD_gdTcells_ab_milo, design = ~ COPD_vs_Ctrl, design.df = COPD_design)
da_results
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)


da_results <- annotateNhoods(COPD_gdTcells_ab_milo, da_results, coldata_col = "gd_cluster")
head(da_results)
da_results$gd_cluster <- ifelse(da_results$gd_cluster_fraction < 0.5, "Mixed", da_results$gd_cluster)
plotDAbeeswarm(da_results, group.by = "gd_cluster",alpha = 0.5)+mytheme


ClusterCompare(COPD_gdTcells, "gdTRM_1", "gdTRM_2")

# gene module -------------------------------------------------------------


sigtable <- openxlsx::read.xlsx('abt/abd5778_Table_S3.xlsx',sheet = 1)

colnames(sigtable) %<>%   str_remove('.\\(.+\\)')

sigtable %<>% as.list() %>% map(~ na.exclude(.x) %>% as.vector)

names(sigtable)

for (i in names(sigtable)) {
  COPD_gdTcells <- AddModuleScore(COPD_gdTcells, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}
colnames(COPD_gdTcells@meta.data) %<>% str_replace("(?<=\\w)1$", '')
COPD_gdTcells@meta.data

Feature_density(COPD_gdTcells,c(names(sigtable)[c(1:9, 14,15)]))

ViolinPlot(COPD_gdTcells,c(names(sigtable)[c(1:9, 14,15)]), colors = umap.colors, box = T, group.by = "gd_cluster")



gene_c_list_ent <-readRDS('/home/big/tanlikai/Human_GDT_2019/Integrated/Genemodule_list_ent_2020AUG.RDS')

gc_name <- gene_c_list_ent$GM %>% unique() %>% sort() %>% as.vector()

gene_c_list_ent$GM
gene_cluster <- map(gc_name, function(x) gene_c_list_ent %>% filter(GM == x) %>% dplyr::select(gene)  %>% pull()) %>% setNames(gc_name)
gene_cluster$GM_B

for (i in gc_name) {
  COPD_gdTcells %<>%  AddModuleScore( features = list(gene_cluster[[i]]), name = i, assay = 'RNA')
  colnames(COPD_gdTcells@meta.data) %<>% str_replace("(?<=\\w)1", '')
  
}

gcanno <- as.vector(c(' Naive or immature T cell',
                      "Innate T cell differentiation",
                      'Proliferating',
                      "Type-3 immunity",
                      "Interferon induced",
                      "CTL response (innate)",
                      "CTL response (adaptive)",
                      
                      "Acute activation",
                      'Treg core module'))


ViolinPlot(COPD_gdTcells, c("GM_C", "GM_D", "GM_F", "GM_G", "GM_H"),
           ncol = 2, box =T,
           colors = umap.colors)



# DEGs --------------------------------------------------------------------





 AllDEGs_COPDpub <-  FindAllMarkers(COPD_gdTcells,logfc.threshold = 0.25,only.pos = T,min.pct = 0.1) 
 AllDEGs_COPDpub  %<>% mutate(pct.diff = pct.1 - pct.2)
 
 
AllDEGs_COPDpub %>%  filter(gene %in% c("CCR6", "DPP4", "RORC","FCGR3A", "NKG7", "EOMES", "TNFRSF18", "CSF1", "GATA3", "AREG" ))

 scgenes <-  unique(COPD_gdTcells@assays$RNA@var.features, unique(AllDEGs_COPDpub$gene))
 
 COPD_gdTcells  %<>% ScaleData( assay = 'RNA', features = scgenes,
                                vars.to.regress = c('nCount_RNA', 'percent.mito'))
 
  view(AllDEGs_COPDpub)
 ClusterCompare(COPD_gdTcells, "Type1_Vd2", "Type3_Vd2")
 
 
 ClusterCompare(COPD_gdTcells, "gdTRM_1", "gdTemra_1")
 
 ClusterCompare(COPD_gdTcells, "gdTRM_1", "gdTRM_2")
 
 Feature_rast(COPD_gdTcells, c("ENTPD1", "PDCD1", "CTLA4"))
 
 
 AllDEGs_COPDpub <-  FindAllMarkers(COPD_gdTcells,logfc.threshold = 0.25,only.pos = T,min.pct = 0.1) 
 
 
top15DEG <- AllDEGs_COPDpub %>% group_by(cluster) %>% top_n(15, avg_log2FC)
 
View(AllDEGs_COPDpub)

top15DEG
 
 DoHeatmap(COPD_gdTcells, top15DEG$gene) %>% heat_theme()
 
 AllDEGs_COPD_TRM_TEMRA_TYPE3 %>% filter(gene == "AREG")
 
 AllDEGs_COPD_TRM_TEMRA_TYPE3 <-  FindAllMarkers(COPD_gdTcells %>% 
                                                 subset(gd_cluster %in%c("gdTRM_1", "gdTemra_1", "Type3_Vd2")),logfc.threshold = 0.25,only.pos = T,min.pct = 0.1)  %>% arrange(cluster, desc(avg_log2FC))
 
 AllDEGs_COPD_TRM_TEMRA_TYPE3
 AllDEGs_COPD_TRM_TEMRA_TYPE3  %>%  filter(cluster == "gdTRM_1")
 
 
Top20_COPD_TRM_TEMRA_TYPE3 <-AllDEGs_COPD_TRM_TEMRA_TYPE3 %>% group_by(cluster) %>% top_n(15, avg_log2FC)


DoHeatmap(COPD_gdTcells, Top20_COPD_TRM_TEMRA_TYPE3$gene) %>% heat_theme()

Dotfeatures = unique(Top20_COPD_TRM_TEMRA_TYPE3$gene) %>% 
  str_replace("MMP25", "AREG") %>% 
  str_replace("MLC1", "EOMES")%>% 
  str_replace("TMEM200A", "GATA3")

DotPlot(COPD_gdTcells,assay = 'RNA',dot.scale = 3.5,
                              features =Dotfeatures)+mytheme+
  heattheme+
  
  theme(text = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.line.y.left = element_line(),
        axis.text.x.bottom = element_text(size = 8, angle = 90),
        legend.box.margin = margin(5,0,0,15,unit = 'mm'),
        legend.box = "horizontal",legend.position = 'bottom',
        axis.title = element_blank())+
  scale_y_discrete(position = 'left')+
  scale_x_discrete(position = 'bottom')+
  xlab(NULL)+ylab(NULL)+
  scale_color_gradient2(low = '#003399', mid = '#ffccff',  high = "#990000")+
  guides(
    color = guide_colorbar(title.position = 'top',direction = 'horizontal',
    ),
    size = guide_legend(title.position = 'top',direction = 'horizontal',label.position = 'bottom'))


TFs <- readLines('HumanTFs_1639.txt') %>% as.vector()
TFs_gd <- intersect(TFs, rownames(COPD_gdTcells))
 
 
 
ClusterCompare(COPD_gdTcells, 'gdTRM_1', 'gdTemra_1', features = TFs_gd)



# protein -----------------------------------------------------------------

citeproteins <-  COPD_gdTcells@assays$Protein@counts %>%  rownames()

grep("CD39", citeproteins, value = T)


Feature_rast(COPD_gdTcells, c("CD357-or-GITR-prot"  ) , assay = "Protein" ,colorgrd = "grd2", ncol = 3, 
             othertheme = coord_fixed())


ClusterCompare(COPD_gdTcells, "gdTRM_1", "gdTemra_1", assay = "Protein")
# GSEA --------------------------------------------------------------------


genelist_c0_c1 <- Genelist_generator(COPD_gdTcells, "gdTRM_1", "gdTemra_1")

genelist_Type1_3_Vd2_c1 <- Genelist_generator(COPD_gdTcells, "Type1_Vd2", "Type3_Vd2")

genelist_TRM_1_TRM_2 <- Genelist_generator(COPD_gdTcells, "gdTRM_1", "gdTRM_2")



ibrary(clusterProfiler)
library(msigdbr)

ALL_msigdb_G  <- rbind(
  # c7
  msigdbr::msigdbr(species = "Homo sapiens", category = "C7"), 
  # hallmarker
  msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  # C2
  msigdbr::msigdbr(species = "Homo sapiens", category = "C2",subcategory = 'CP:KEGG'),
  # GOBP
  msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')) %>% 
  dplyr::select(gs_name, gene_symbol)


GOBP <-  msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP') %>% 
  dplyr::select(gs_name, gene_symbol)

GSEA_c0_c1 <- GSEA(geneList = genelist_c0_c1, TERM2GENE=ALL_msigdb_G,
                 # minGSSize    = 10,
                 pvalueCutoff = 0.05, pAdjustMethod = "BH") 

GSEA_c0_c1@result  <- GSEA_c0_c1@result %>% arrange(desc(NES))  %T>% view()



GSEA_c0_c1_GO <- GSEA(geneList = genelist_c0_c1, TERM2GENE=GOBP,
                   # minGSSize    = 10,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH") 


view(GSEA_c0_c1_GO@result)

GSEA_multipplot(GSEA_c0_c1, c("GSE7852_TREG_VS_TCONV_UP",
                            'GOBP_BRANCHING_MORPHOGENESIS_OF_AN_EPITHELIAL_TUBE'
),title = 'GSEA tissue regulatory', 
c1 = "gdTRM_1", c2 =  "gdTemra_1", base_size = 8)




 GSEA_multipplot(GSEA_c0_c1, c("GOBP_INNATE_IMMUNE_RESPONSE",
                                   'GOBP_CELL_KILLING' ),title = 'GSEA cytotoxcity', 
                 c1 = "gdTRM_c0", c2 =  "gdTemra_c1", base_size = 8 )

 
 
 GSEA_Vd2 <- GSEA(geneList = genelist_Type1_3_Vd2_c1 , TERM2GENE=ALL_msigdb_G,
                    # minGSSize    = 10,
                    pvalueCutoff = 0.05, pAdjustMethod = "BH") 
 
 
 view(GSEA_Vd2@result)
 


 Feature_rast(COPD_gdTcells, c("ITGAE", "ITGA1", "ZNF683", "KLRG1", "ITGB2",
                               "KLF2", 'TRDV1' , 'TRDV2', "TRDV3", "TRGV9"))
 
 
 Feature_rast(COPD_gdTcells, c("CD45RA-prot", "CD27-prot", "CD103-or-IntegrinalphaE-prot"  , "KLRG1-or-MAFA-prot", "CD26-prot" ,
                               "TCRVdelta2-prot"  ) , assay = "Protein" ,colorgrd = "grd2", ncol = 3, 
              othertheme = coord_fixed())
 
 
 F7A <- 
   Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", 
                noaxis = F, axis.number = T, do.label = F, g = "gdTcells", colorset = c("lightgrey", "red"))+
   geom_hline(yintercept = 0.23, linewidth = 0.2)+
   geom_vline(xintercept = 0.018, linewidth = 0.2)
 
 
 Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", 
              noaxis = F, axis.number = T, do.label = F, g = "gdTcells", colorset = c("lightgrey", "red"))
 
 
 
 F7B <-  Feature_rast(COPD_gdTcells, c("CD45RA-prot", "CD27-prot", "CD103-or-IntegrinalphaE-prot"  , "KLRG1-or-MAFA-prot", "CD26-prot" ,
                                       "TCRVdelta2-prot"  ,"CD16-prot" ) , assay = "Protein" ,colorgrd = "grd2", ncol = 3, 
                      othertheme = coord_fixed()) %T>% print()

 
 F7C <-  Feature_rast(COPD_gdTcells, c( "gd_cluster", "disease" ) , assay = "Protein" , ncol = 2, do.label = F,
                      othertheme = coord_fixed()) %T>% print()
 
 COPD_gdTcells$seurat_clusters
 
 
 
 
 
 
 
 
# Fig4 --------------------------------------------------------------------

F6A1 <- Feature_rast(COPD_pub, colorset = "gg", noaxis = T, do.label = F)+NoLegend()+ggtitle("COPD Cell Atlas")

 COPD_pub$Level_1
 
 COPD_pub$all <-  "All cells"
 
F6A2 <-  (ViolinPlot(COPD_pub, group.by =   "all",sz = 0,
                     g= "CD3_score1", x.angle = 90) + 
            geom_hline(yintercept = 0.3, color = "red") +
            theme_classic()+mytheme+
            ggtitle("T Cell Selection")+NoLegend() )%T>%  print()


F6A3  <-  (Feature_rast(COPD_Tcells, d1 = "TRAB_score1", d2 = "TRD_score1", 
                        noaxis = F, axis.number = T, do.label = F, g = "gdTcells", 
                        colorset = c("lightgrey", "red"))+
             geom_hline(yintercept = 0.23, linewidth = 0.2)+
             geom_vline(xintercept = 0.018, linewidth = 0.2) + NoLegend() +
             ggtitle("gdT Cell Selection")) %T>%  print()
  
 

F6A <-  PG(list(F6A1, F6A2, F6A3), 
           ncol = 3, rw = c(0.8,0.5,0.8), align = "hv",axis = 'tb',
           labels = c("A", NA,NA))  %T>%  print()
 



F6B <-  Feature_rast(COPD_gdTcells, c( "gd_cluster" ) , assay = "Protein" , 
                     noaxis = F, do.label = F,
                    othertheme = list(coord_fixed(), theme(legend.position = "bottom"))) %T>% print()
 


copdratio <-  (table(COPD_gdTcells$disease, COPD_gdTcells$gd_cluster) %>%  
                 prop.table(margin = 1)*100) %>% data.frame() %>% filter(!is.na(Freq)) %>% 
  `colnames<-`(c("disease", "gd_cluster","Percent"))


library(ggalluvial)

F6C <-  copdratio %>% ggplot(aes(x = disease, y = Percent, fill = gd_cluster,stratum = gd_cluster)) +
  geom_stratum(size = 0.2)+fill_m()+theme_bw()+mytheme+NoLegend()+
  theme(axis.text.x = element_text(angle = 90, size = 8,hjust = 1, vjust = 0.5))
 
F6C

F6BC <- PG(list(F6B, F6C), ncol = 2, rw = c(1, 0.4),labels = c("B", "C")) %T>%  print()



F6ABC <- PG(list(F6A,F6BC), ncol = 1, rh = c(1,1.5)) %T>%  print()

COPD_gdTcells$disease %>% table()

table(COPD_gdTcells$orig.ident,COPD_gdTcells$disease )

COPD_gdTcells$gd_cluster %>% table()

Idents(COPD_gdTcells) <-  "gd_cluster"

ClusterCompare(COPD_gdTcells, assay = "Protein", id1 = "gdTRM_1", id2 = "gdTemra_1", test = "wilcox")
ClusterCompare(COPD_gdTcells, assay = "Protein", id1 = "Type1_Vd2", id2 = "Type3_Vd2", test = "wilcox")


F6D <-  Feature_rast(COPD_gdTcells, c(   "CD103-or-IntegrinalphaE-prot"  , "CD49a-prot",
                                      "CD45RA-prot", "CD45RO-prot",    "CD27-prot",  
                                    "TCRVdelta2-prot" ,
                                      
                             "KLRG1-or-MAFA-prot",  
                              "CD26-prot" ,"CD16-prot" , "CD57Recombinant-prot"
                              ) , assay = "Protein" , 
             colorgrd = "grd2", ncol = 2, sz = 0.4,
             othertheme = list(coord_fixed(),theme(
               plot.margin = margin(0,0,-2,-3, "pt"),
               legend.margin = margin(0,0,0,-10, "pt"))  )   ) %T>% print()



Feature_rast(COPD_gdTcells, "CD196-or-CCR6-prot", assay = "Protein")
Feature_rast(COPD_gdTcells, "CCR6")

F6A_D <- PG(list(F6ABC, F6D), rw = c(1.4,1), labels = c(NA, "D"))


ViolinPlot(COPD_gdTcells,c(names(sigtable)[c(1:9, 14,15)]), colors = umap.colors, box = T, group.by = "gd_cluster")

DotPlot(COPD_gdTcells %>% subset(gd_cluster %in% c("gdTRM_1", "gdTemra_1", "Type3_Vd2")),assay = 'RNA',dot.scale = 3.5,
        features =Dotfeatures)+heattheme

F6E <- (DotPlot(COPD_gdTcells %>% subset(gd_cluster %in% c("gdTRM_1", "gdTemra_1", "Type3_Vd2")),assay = 'RNA',dot.scale = 3.5,
               features =Dotfeatures)+
mytheme+
  theme(text = element_text(size = 8), 
        axis.text.y = element_text(size = 8),
        axis.line.y.left = element_line(),
        axis.text.x.bottom = element_text(size = 8, angle = 90, face = "italic", vjust = 0.5, hjust = 1),
        legend.key.height     = unit(3,'mm'),
        legend.key.width   = unit(5,'mm'),
        legend.box.margin = margin(5,5,0,0,unit = 'mm'),
        legend.box = "vertical",legend.position = 'right',
        axis.title = element_blank())+
  scale_y_discrete(position = 'left')+
  scale_x_discrete(position = 'bottom')+
  xlab(NULL)+ylab(NULL)+
  scale_color_gradient2(low = '#003399', mid = 'white',  high = "#990000")+
  guides(
    color = guide_colorbar(title.position = 'top',direction = 'horizontal',
    ),  size = guide_legend(title.position = 'top',direction = 'vertical',
                        label.position = 'right'))+
     NULL)%T>% print()

F6F <- Feature_rast(COPD_gdTcells, c("STAT1", "NKG7", "EOMES", "NR3C1", "ENTPD1", "GATA3","RORC","CCR6", "DPP4"  ), ncol = 3, sz = 0.4,
                    othertheme = list(coord_fixed(),theme(
                      plot.margin = margin(0,0,-2,-3, "pt"),
                      legend.margin = margin(0,0,0,-10, "pt"))  )  ) %T>% print()


Feature_rast(COPD_gdTcells, c("IRF2", "IRF4", "BATF"))

F6F_des <- Feature_density(COPD_gdTcells, c( "EOMES", "GATA3","RORC"), ncol = 3, 
                # sz = 0.4,
             othertheme = list(coord_fixed(),theme(
               plot.margin = margin(0,0,-2,-3, "pt"),
               legend.position =  "none"  )  )  ) %T>% print()




library(ggpubr)

GM_DM <- ViolinPlot(COPD_gdTcells,"GM_D", 
                    othertheme =   list(theme(axis.text.x =element_text(angle = 90, size = 8,hjust = 1, vjust = 0.5),
                                              plot.title = element_text(size = 8)), 
                                        stat_compare_means(method = "kruskal.test", size = gs(6) ),
                              
                                        ylab(NULL), ggtitle("Vd2 Type3 immunity")) ,colors = umap.colors, box = T)%T>% print()




# GM_FM <- ViolinPlot(COPD_gdTcells,"GM_F", 
#                     othertheme =   list(theme(axis.text.x = element_blank(),stat_compare_means(method = "kruskal.test", size = gs(6)),
#                                               plot.title = element_text(size = 6)),
#                                         ylab(NULL), ggtitle("CTL response (Vd2)")) ,colors = umap.colors, box = T)%T>% print()




GM_GM <-  ViolinPlot(COPD_gdTcells,"GM_G", 
                               othertheme =   list(theme(axis.text.x = element_text(angle = 90, size = 8,hjust = 1, vjust = 0.5),
                                                         plot.title = element_text(size = 6)),stat_compare_means(method = "kruskal.test", size = gs(6)),
                                                   ylab(NULL), ggtitle("CTL response ")) ,colors = umap.colors, box = T)%T>% print()


TregM <-  ViolinPlot(COPD_gdTcells,"Tregs", 
           othertheme =   list(theme(axis.text.x = element_text(angle = 90, size = 8,hjust = 1, vjust = 0.5),
                                     plot.title = element_text(size = 6)),stat_compare_means(method = "kruskal.test", size = gs(6)),
                               ylab(NULL), ggtitle("Treg module")) ,colors = umap.colors, box = T)%T>% print()




F6G <-  PG(list(GM_DM,  GM_GM, TregM), ncol = 3, labels = "Gene Modular Score",label_fontface = "plain",label_x = 0.8
         ) %T>% print()


  F6FG <- 
PG(list(F6F_des, F6G), ncol = 2, rw = c(1,1.2), labels = c("F", "G")) %T>% print()



F6E


gseatrm <-  GSEA_multipplot(GSEA_c0_c1, c("GOBP_INNATE_IMMUNE_RESPONSE",
                                          'GOBP_CELL_KILLING' ,
                                          "GSE7852_TREG_VS_TCONV_UP",
                                          'GOBP_BRANCHING_MORPHOGENESIS_OF_AN_EPITHELIAL_TUBE'),title = 'GSEA: gdTRM_1 vs gdTemra_1', legendpvalue = T, legend.position = "no",
                            plots = 1:2, rel_h = c(1,0.25),
                            
                            c1 = "gdTRM_1", c2 =  "gdTemra_1", base_size = 8 ) %T>% print()


F6I <-   PG(list(gseatrm, NA), ncol =2, rw = c(1, 1), labels = c("I", NA)) %T>% print()
 
F6E_I <- PG(list(F6E, F6FG, F6I), ncol = 1, labels = "E", rh = c(1.2,1,1.2)) %T>% print()

newF6 <-  PG(list(F6A_D, F6E_I), ncol = 1, rh = c(1,1.3)) %T>% figsave( path = figpath_ni, "Fig4_PublicCOPD_2025_new.pdf", 200, 250) 

proteinlist


table(
COPD_gdTcells@meta.data$disease,
COPD_gdTcells@meta.data$orig.ident)

Feature_rast(COPD_gdTcells, c("gd_cluster", "disease"))

COPD_gdTcells$orig.ident <-  paste0(COPD_gdTcells$disease, "_", COPD_gdTcells$orig.ident)
 
Feature_rast(COPD_gdTcells, c( "disease"), facets = "orig.ident", do.label = F, sz = 1.5)
Feature_rast(COPD_gdTcells, c( "gd_cluster"), facets = "orig.ident", do.label = F, sz = 1.5)


Feature_rast(COPD_Tcells, colorgrd = "gg")
Feature_rast(COPD_gdTcells, c("ident", "NR3C1", "RBPJ", "STAT1", "IRF4"), ncol = 2 )
