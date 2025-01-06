library(sceasy)
library(anndata)
library(SeuratData)
library(SeuratDisk)


ad <- anndata::read_h5ad('public/COVID_onlyT.h5ad')
ad_path <- "public/COVID_onlyT.h5ad"


convertFormat(ad_path, from="anndata", to="seurat", outFile='public/COVID_onlyT.rds')




library(reticulate)
ad <- import("anndata", convert = FALSE)
covidT_ad <- ad$read_h5ad('public/COVID_onlyT.h5ad')
Convert(ad_path,  dest = "h5seurat", overwrite = TRUE)
COVID_T <-LoadH5Seurat("public/COVID_onlyT.h5seurat")
COVID_T


library(zellkonverter)


ad <- readH5AD('public/COVID_onlyT.h5ad')
# Then you can use Seurat's function as.Seurat() to convert your object to Seurat. I also had to specify the default parameter counts and data to fit my data. E.g. i had to specify
COVID_T <- as.Seurat(ad, counts = "X", data = NULL)
COVID_T@assays$originalexp
COVID_T@assays$originalexp@counts

DefaultAssay(COVID_T) <-  "originalexp"

COVID_T@reductions$X_tsne


COVID_T  %<>% AddModuleScore(features = TRDlist, name = "TRD_score")
COVID_T  %<>% AddModuleScore(features = TRABlist, name = "TRAB_score")
COVID_T  %<>% AddModuleScore(features = CD3list, name = "CD3_score")



Feature_rast(COVID_T, c( "TRD_score1", "TRAB_score1") , d1 = "Xtsne_1", d2 = "Xtsne_2")
COVID_T@meta.data

Feature_rast(COVID_T, d1 = "TRAB_score1", d2 = "TRD_score1", sz = 0.2,
             noaxis = F, axis.number = T, do.label = F, g = "highlight",
             colorset = c(alpha("lightgrey",0.5), "red") )+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0.018, linewidth = 0.2)

Feature_rast(COVID_T, d1 = "TRAB_score1", d2 = "TRD_score1", sz = 0.2,
             noaxis = F, axis.number = T, do.label = F, g = "majorType" )+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0.018, linewidth = 0.2)


ViolinPlot(COVID_T, g = "CD3_score1",group.by = "majorType")+
  geom_hline(yintercept = -0.2, linewidth = 0.2)



ViolinPlot(COVID_T, g = "CD3_score1",group.by = "celltype", x.angle = 90)+
  geom_hline(yintercept = -0.35, linewidth = 0.2)



# purification
COVID_T <- subset(COVID_T, CD3_score1 > -0.35)

Feature_rast(COVID_T, d1 = "TRAB_score1", d2 = "TRD_score1", sz = 0.2,
             noaxis = F, axis.number = T, do.label = F, g = "TRDV2" )+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0.018, linewidth = 0.2)



Feature_rast(COVID_T, d1 = "TRAB_score1", d2 = "TRD_score1", sz = 0.2,
             noaxis = F, axis.number = T, do.label = F, g = "highlight",
             colorset = c(alpha("lightgrey",0.5), "red") )+
  geom_hline(yintercept = 0.23, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)


COVID_T@meta.data   %<>% mutate(gdTcells = case_when( 
  TRAB_score1 <= 0& TRD_score1 >=  0.2  ~ TRUE, TRUE ~ FALSE
))


Feature_rast(COVID_T, d1 = "TRAB_score1", d2 = "TRD_score1",
             noaxis = F, axis.number = T, do.label = F, g = "ITGA1",
             colorset = c(alpha("lightgrey",0.5), "red") )+
  geom_hline(yintercept = 0.2, linewidth = 0.2)+
  geom_vline(xintercept = 0, linewidth = 0.2)


colnames(COVID_T@meta.data)

COVID_T$gdTcells %>% table()
  
table(COVID_T$gdTcells, COVID_T$highlight)


ViolinPlot(COVID_T, g = "TRD_score1",group.by = "highlight", x.angle = 90)+
  geom_hline(yintercept = 0.2, linewidth = 0.2)



saveRDS(COVID_T, "public/COVID_ZeminZhang_allT.rds")


# select gdT  -------------------------------------------------------------




COVID_GDT <- subset(COVID_T,   gdTcells == TRUE)

saveRDS(COVID_GDT, "public/COVID_ZeminZhang_GDT.rds")

rm(COVID_GDT)
table(COVID_GDT$datasets)

gc(full = T , reset = T)


COVID_GDT <- subset(COVID_GDT,   gdTcells == TRUE & datasets != "d12" & datasets != "d10")


COVID_GDT <- readRDS("public/COVID_ZeminZhang_GDT.rds")


#integration amd reclustering  -----------------------------------------------------------
COVID_GDT[["originalexp"]] <- split(COVID_GDT[["originalexp"]], f = COVID_GDT$datasets)

COVID_GDT  %<>%   ScaleData( assay = 'originalexp',
                             vars.to.regress = c('nFeature_originalexp',"datasets")) %>%  
  FindVariableFeatures(assay = 'originalexp',nfeatures = 3000, selection.method = 'vst')

COVID_GDT@assays$originalexp@meta.data


COVID_GDT@assays$originalexp@var.features <- COVID_GDT@assays$originalexp@var.features%>%  
  str_subset('^RP|^MT|TRAV|TRBV|NO-NAME|IGL|IGH|LINC|MIR', negate = T)


COVID_GDT <- RunPCA(COVID_GDT, npcs = 100, verbose =
                      T, nfeatures.print = 40)


table(COVID_GDT$SARS.CoV.2)




COVID_GDT <- IntegrateLayers(object = COVID_GDT, method = CCAIntegration,  assay = "originalexp", 

                             orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
COVID_GDT[["originalexp"]] <- JoinLayers(COVID_GDT[["originalexp"]])

COVID_GDT <- FindNeighbors(COVID_GDT, reduction = "integrated.cca", dims = 1:50)
COVID_GDT <- FindClusters(COVID_GDT, resolution = 0.5)
COVID_GDT <- RunUMAP(COVID_GDT, dims = 1:50, reduction.name = "UMAP",reduction.key = "UMAP_",
                     reduction = 'integrated.cca')



Feature_rast(COVID_GDT, facets = "CoVID.19.severity")
Feature_rast(COVID_GDT, c("ident", "datasets", "Sample.type", "ITGAE"), ncol = 2)


Feature_rast(COVID_GDT, c("ident", "TRDV1", "TRDV2", "CoVID.19.severity", "Outcome", "SARS.CoV.2"), ncol = 3)
COVID_GDT$SARS.CoV.2

COVID_GDT$nCount_originalexp

COVID_GDT$CoVID.19.severity
COVID_GDT$Outcome

table(COVID_GDT$datasets, COVID_GDT$Sample.type)

COVID_GDT  %<>%   ScaleData( assay = 'originalexp',
                                 vars.to.regress = c('nFeature_originalexp',"datasets")) %>%  
  FindVariableFeatures(assay = 'originalexp',nfeatures = 3000, selection.method = 'vst')
COVID_GDT@assays$originalexp@var.features
COVID_GDT@assays$originalexp@var.features <- COVID_GDT@assays$originalexp@var.features%>%  
  str_subset('^RP|^MT|TRAV|TRBV|NO-NAME|IGL|IGH|LINC|MIR', negate = T)


COVID_GDT <- RunPCA(COVID_GDT, npcs = 100, verbose =
                          T, nfeatures.print = 40)

COVID_GDT$SARS.CoV.2

ElbowPlot(COVID_GDT, ndims = 100)

COVID_GDT  <- JackStraw(COVID_GDT, num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:50) 

JackStrawPlot(COVID_GDT, dims = 1:50 )

COVID_GDT <- RunUMAP(COVID_GDT, dims = 1:50, reduction.name = "UMAP",reduction.key = "UMAP_",
                         reduction = 'pca') %>%
  FindNeighbors(dims = c(1:50))

COVID_GDT <- FindClusters(COVID_GDT, resolution = 0.7)
Feature_rast(COVID_GDT, c( "Outcome"), facets = "Sample.type")

COVID_GDT$Outcome
COVID_GDT$SARS.CoV.2

ViolinPlot(COVID_GDT, "IFNG", group.by = "Outcome")
COVID_GDT$CoVID.19.severity

table(COVID_GDT$Sample.type, COVID_GDT$datasets)

COVID_GDT@meta.data %>% filter(datasets== "d07") %>% 
  count(PatientID, Sample.type)
COVID_GDT$CoVID.19.severity %>% unique()

ViolinPlot(COVID_GDT %>%  subset(Sample.type == "fresh BALF"), g = c("AREG", "IFNG", "ITGAE", "GATA3","GZMA",
                                                                     "ITGA1", "KLRG1", "CSF1"), 
           group.by = "CoVID.19.severity", sz = 2, colors = umap.colors)



ViolinPlot(COVID_GDT %>%  subset(datasets == "d07"), g = c("AREG", "IFNG", "ITGAE", 
                                                                     "ITGA1", "KLRG1", "CSF1"), 
           group.by = "CoVID.19.severity", colors = umap.colors, sz= 1)

subset(COVID_GDT, Sample.type == "fresh BALF")$PatientID %>% table()


COVID_M_S_heat <-  ClusterCompare(subset(COVID_GDT, Sample.type == "fresh BALF"),id1 = "mild/moderate",
               id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp")

COVID_M_S_heat$plot

view(COVID_M_S_heat$table)

ClusterCompare(subset(COVID_GDT, Sample.type == "fresh BALF"),id1 = "mild/moderate",
               id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp")


table(COVID_GDT$Sample.type, COVID_GDT$datasets)


ClusterCompare(subset(COVID_GDT, datasets == "d07"),id1 = "mild/moderate",
               id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp")

d07_pbmc <-  
ClusterCompare(subset(COVID_GDT, datasets == "d07" & Sample.type == "fresh PBMC"),id1 = "mild/moderate", log2fc = 0,
               id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp")

d07_pbmc$table
d07_balf <-  
  ClusterCompare(subset(COVID_GDT, datasets == "d07" & Sample.type == "fresh BALF"),id1 = "mild/moderate", log2fc = 0,
                 id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp")


d07_balf

d07_balf$table
ClusterCompare(COVID_GDT,id1 = "mild/moderate",
               id2 = "severe/critical", group.by = "CoVID.19.severity", assay = "originalexp",genetoshow = 5)
COPD_gdTcells$disease

ClusterCompare(COPD_gdTcells, id1 = "Healthy", id2 =  "Donor", group.by = "disease")
Idents(COVID_GDT) <-  "CoVID.19.severity"

Feature_rast(COVID_GDT)


# GSEA

genelist_COVID_M_S <- Genelist_generator(subset(COVID_GDT, Sample.type == "fresh BALF"),
                                     "mild/moderate", "severe/critical", group.by = "CoVID.19.severity")

# genelist_Type1_3_Vd2_c1 <- Genelist_generator(COPD_gdTcells, "Type1_Vd2", "Type3_Vd2")

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



GSEA_COVID_M_S <- GSEA(geneList = genelist_COVID_M_S, TERM2GENE=ALL_msigdb_G,
                   # minGSSize    = 10,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH") 

view(GSEA_COVID_M_S@result)



GSEA_multipplot(GSEA_COVID_M_S, description_to_show = 
c("GOBP_GRANULOCYTE_CHEMOTAXIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), c1 = "COVID_moderate", c2 = "COVID_severe")


Feature_rast(subset(COVID_GDT, datasets == "d07" ), c("Sample.type", "CoVID.19.severity", "TRDV1", "TRDV2"))

