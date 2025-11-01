
# IL32 PBMC stimulation ---------------------------------------------------



library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(tidyr)

library(reticulate)
library(anndata)

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
source('/home/big/tanlikai/script/rscripts/funcs.r')

file.edit("/home/big/tanlikai/script/rscripts/funcs.r")
# IL32_PBMC_T cells

IL32_T_H5AD <-  "/home/big/public_datasets/parse_10m/Parse_10M_PBMC_cytokines_t_pbs_il32_beta_ds10000.h5ad"
IL32_PB_H5AD <- "/home/big/public_datasets/parse_10m/Parse_10M_PBMC_cytokines_PBMC_pbs_il32_beta_ds10000.h5ad"

IL32_T_RDS <-  "/home/big/public_datasets/parse_10m/Parse_10M_PBMC_cytokines_t_pbs_il32_beta_ds10000.RDS"

IL32_PB_RDS <- "/home/big/public_datasets/parse_10m/Parse_10M_PBMC_cytokines_PBMC_pbs_il32_beta_ds10000.RDS"

T32co <-readRDS(IL32_T_RDS)
T32co$RNA_snn_res.0.8
# read H5AD T cells ---------------------------------------------------------------


py_install("anndata", upgrade = TRUE)

# counts_matrix <- t(ad$X)
# head(counts_matrix)


ad$layers

ad <-  read_h5ad(IL32_T_H5AD)
cell_metadata <- ad$obs %T>% head(.) %T>% print()
gene_metadata <- ad$var %T>% head(.) %T>% print()


genes <-  ad$var_names
barcodes <-  ad$obs_names

barcodes

# read raw counts

IL32_T_H5AD <-  "/home/big/public_datasets/parse_10m/Parse_10M_PBMC_cytokines_t_pbs_il32_beta_ds10000.h5ad"

mtx_data <- data.table::fread("/home/big/public_datasets/parse_10m/IL32Tcounts_ds.mtx", 
                              skip = 3, col.names = c("i_row_idx", "j_col_idx", "count"))

head(mtx_data)

dims <- as.numeric(strsplit(readLines("/home/big/public_datasets/parse_10m/IL32Tcounts_ds.mtx",
                                      n = 3)[3], " ")[[1]])
n_cells <- dims[1]
n_genes <- dims[2]
n_genes
n_cells
dims
n_genes
raw_counts_matrix <- sparseMatrix(
  i = mtx_data$i_row_idx,
  j = mtx_data$j_col_idx,
  x = mtx_data$count,
  dims = c(n_cells, n_genes),
  dimnames = list(barcodes, genes)
)




head(counts_matrix)


T32co<- CreateSeuratObject(counts = t(raw_counts_matrix),   assay = "RNA",  
                             min.cells = 200,
                             meta.data = cell_metadata)


T32co  %<>% NormalizeData()

ViolinPlot(T32co, "ITGA1", group.by = "cytokine", box = T)

DEG_allT_IL32  <- ClusterCompare(T32co, group.by = "cytokine", id1 = "IL-32-beta", id2 = "PBS")



# read H5AD PBMC cells ---------------------------------------------------------------


ad <-  read_h5ad(IL32_PB_H5AD)

counts_matrix <- t(ad$X)
head(counts_matrix)

cell_metadata <- ad$obs %T>% head(.) %T>% print()
gene_metadata <- ad$var %T>% head(.) %T>% print()
head(gene_metadata)
# Create a Seurat V5 object
# For a basic conversion, you can just use the counts.
# For a more complete conversion, you can add the metadata.
PB_32co<- CreateSeuratObject(counts = counts_matrix,
                                meta.data = cell_metadata)

saveRDS(PB_32co, IL32_PB_RDS)

# Add gene metadata to the Seurat object
# PB_32co[['RNA']]$data <- counts_matrix


rm(PB_32co)




# QC ----------------------------------------------------------------------

T32co[["percent.mt"]] <- PercentageFeatureSet(T32co, pattern = "^MT-")


T32co$nCount_RNA

T32co$sample

Feature_rast(T32co, d1 = "nFeature_RNA", d2 = "nCount_RNA", g=  "percent.mt", noaxis = F, axis.number = T, colorgrd = "grd2")

ViolinPlot(T32co, "nCount_RNA", group.by = "cytokine", box = T)


Feature_rast(T32co, d1 = "nCount_RNA", d2 = "CD4", g=  "cytokine", noaxis = F, axis.number = T, colorset = "gg")

T32co <- subset(T32co, subset = nFeature_RNA > 800 & nFeature_RNA < 4500 & percent.mt < 5) 
T32co  %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) 

variable_genes_filtered <- T32co[["RNA"]]@meta.data$var.features %>%  na.omit() %>% as.vector() %T>% print()


# Identify gene sets to remove from downstream analysis (PCA, clustering)
# Mitochondrial genes
mt_genes <- grep(pattern = "^MT-", x = rownames(T32co), value = TRUE)
# Ribosomal genes
ribo_genes <- grep(pattern = "^RP[SL]", x = rownames(T32co), value = TRUE)
# T-cell receptor genes (example, may need adjustment for your data)
tcr_genes <- grep(pattern = "^TR[ABVJD]", x = rownames(T32co), value = TRUE)

IG_genes <- grep(pattern = "^IG[HLSK]", x = rownames(T32co), value = TRUE) %T>% print()


ENS_genes <- grep(pattern = "^ENS|LINC|HIST|HSP|AS1", x = rownames(T32co), value = TRUE) %T>% print()



# Remove these unwanted genes from the variable features list

variable_genes_filtered <- variable_genes_filtered[!variable_genes_filtered %in% c(mt_genes, ribo_genes, tcr_genes, IG_genes,ENS_genes)]

variable_genes_filtered

T32co$donor

# This step is performed ONLY on the filtered list of variable genes.
T32co <- ScaleData(
  T32co,
  features = c(variable_genes_filtered,DEG_allT_IL32$table$gene) %>%  unique(),
  vars.to.regress = c("percent.mt", "nFeature_RNA", "donor")
)

DEG_allT_IL32  <- ClusterCompare(T32co, group.by = "cytokine", id1 = "IL-32-beta", id2 = "PBS")
DEG_allT_IL32$table %>% filter(gene %in% c("ITGAE", "ITGA1", "ZNF683", "EOMES", "CXCR6", "CTLA4", "AREG", "GATA3", "BATF", "CSF1", "CSF2", "GZMB", "GZMA", "JUNB", "FOXP3", "AREG"))
DEG_allT_IL32$plot
DEG_allT_IL32$table$gene
## Step 4: Perform Dimensionality Reduction (PCA) 📉
## -----------------------------------------------------------
saveRDS(T32co, IL32_T_RDS)

# Run PCA using only the clean set of variable genes
T32co <- RunPCA(
  T32co,
  features = variable_genes_filtered, npcs = 100, ndims.print = 30, seed.use = 100
)

print(T32co[["pca"]], dims = 1:30, nfeatures = 30)


  # Visualize PCA results to determine the number of dimensions for clustering
ElbowPlot(T32co, ndims = 100)
T32co %>% 
  JackStraw( num.replicate = 100, dims = 50)%>%
  ScoreJackStraw(dims = 1:80) 


JackStrawPlot(T32co, dims = 1:50 ) %>%
  figsave('T32co.jackstraw.pdf' , w = 400, h = 400)

# dimensional reduction & clustering --------------------------------------

T32co  %<>%  RunUMAP( dims = 1:40, 
                          reduction = 'pca', reduction.key = "UMAP",seed.use = 1) %>%
  FindNeighbors(dims = c(1:40))




for (i in seq(0.8,1.5,0.1) ) {
  T32co  %<>% FindClusters( resolution = i)
}
T32co  %<>% FindClusters( resolution = 0.6)


Feature_rast(T32co, colorset = "gg", facets = "cytokine")

Feature_rast(T32co, c("TRDV1", "TRDV2",  "TRDC", "TRGC1", "TRGC2", "ident"), colorset = "gg")

Feature_rast(T32co, c("CD4", "CD8A","FOXP3", "GZMB", "GZMA", "EOMES", "IFNG","ITGA1", "ITGAE",  "KLRG1", "KLRD1", "KLF2","TRDV1", "TRDV2", "CCR6", "RORC", "CCR4", "CCR7", "CD27"), sz = 0.2)


Feature_rast(T32co, c("CD4", "CD8A", "TRDC", "ITGA1", "ITGAE", "ZNF683", "CXCR6" , "JUNB", "BATF", "CTLA4"), facets = "cytokine", sz = 0.1)

Feature_rast(T32co, "ITGAE", facets = c("donor","cytokine" ), sz = 0.3, facetcol = 6)

T32co$RNA_snn_res.1.5

T32co <- T32co %>%  SetIdent(value = "RNA_snn_res.0.6")

Feature_rast(T32co, "RNA_snn_res.0.6", colorset = "gg")


Feature_rast(T32co, facets = "cytokine")

saveRDS(T32co, IL32_T_RDS)

ClusterCompare(T32co, "15", "10")

# gene modules  -----------------------------------------------------------



sigtable <- openxlsx::read.xlsx('abt/abd5778_Table_S3.xlsx') %>% 
  `colnames<-`(str_remove(colnames(.), '.\\(.+\\)') ) %>%  as.list() %>% 
  map(~ na.exclude(.x) %>% as.vector) %T>% print()

names(sigtable)
sigtable$Tissue.resident <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "RBPJ", "GPR25")

TRM <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "RBPJ", "GPR25")

#caculate scores
for (i in names(sigtable)) {
  T32co <- AddModuleScore(T32co, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}


T32co <- AddModuleScore(T32co, features = list(TRM), name = "TRM_score", assay = 'RNA')


colnames(T32co@meta.data) %<>% str_replace("(?<=\\w)1$", '')
T32co@meta.data


Feature_rast(T32co,c(names(sigtable), "TRM_score1"), colorgrd = "grd2", sz = 0.2, ncol = 4)

ViolinPlot(T32co, "TRM_score1", group.by = "cytokine", box = T,sz = 0)

ViolinPlot(T32co, "TRM_score",  box = T,sz = 0)

VlnPlot(T32co, "ITGA1", pt.size = 0,split.by = "cytokine")

VlnPlot(T32co, "TRM_score1", pt.size = 0,split.by = "cytokine", group.by = "RNA_snn_res.0.8")
table(T32co$cytokine, T32co$RNA_snn_res.0.8) %>% prop.table(margin = 1)*100




Feature_rast(T32co, d1 = "ITGA1", d2 = "TRM_score1", g = "ITGAE", noaxis = F, axis.number = T)

Feature_rast(T32co, TRM)

T32co$bc <-  colnames(T32co)

bc_TRM <- colnames(   
  subset(T32co, ITGA1 > 0.6 &  ITGAE >  0.6 )
)

T32co@meta.data   %<>% mutate(CD103_CD49a_TRM = case_when(bc %in% bc_TRM ~ "CD103_CD49a_TRM", !is.na(bc) ~ "other")) 

T32co@meta.data   %<>% mutate(TRMs =  case_when(TRM_score > 0.4 ~ "TRM", TRM_score <= 0.4 ~ "NonTRM")) 


Feature_rast(T32co, "CD103_CD49a_TRM", facets = "cytokine", navalue = "lightgrey", sz = 0.3)

table(T32co$sample, T32co$CD103_CD49a_TRM) %>% prop.table(margin = 1)*100

ClusterCompare(T32co, id1 = "CD103_CD49a_TRM", id2 = "other", group.by = "CD103_CD49a_TRM")

ClusterCompare(T32co, id1 = "TRM", id2 = "NonTRM", group.by = "TRMs")

table(T32co$sample, T32co$TRMs) %>% prop.table(margin = 1)*100

Feature_rast(T32co, facets = c("CD103_CD49a_TRM", "cytokine"))


T32co$TRM_Cytokin <-  paste0(T32co$TRMs, "_",T32co$cytokine)
T32co$ITGAE1_Cytokine  <-  paste0(T32co$CD103_CD49a_TRM, "_",T32co$cytokine)


T32co$cytokine
saveRDS(T32co, IL32_T_RDS)


ClusterCompare(T32co, "TRM_IL-32-beta", "NonTRM_PBS", group.by = "TRM_Cytokin")


ClusterCompare(T32co, "TRM_IL-32-beta", "TRM_PBS", group.by = "TRM_Cytokin")


ClusterCompare(T32co, "10", "12", group.by = "RNA_snn_res.1")


ClusterCompare(T32co, "10", "3")

ViolinPlot(T32co,c("DPP4", "KLRB1", "RORC", "CCR6"), group.by = "RNA_snn_res.1", colors = ggplotColours(28))

ViolinPlot(T32co,c("Th1"), colors = ggplotColours(28))

T32co$RNA_snn_res.1

Feature_rast(T32co, c("CD4", "CD8A","FOXP3", "GZMB", "GZMA", "EOMES", "IFNG","ITGA1", "ITGAE",  "KLRG1", "KLRD1", "KLF2","TRDV1", "TRDV2", "TRDC", "CCR6", "RORC", "CCR4", "CCR7", "CD27"), sz = 0.2)

multicores()

T32co$CD103_CD49a_TRM



Feature_rast(T32co, c("RNA_snn_res.1", "CCR6", "DPP4", "KLRB1", "RORC", "CD4", "CD8A", "CCR4", "GATA3"), colorset = 'gg', ncol = 3)




T32co$RNA_snn_res.0.6 <-  factor(T32co$RNA_snn_res.0.6,  levels= 0:19)
Feature_rast(T32co, "CD8B", d1 = "CD4", d2 = "CD8A")

Feature_rast(T32co, "RNA_snn_res.0.6")
annoframe <-   data.frame(
  RNA_snn_res.0.6 = 0:19,
  Cell_pheno = c(
    "Naïve_CD4", "Navie_CD8", "Th2", "TEM", "Naïve_CD4", 
    "Th2", "Th2", "Treg", "TCM_CD8", "CD8pos_gdT",
    "TCM", "Naïve_CD4", "Naïve_CD4", "gdT", "Tc17",
    "exTRM_CD8", "Teminal_CD8", "Navie_CD8", "Uidentiffied", "Uidentiffied"
  )
)


Cell_pheno = c(
  "Naïve_CD4", "Navie_CD8", "Th2", "TEM", "Naïve_CD4", 
  "Th2", "Th2", "Treg", "TCM_CD8", "CD8pos_gdT",
  "TCM", "Naïve_CD4", "Naïve_CD4", "gdT", "Tc17",
  "exTRM_CD8", "Teminal_CD8", "Navie_CD8", "Uidentiffied", "Uidentiffied"
) %>%  `names<-`(as.character(0:19))

T32co$barcodes <-  colnames(T32co)


T32co$Cell_pheno <- Cell_pheno[T32co$RNA_snn_res.0.6] %>% as.vector()


T32co@meta.data  %<>%  mutate(Cell_pheno = case_when(RNA_snn_res.1 == "10" ~ "Th17", !is.na(Cell_pheno) ~ Cell_pheno
                                                      ))



Feature_rast(T32co, c("Cell_pheno","RNA_snn_res.0.6"))

T32co %<>% SetIdent(value = "Cell_pheno")
saveRDS(T32co, IL32_T_RDS)

# GSEA --------------------------------------------------------------------



library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(GSEABase)


# Reference  datasets
ALL_msigdb_G  <- rbind(
  # c7 Immunology 
  msigdbr::msigdbr(species = "Homo sapiens", category = "C7"), 
  # hallmarker
  msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  # C2 KEGG
  msigdbr::msigdbr(species = "Homo sapiens", category = "C2",subcategory = 'CP:KEGG'),
  # GOBP
  msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')) %>% 
  dplyr::select(gs_name, gene_symbol)
ALL_msigdb_G


HALLMARKERS <-    msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% 
dplyr::select(gs_name, gene_symbol)


H_K_G  <- rbind(
  # # c7 Immunology 
  # msigdbr::msigdbr(species = "Homo sapiens", category = "C7"), 
  # hallmarker
  msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
  # C2 KEGG
  msigdbr::msigdbr(species = "Homo sapiens", category = "C2",subcategory = 'CP:KEGG'),
  # GOBP
  msigdbr::msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'GO:BP')) %>% 
  dplyr::select(gs_name, gene_symbol)


ALL_msigdb_G$gs_name %>% unique()   %>%  length()

T32co$CD103_CD49a_TRM
TRMvsNoTRM <- Genelist_generator(T32co, 'TRM', 'NonTRM',   group.by = "TRMs")
TRMvsNoTRM_E1 <- Genelist_generator(T32co, "CD103_CD49a_TRM" , 'other',   group.by = "CD103_CD49a_TRM")



TRMvsNoTRM_E1 <- Genelist_generator(T32co, "CD103_CD49a_TRM" , 'other',   group.by = "CD103_CD49a_TRM")


IL32_ctrl <- Genelist_generator(T32co, "IL-32-beta" , 'PBS',   group.by = "cytokine")




GSEA_TRMvsTEMRA_allref<-GSEA(geneList = TRMvsNoTRM, TERM2GENE=ALL_msigdb_G,  
                             pvalueCutoff = 0.05, pAdjustMethod = "BH") 


GSEA_TRMvs_H<-GSEA(geneList = TRMvsNoTRM, TERM2GENE=HALLMARKERS,  
                             pvalueCutoff = 0.05, pAdjustMethod = "BH")  %>% arrange(.@result, desc(NES))

GSEA_TRMvs_HKG<-GSEA(geneList = TRMvsNoTRM, TERM2GENE=H_K_G,  
                   pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% arrange(.@result, desc(NES))


GSEA_TRMvs_H@result

GSEA_IL32_H<-GSEA(geneList = IL32_ctrl, TERM2GENE=HALLMARKERS,  
                   pvalueCutoff = 0.05, pAdjustMethod = "BH")  %>% arrange(.@result, desc(NES))

GSEA_IL32_H@result


GSEA_multipplot(GSEA_TRMvs_H,description_to_show = c( "HALLMARK_TGF_BETA_SIGNALING","HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),legendpvalue = T, c1 = "TRM", c2 = "other" )


gseaplot2(GSEA_IL32_H,c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE"))

GSEA_multipplot(GSEA_IL32_H,description_to_show = c( "HALLMARK_TGF_BETA_SIGNALING"),legendpvalue = T, c1 = "IL32b", c2 = "PBS" )


GSEA_TRMvs_H@result %>%  filter(NES > 0)


GSEA_TRMvs_HKG@result  %>%  filter(NES > 0) %>% view()


GSEA_TRMvsTEMRA_allref@result %>%  filter(NES > 0) %>% view()


# propostion change -------------------------------------------------------
library(rstatix)

T32co$donor


cell_proportions <- T32co@meta.data %>%
  # Count cells for each group
  dplyr::count(donor, cytokine, Cell_pheno, name = "n") %>%
  complete(nesting(donor, cytokine), Cell_pheno, fill = list(n = 0)) %>%
  # Calculate the total number of cells for each donor/cytokine pair
  group_by(donor, cytokine) %>%
  
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  # Calculate the proportion as a percentage
  mutate(proportion = (n / total_cells) * 100)


cell_proportions



stat_test <- cell_proportions %>%
  group_by(Cell_pheno) %>%
  wilcox_test(proportion ~ cytokine, paired = TRUE) %>%
  # Adjust p-values for multiple comparisons
  adjust_pvalue(method = "fdr",output.col = "fdr") %>%
  add_significance("fdr") %>%
  # Add coordinates for placing p-values on the plot
  add_xy_position(x = "cytokine")

stat_test 



ggplot(cell_proportions, aes(x = cytokine, y = proportion)) +
  geom_line(aes(group = donor), color = "grey70", alpha = 0.8) +
  geom_point(aes(color = cytokine), size = 3.5, alpha = 0.9) +
  ggpubr::stat_pvalue_manual(
    stat_test,
    label = "fdr",
    tip.length = 0.01,
    bracket.nudge.y = 0.5,
    hide.ns = TRUE
  ) +
  facet_wrap(~ Cell_pheno, scales = "free_y", ncol = 5) +
  scale_color_manual(values = c("PBS" = "#0072B2", "IL-32-beta" = "#D55E00")) +
  labs(
    title = "Cell Type Proportion Changes Upon Cytokine Stimulation",
    x = "Condition",
    y = "Proportion of Cells (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+mytheme

ggpubr::stat_pvalue_manual()

# TRM proportion ----------------------------------------------------------

T32co$CD103_CD49a_TRM
# This creates a summary data frame. We assume the 'TRMs' column
# contains values like "TRM" and "Non_TRM" or is a logical TRUE/FALSE.
trm_proportions <- T32co@meta.data %>%

  # Group by sample and cell type to get counts
  group_by(donor, cytokine, Cell_pheno) %>%
  # Calculate the number of TRM cells and the total cells in each group
  summarise(
    trm_count = sum(CD103_CD49a_TRM == "CD103_CD49a_TRM", na.rm = TRUE), # Counts cells where TRMs is "TRM"
    total_cells = n(),
    .groups = "drop"
  ) %>%
  # Calculate the proportion as a percentage
  mutate(proportion = (trm_count / total_cells) * 100) %>% 
  complete(nesting(donor, cytokine), Cell_pheno, fill = list(proportion = 0))


totlaT <-  T32co@meta.data %>%
  
  # Group by sample and cell type to get counts
  group_by(donor, cytokine) %>%
  # Calculate the number of TRM cells and the total cells in each group
  summarise(
    trm_count = sum(CD103_CD49a_TRM == "CD103_CD49a_TRM", na.rm = TRUE), # Counts cells where TRMs is "TRM"
    total_cells = n(),
    .groups = "drop"
  ) %>%
  # Calculate the proportion as a percentage
  mutate(proportion = (trm_count / total_cells) * 100) 

totlaT$Cell_pheno <-  "TotalT"

totlaT <-  totlaT[, colnames(trm_proportions)]
totlaT  
trm_proportions <-  rbind(totlaT, trm_proportions)
allp <-  trm_proportions$Cell_pheno %>% unique()

trm_proportions$Cell_pheno <-  factor(trm_proportions$Cell_pheno, levels =  allp[c(1,4, 7, 2,3,5,6,8:15)] )

stat_test_2 <- trm_proportions %>%
  group_by(Cell_pheno) %>%
  wilcox_test(proportion ~ cytokine, paired = TRUE) %>%
  adjust_pvalue(method = "fdr", output.col = "fdr") %>%
  add_significance("fdr") %>%
  add_xy_position(x = "cytokine")

stat_test_2
stat_test_2$fdr <-  round(stat_test_2$fdr, digits = 4)
stat_test_2$y.position <-  stat_test_2$y.position - 3

ggplot(trm_proportions, aes(x = cytokine, y = proportion)) +
  geom_boxplot(fll = 'transparent')+
  # Draw lines connecting the same donor between conditions
  geom_line(aes(group = donor), color = "grey70", alpha = 0.8) +
  # Add points for each sample, colored by condition
  geom_point(aes(color = cytokine), size = 2, alpha = 0.9) +
  # Add the adjusted p-values from our statistical test
  ggpubr::stat_pvalue_manual(
    stat_test_2,
    label = "p",
    tip.length = 0.01,
    bracket.nudge.y = -1,
    hide.ns = T
  ) +
  # Create a separate panel for each cell type
  facet_wrap(~ Cell_pheno, scales = "free_y", ncol = 5) +
  # Set colors and labels
  scale_color_manual(values = c("PBS" = "#0072B2", "IL-32-beta" = "#D55E00")) +
  labs(
    title = "TRM Proportion Changes Upon IL-32B Stimulation",
    x = "Condition",
    y = "Proportion of TRM Cells (%)"
  ) +
  # Apply a clean theme
  theme_bw() + mytheme

ViolinPlot(T32co, "KLF2", group.by = "CD103_CD49a_TRM")


Feature_rast(T32co, c('ident', "cytokine", "ITGAE", "ITGA1", "ZNF683","CXCR6",  "KLRG1", "ITGAM", "KLF2"), ncol = 3, sz = 0.2)


Feature_rast(T32co, "CD103_CD49a_TRM", colorset = c('red', alpha("lightgreen", 0.3)))

# Treg and IL32 -----------------------------------------------------------


DEG_Treg <- ClusterCompare(T32co %>%  subset(Cell_pheno == "Treg"), id1 = "IL-32-beta", id2  = "PBS",group.by = "cytokine"    )
DEG_Treg$table
top50  <- DEG_Treg$table %>%  filter(avg_log2FC  >  0.5) %>%  pull(gene)
DEG_Treg$plot
T32co$donor
DEG_Treg$table %>%  filter(avg_log2FC  < -0.5) 
DoMultiBarHeatmap(T32co %>%  subset(Cell_pheno == "Treg", downsample = 200), features = top50, additional.group.by = "donor" , group.by = "cytokine") + scale_fill_gradient2(low = "blue", mid = "white", high = "red")

Treg_genelist <-  Genelist_generator(T32co %>%  subset(Cell_pheno == "Treg"), id1 = "IL-32-beta", id2  = "PBS",group.by = "cytokine" )

GSEA_TREG_HKG<-GSEA(geneList = Treg_genelist, TERM2GENE=H_K_G,  
                     pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% arrange(.@result, desc(NES))

GSEA_TREG_H<-GSEA(geneList = Treg_genelist, TERM2GENE=HALLMARKERS,  
                    pvalueCutoff = 0.05, pAdjustMethod = "BH") %>% arrange(.@result, desc(NES))

GSEA_TREG_HKG@result %>%  filter(NES > 0)

GSEA_TREG_H@result %>%  filter(NES < 0)


GSEA_multipplot(GSEA_TREG_H,description_to_show = c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),legendpvalue = F, c1 = "IL32b", c2 = "PBS" )

view(GSEA_TREG_H@result)


# IL32 stimulation FACS ---------------------------------------------------

IL32_treatment <-  read.csv("FACSdata/long_format_data_with_treatment.csv")
IL32_treatment  <-  read.csv("FACSdata/long_format_data_with_treatment.csv")  %>%
  mutate(treatment  = case_when(treatment == 1 ~ "blank",
                                treatment == 2 ~ "TGFb/IL15",
                                treatment == 3 ~ "IL32b", 
                                treatment == 4 ~ "IL32g", 
                                treatment == 5 ~ "IL32b/TGFb/IL15",
                                treatment == 6 ~ "IL32g/TGFb/IL15"
  )  ) %>%  mutate(treatment = 
                     factor(treatment, levels = 
                              c("blank","IL32b", "IL32g","TGFb/IL15", "IL32b/TGFb/IL15", "IL32g/TGFb/IL15" )))

IL32_treatment



IL32_treatment %>%  filter(day == 'D7' & Celltype == "CD4") %>% 
  ggplot(aes(x = treatment, y = CD103_percent)) + 
  ggtitle("CD103+ of CD4 T cells")+
  geom_line(aes(group = Donor), color = "grey")+
  geom_boxplot(aes(color = treatment), fill = 'transparent')+
  theme_classic()


IL32_treatment %>%  filter(day == 'D7' & Celltype == "CD8") %>% 
  ggplot(aes(x = treatment, y = CD103_percent)) +
  ggtitle("CD103+ of CD8 T cells")+
  geom_line(aes(group = Donor), color = "grey")+
  geom_boxplot(aes(color = treatment), fill = 'transparent')+
  theme_classic()

IL32_treatment %>%  filter(day == 'D7' & Celltype == "CD8") 



IL32_treatment %>%  filter(day == 'D7' & Celltype == "Vd1") %>% 
  ggplot(aes(x = treatment, y = CD103_percent)) + 
  ggtitle("CD103+ of Vd1 gdT cells")+
  geom_line(aes(group = Donor), color = "grey")+
  geom_boxplot(aes(color = treatment), fill = 'transparent')+
  theme_classic()

FC_IL32 %>% filter(Celltype == "Vd1")

FC_IL32 <-  IL32_treatment %>%
  filter(day == 'D7'  & treatment %in% c("blank", "TGFb/IL15", "IL32b/TGFb/IL15")) %>% 
  group_by(Donor, Celltype) %>% mutate( CD103_increase = CD103_percent - (CD103_percent[treatment == "blank"])) %>% 
  filter(treatment != "blank")


ggplot(FC_IL32 %>%  filter(Celltype %in% c('CD8', "Vd1")), aes(x = treatment, y = CD103_increase))+
  geom_line(aes(group = Donor), color = "grey")+
  geom_point(aes(color = treatment))+
  # geom_boxplot(aes(color = treatment), fill = 'transparent')+
  facet_wrap(~~Celltype, scale = "free_y")+
  ggpubr::stat_compare_means(paired = T, method = "t.test",label = "p.format")+
  theme_minimal()



# figures -----------------------------------------------------------------


 GSEA_IL32_TGFB <-  GSEA_multipplot(GSEA_IL32_H,description_to_show = c( "HALLMARK_TGF_BETA_SIGNALING"),legendpvalue = F, c1 = "IL32b", c2 = "PBS", plots =c(1,2)  )

GSEA_IL32_TGFB

IL32_umap <-  Feature_rast(T32co, "Cell_pheno", labelsize = 6)  + NoLegend()


Feature_rast(T32co,labelsize = 4)+ NoLegend()
Cellpop
Cellpop <-  ggplot(trm_proportions %>% filter(Cell_pheno %in% c("TotalT", 'gdT', 'Navie_CD8', "Tc17")), aes(x = cytokine, y = proportion)) +
  geom_line(aes(group = donor), color = "grey70", alpha = 0.8) +
  geom_point(aes(color = cytokine), size = 1.2, alpha = 0.9) +
  ggpubr::stat_pvalue_manual(
    stat_test_2 %>% filter(Cell_pheno %in% c("TotalT", 'gdT', 'Navie_CD8', "Tc17")),
    label = "p",
    tip.length = 0.01,
    # bracket.nudge.y = 0.5,
    hide.ns = TRUE
  ) +
  facet_wrap(~ Cell_pheno, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("PBS" = "#0072B2", "IL-32-beta" = "#D55E00")) +
  labs(
    title = "Cell Type Proportion Changes Upon Cytokine Stimulation",
    x = "Condition",
    y = "Proportion of Cells (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+mytheme

ggplot(trm_proportions %>% filter(Cell_pheno %in% c('gdT')), aes(x = cytokine, y = proportion)) +
  geom_line(aes(group = donor), color = "grey70", alpha = 0.8) +
  geom_point(aes(color = cytokine), size = 1.2, alpha = 0.9) +
  ggpubr::stat_pvalue_manual(
    stat_test_2 %>% filter(Cell_pheno %in% c('gdT' )),
    label = "p",
    tip.length = 0.01,
    # bracket.nudge.y = 0.5,
    hide.ns = TRUE
  ) +
  # facet_wrap(~ Cell_pheno, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("PBS" = "#0072B2", "IL-32-beta" = "#D55E00")) +
  labs(
    title = "TRMs in gdT cells ",
    x = "Condition",
    y = "Proportion of Cells (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




IL32_sti_facs <- ggplot(FC_IL32 %>%  filter(Celltype %in% c('CD8', "Vd1")), aes(x = treatment, y = CD103_increase))+
  geom_line(aes(group = Donor), color = "grey")+
  geom_point()+
  # geom_boxplot(aes(color = treatment), fill = 'transparent')+
  facet_wrap(~~Celltype, scale = "free_y" , ncol = 1)+
  ggpubr::stat_compare_means(paired = T, method = "t.test",label = "p.format")+
  theme_classic()+mytheme+theme(axis.text.x = element_text(size = 6 , angle = 45))

IL32_sti_facs

figsave(IL32_sti_facs, "IL32_sti_facs.pdf", path= figpath_ni, w  = 30, h = 70)


GSEA_IL32_TGFB

figsave(Cellpop, "IL32_TRM_sc.pdf", path= figpath_ni, w  = 30, h = 70)





# normal lung -------------------------------------------------------------
lungH5AD <-  "normallungT.h5ad"

mtx_data_lung <- data.table::fread("normallung.mtx", 
                              skip = 3, col.names = c("i_row_idx", "j_col_idx", "count"))

dims <- as.numeric(strsplit(readLines("normallung.mtx",
                                      n = 3)[3], " ")[[1]])


ad <-  read_h5ad(lungH5AD)
genes <-  ad$var_names
barcodes <-  ad$obs_names



cell_metadata <- ad$obs %T>% head(.) %T>% print()
gene_metadata <- ad$var %T>% head(.) %T>% print()

length(genes)

gene_translation <- bitr(
  geneID = genes,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)
length(gene_translation)
uniquegene <- gene_translation$SYMBOL %>% dis
gene_translation$ENSEMBL
geneID <-  gene_translation %>% count(SYMBOL) %>% filter(n == 1) %>% pull(SYMBOL)

esid_keep <-  gene_translation %>%  filter(SYMBOL %in% geneID) 
esid_keep
gene_translation %>% count(SYMBOL) %>% filter(n == 1) 

gene_translation %>%  filter(SYMBOL %in% geneID) 

dplyr::pull()
gene_translation %>% filter(!is.na(SYMBOL))



n_cells <- dims[1]
n_genes <- dims[2]
n_genes
n_cells
dims
n_genes
raw_counts_matrix <- sparseMatrix(
  i = mtx_data_lung$i_row_idx,
  j = mtx_data_lung$j_col_idx,
  x = mtx_data_lung$count,
  dims = c(n_cells, n_genes),
  dimnames = list(barcodes, genes)
)

geneID
gene_translation %>%  filter(SYMBOL %in% geneID) 
raw_counts_matrix <-  raw_counts_matrix[,esid_keep$ENSEMBL] %>%  `colnames<-`(esid_keep$SYMBOL)
raw_counts_matrix

nromallung <-  CreateSeuratObject(counts = t(raw_counts_matrix),   assay = "RNA",  
                           min.cells = 100,
                           meta.data = cell_metadata)
saveRDS(nromallung, "public/nromallung_T.rds")

nromallung@assays$RNA@data

nromallung$development_stage %>%  table()
 

CD3list <- list(c("CD3D", "CD3E", "CD3G"))

TRDlist <- list(c("TRDV1", "TRDV2", "TRDV3", "TRDV4", 'TRDC'))

TRABlist <- list(
  grep("^TRAV|^TRBV|TRBC1|TRBC2|^TRAC", rownames(nromallung), value = T)
  
)
TRABlist

nromallung@meta.data
nromallung  %<>% AddModuleScore(features = CD3list, name = "CD3_score") %>% 
  AddModuleScore(features = TRDlist, name = "TRD_score") %>% 
  AddModuleScore(features = TRABlist, name = "TRAB_score")

ViolinPlot(nromallung, group.by =   "sex", g= c("CD3_score1", "TRD_score1","TRAB_score1"),
           x.angle = 90, colors = ggplotColours(60), sz = 0.5, ncol = 3  ) 



ViolinPlot(nromallung, group.by =   "cell_type", g= c( "TRD_score1"), mythe = F, othertheme = list(theme(axis.text.x = element_text(size = 11, angle = 90))),
           x.angle = 90, sz =1, ncol = 1 , colors = umap.colors ) 
nromallung$cell_type


Feature_rast(nromallung, d1 = "TRAB_score1", d2 = "TRD_score1", 
            othertheme = list( geom_hline(yintercept = 0.20
                                          , linewidth = 0.5, linetype = "dashed"),
                               geom_vline(xintercept = -0.018, linewidth = 0.5,linetype = "dashed")),
            noaxis = F, axis.number = T, do.label = F, g = c(  "TRDV2", "TRDV1"),  navalue = "lightgrey") 



nromallung@meta.data   %<>% mutate(gdTcells_infered = case_when( 
  TRAB_score1 <= 0 & TRD_score1 >=  0.2 ~ "gdT", TRUE ~ "nongdT"
))


Feature_rast(nromallung, d1 = "TRAB_score1", d2 = "TRD_score1", 
             othertheme = list( geom_hline(yintercept = 0.2, linewidth = 0.5, linetype = "dashed"),
                                geom_vline(xintercept = 0, linewidth = 0.5,linetype = "dashed")), mythe =F, colorset = c("red", "lightgrey"),
             noaxis = F, axis.number = T, do.label = F, g = c(   "gdTcells_infered"),  navalue = "lightgrey") 

nromallung@meta.data
nromallung@assays$RNA@data
nromallung_gdT <- subset(nromallung,   gdTcells_infered == "gdT" & CD4 < 0.2) 


dim(nromallung_gdT)

nromallung_gdT@assays$RNA@scale.data

nromallung_gdT@meta.data
nromallung_gdT <- FindVariableFeatures(nromallung_gdT, selection.method = "vst", nfeatures = 1500)


var_genes <- VariableFeatures(nromallung_gdT) %>% 
  str_subset("^(MT|RP[SL]|HIST|^AC|^AL|^AF|-AS1$|^AP|^TRA|^TRB|^IG|LINC|LOC|^MIR|$DT)", negate = TRUE)  %T>% print()
var_genes <-  c(var_genes, c('CD4', "CD8A", "CD8B", "ZNF683", "ITGAE", "CXCR6", "AREG", "CSF1", "CSF2", "GZMB", "GZMA", "RORC"))

nromallung_gdT %<>% ScaleData(vars.to.regress = c("nFeature_RNA", "institute", "donor_age","self_reported_ethnicity"), verbose = FALSE)
nromallung_gdT@meta.data 

nromallung_gdT <- RunPCA(nromallung_gdT, features = var_genes, verbose = FALSE, npcs = 100)
ElbowPlot(nromallung_gdT, ndims = 50) 

nromallung_gdT <- JackStraw(nromallung_gdT, num.replicate = 100, dims = 80)
# ScoreJackStraw(sepsis_2, dims = 1:100) %>% 
#   ggplot(aes(x = dims, y = score)) + 
#   geom_point() + 
#   geom_line() + 
#   labs(title = "Jackstraw scores for Sepsis Batch2 PCA") +
#   theme_minimal()
nromallung_gdT <- ScoreJackStraw(nromallung_gdT, dims = 1:80) 
JackStrawPlot(nromallung_gdT, 1:80) + 
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") +
  labs(title = "Jackstraw plot for gdT lung")


nromallung_gdT <- RunUMAP(nromallung_gdT, dims = 1:19, verbose = FALSE, reduction.key = "UMAP_", min.dist = 0.00) %>% 
  FindNeighbors(dims = 1:19, verbose = FALSE) %>% 
  FindClusters(resolution = 1, verbose = FALSE)

nromallung_gdT <-RunUMAP(nromallung_gdT, dims = 1:20, verbose = FALSE, reduction.key = "UMAP_", min.dist = 0.001)

nromallung_gdT  %<>%    FindClusters(resolution = 1.5, verbose = FALSE)


ViolinPlot(nromallung_gdT, c("TRDV1", "TRDV2", "TRDV3"), group.by = "gdTCR",colors = umap.colors)
nromallung_gdT$bc <-  colnames(nromallung_gdT)

Feature_rast(nromallung_gdT, d1 = "TRDV3", d2 = "TRDV2", g = "gdTCR", noaxis = F, axis.number = T)
nromallung_gdT
Vd1bc <-  nromallung_gdT %>% subset(TRDV1 > 0.5) %>% colnames()
Vd2bc <-  nromallung_gdT %>% subset(TRDV2 > 0.5) %>% colnames()
Vd3bc <-  nromallung_gdT %>% subset(TRDV3 > 0.5) %>% colnames()
nromallung_gdT@meta.data  %<>% mutate(
  TRDVs = case_when(
                    (bc %in% Vd1bc) & ( bc %in% Vd2bc) ~ "DP",
                   bc %in% Vd1bc ~ "Vd1",
                    bc %in% Vd2bc ~ "Vd2",
                    bc %in% Vd3bc ~ "Vd3")
  
)



nromallung_gdT@meta.data  %<>% mutate(gdTCR = case_when(seurat_clusters == "7" ~ "Vd2", TRUE~  "Vd1/Vd3"))
nromallung_gdT$donor_age %>% unique()
Feature_rast(nromallung_gdT, c("ident", "ZNF683", "gdTCR"))
Feature_rast(nromallung_gdT, c("ident", "donor_age"))
Feature_rast(nromallung_gdT, c("ident", "ITGA1", "TRDV2"))

ViolinPlot(nromallung_gdT, "ITGA1", group.by = "seurat_clusters", sz = 1.2)



nromallung_gdT$seurat_clusters
Feature_density(nromallung_gdT, c("TRDV1", "TRDV2", "TRDV3", "ZNF683", "ITGA1", "ITGAE", "TCF7", "CCR7", "KLRB1", "CCR6","DPP4" , "RORC", 
                                  "KLRG1", "AREG", "CSF2", "TBX21", "TNF", "GZMB", "EOMES", "GATA3"), sz = 0.4)

nromallung_gdT@meta.data %>%  group_by(donor_age, seurat_clusters) %>% 
  summarise(n = n() )%>% group_by(donor_age) %>%  mutate(percent = n/sum(n)*100, total = sum(n)) %>%  filter(seurat_clusters == 8) 


nromallung_gdT@meta.data %<>% mutate(age_stage = case_when(  donor_age %in% c("20 years", "20-24 years", "23 years", "25 years", "25-30 years", "27 years","33 years") ~ "young adult",
   donor_age %in% c("35-39 years", "48 years", "40-44 years", "50-54 years") ~ "mid age",
   donor_age %in% c("64 years",  "68 years",    "75 years",   "73 years"  , "70-74 years") ~ "old"
  
  ))

nromallung$donor_age %>% table() %>% prop.table()*100

nromallung_gdT$donor_age %>% table() %>% prop.table()*100
table(nromallung_gdT$donor_age, nromallung_gdT$donor_id )

nromallung_gdT$donor_smoking_status

Feature_rast(nromallung_gdT, "age_stage", facets = "donor_smoking_status")

table(nromallung_gdT$donor_id, nromallung_gdT$age_stage)
saveRDS(nromallung_gdT, "public/nromallung_gdT.rds")

Feature_rast(nromallung_gdT, "TRDVs", noaxis = F, sz = 1.5,
             othertheme = list(coord_fixed()))

Feature_rast(nromallung_gdT, c("ident", "TRDVs"), noaxis = F, sz = 1.5,
             othertheme = list(coord_fixed()))

Feature_density(nromallung_gdT, reduction = "ref.umap",
                othertheme = list(coord_fixed(), NoLegend()),sz = 0.7,
                c("CCR6", "RORC", "IL23R", "KLRB1"))


Feature_density(nromallung_gdT, 
                othertheme = list(coord_fixed(), NoLegend()),sz = 0.7,
                c("AREG", "CSF1", "CSF2", "CTLA4"))


Feature_rast(nromallung_gdT, 
                othertheme = list(coord_fixed(), NoLegend()),sz = 0.7,
                c("ZNF683", "ITGA1", "ITGAE", "CXCR6"))



# intergrate normalllung_gdT to gdT ref

GDTlung_s[["RNA"]]

GDTlung_s_ref@assays
GDTlung_s_ref <- GDTlung_s
GDTlung_s_ref[[c("GM")]] <-  NULL  
GDTlung_s_ref[[c("CITE")]] <-  NULL  
GDTlung_s_ref[[c("HTO")]] <-  NULL  

GDTlung_s_ref[[c("AUC")]] <-  NULL  

GDTlung_s_ref  %<>%  RunUMAP( dims = 1:35, seed=1955,umap.method = "uwot", 
                              return.model = T,assay = "integrated",
                              reduction = 'pca', min.dist = 0.2 ,reduction.key = "umap_" )
Feature_rast(GDTlung_s_ref, "pheno", d1 = "umap_1", d2 = "umap_2")


multicores()
gdTlung.anchors <- FindTransferAnchors(reference = GDTlung_s_ref, query = nromallung_gdT, dims = 1:30,
                                        reference.reduction = "pca")


predictions <- TransferData(anchorset = gdTlung.anchors, refdata = GDTlung_s_ref$pheno, dims = 1:30)

predictions

table(predictions$predicted.id)
nromallung_gdT <- AddMetaData(nromallung_gdT, metadata = predictions)
GDTlung_s_ref@reductions$umap

nromallung_gdT <- MapQuery(anchorset = gdTlung.anchors, reference = GDTlung_s_ref, query = nromallung_gdT,
                           refdata = list(celltype = "pheno"), reference.reduction = "pca", reduction.model = "umap")

nromallung_gdT@reductions$ref.umap
Feature_rast(nromallung_gdT, "predicted.id", d1 = "refUMAP_1", d2 = "refUMAP_2")
Feature_rast(nromallung_gdT, "predicted.id")
umapss <-  FetchData(nromallung_gdT, c("refUMAP_1", "refUMAP_2"))

nromallung_gdT$UMAP_1 <- umapss$refUMAP_1
nromallung_gdT$UMAP_2 <-  -(umapss$refUMAP_2)
nromallung_gdT@reductions$ref.umap@cell.embeddings[,'refUMAP_2']  <- 

-(nromallung_gdT@reductions$ref.umap@cell.embeddings[,'refUMAP_2'])

Feature_density(nromallung_gdT, "ZNF683" ,reduction = "ref.umap" )

ViolinPlot(nromallung_gdT, c("ZNF683", "ITGA1", "ITGAE", "CXCR6"), group.by = "predicted.id", colors = umap.colors)

ClusterCompare(nromallung_gdT, group.by = "predicted.id", 
                 id1 = "Lung_TRM", id2= "Lung_circ")

# numberts 
nromallung_gdT$CD8.Cytotoxictiy1
nromallung_gdT$Tissue.resident1
#caculate scores
for (i in names(sigtable)) {
  nromallung_gdT <- AddModuleScore(nromallung_gdT, features = list(sigtable[[i]]), name = i, assay = 'RNA')
  
}

ViolinPlot(nromallung_gdT, c("CD8.Cytotoxictiy1", "Tissue.resident1"), colors = umap.colors, box = T)
j
TRM_nm <- nromallung_gdT@meta.data %>%  group_by(donor_age, seurat_clusters) %>% 
  summarise(n = n() )%>% group_by(donor_age) %>%  mutate(percent = n/sum(n)*100, total = sum(n)) %>%  filter(seurat_clusters == 8) %>% filter(total > 25)



TRM_nm2 <- nromallung_gdT@meta.data %>%  group_by(donor_age, predicted.id) %>% 
  summarise(n = n() )%>% group_by(donor_age) %>%  mutate(percent = n/sum(n)*100, total = sum(n)) %>%  filter(predicted.id == "Lung_TRM") 


TRM_nm
TRM_nm2

TRM_nm$disease <- "healthy"
TRM_nm2$disease <- "healthy"

TRM_rm_c  <- data.frame(disease = "COPD",
                        percent = c(
                          65.51899, 45.04182, 11.11111, 70.89305, 60.33479, 29.16667, 10.37118
                        ))

TRM_rm_cn <-   data.frame(disease = c(TRM_nm2$disease, TRM_rm_c$disease), percent =c(TRM_nm2$percent, TRM_rm_c$percent)  )
TRM_rm_cn 
ggplot(TRM_rm_cn, aes(x = disease, y = percent)) + geom_boxplot()+geom_point(size = 2)+  ggpubr::stat_compare_means(  paired = F, method = 'wilcox') + theme_classic()+ylab("gdTRM%")

  




