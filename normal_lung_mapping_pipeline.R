#!/usr/bin/env Rscript

# Normal lung γδT mapping workflow (script version)
# --------------------------------------------------
# This script reorganizes the original notebook-style code into clearly
# delineated sections so the overall mapping logic is easier to follow.  The
# analysis covers two public datasets (nromallung_gdT and gdt_integrated) that
# are both mapped to the internal GDTlung reference.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(magrittr)
  library(tibble)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(openxlsx)
})

source("funcs.r")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
config <- list(
  public_files = list(
    h5ad = "normallungT.h5ad",
    mtx = "normallung.mtx"
  ),
  marker_sets = list(
    CD3 = c("CD3D", "CD3E", "CD3G"),
    TRD = c("TRDV1", "TRDV2", "TRDV3", "TRDV4", "TRDC"),
    TRAB_pattern = "^TRAV|^TRBV|TRBC1|TRBC2|^TRAC"
  ),
  signature_xlsx = "abt/abd5778_Table_S3.xlsx",
  reference = list(
    input_rds = "GDTlung_6p_Seurat_07062021.rds",
    output_rds = "GDTlung.ref.rds",
    dims = 1:35,
    min_dist = 0.2,
    seed = 1955
  ),
  query_reference_rds = "Lung_gdT_subset_seurat_ncRNA_rm.rds",
  outputs = list(
    normal_lung_rds = "public/nromallung_T.rds",
    normal_lung_gdt_rds = "public/nromallung_gdT.rds"
  )
)

# -----------------------------------------------------------------------------
# Section 1: Load the public normal lung dataset
# -----------------------------------------------------------------------------
stopifnot(file.exists(config$public_files$h5ad), file.exists(config$public_files$mtx))

mtx_data_lung <- data.table::fread(
  config$public_files$mtx,
  skip = 3,
  col.names = c("i_row_idx", "j_col_idx", "count")
)
dims <- as.numeric(strsplit(readLines(config$public_files$mtx, n = 3)[3], " ")[[1]])

ad <- read_h5ad(config$public_files$h5ad)
genes <- ad$var_names
barcodes <- ad$obs_names
cell_metadata <- ad$obs
gene_metadata <- ad$var

message("Translating gene identifiers ...")
gene_translation <- bitr(
  geneID = genes,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  filter(n() == 1) %>%
  ungroup()

raw_counts_matrix <- sparseMatrix(
  i = mtx_data_lung$i_row_idx,
  j = mtx_data_lung$j_col_idx,
  x = mtx_data_lung$count,
  dims = c(dims[1], dims[2]),
  dimnames = list(barcodes, genes)
)
raw_counts_matrix <- raw_counts_matrix[, gene_translation$ENSEMBL]
colnames(raw_counts_matrix) <- gene_translation$SYMBOL

# Parameters chosen during exploratory analysis
gdT_preprocessing <- list(
  min_cells = 100,
  nfeatures = 1500,
  dims_to_use = 20,
  resolution = 1.5,
  cd4_cutoff = 0.2,
  trd_threshold = 0.2,
  trab_threshold = 0
)

nromallung <- CreateSeuratObject(
  counts = t(raw_counts_matrix),
  assay = "RNA",
  min.cells = gdT_preprocessing$min_cells,
  meta.data = cell_metadata
)

if (!is.null(config$outputs$normal_lung_rds)) {
  saveRDS(nromallung, config$outputs$normal_lung_rds)
}

# -----------------------------------------------------------------------------
# Section 2: Score γδT-associated modules and infer gdT cells
# -----------------------------------------------------------------------------
message("Scoring γδT modules ...")
trab_genes <- grep(config$marker_sets$TRAB_pattern, rownames(nromallung), value = TRUE)

nromallung <- nromallung %>%
  AddModuleScore(features = list(config$marker_sets$CD3), name = "CD3_score") %>%
  AddModuleScore(features = list(config$marker_sets$TRD), name = "TRD_score") %>%
  AddModuleScore(features = list(trab_genes), name = "TRAB_score")

nromallung$gdTcells_infered <- if_else(
  nromallung$TRAB_score1 <= gdT_preprocessing$trab_threshold &
    nromallung$TRD_score1 >= gdT_preprocessing$trd_threshold,
  "gdT",
  "nongdT"
)

# -----------------------------------------------------------------------------
# Section 3: Preprocess inferred gdT population
# -----------------------------------------------------------------------------
message("Preprocessing inferred γδT subset ...")

nromallung_gdT <- subset(
  nromallung,
  gdTcells_infered == "gdT" & CD4 < gdT_preprocessing$cd4_cutoff
)

nromallung_gdT <- FindVariableFeatures(
  nromallung_gdT,
  selection.method = "vst",
  nfeatures = gdT_preprocessing$nfeatures
)
var_genes <- VariableFeatures(nromallung_gdT) %>%
  str_subset("^(MT|RP[SL]|HIST|^AC|^AL|^AF|-AS1$|^AP|^TRA|^TRB|^IG|LINC|LOC|^MIR|$DT)", negate = TRUE) %>%
  union(c("CD4", "CD8A", "CD8B", "ZNF683", "ITGAE", "CXCR6", "AREG", "CSF1", "CSF2", "GZMB", "GZMA", "RORC"))

nromallung_gdT <- ScaleData(
  nromallung_gdT,
  vars.to.regress = c("nFeature_RNA", "institute", "donor_age", "self_reported_ethnicity"),
  verbose = FALSE
)

nromallung_gdT <- RunPCA(nromallung_gdT, features = var_genes, npcs = 100, verbose = FALSE)
nromallung_gdT <- RunUMAP(
  nromallung_gdT,
  dims = 1:gdT_preprocessing$dims_to_use,
  reduction.key = "UMAP_",
  min.dist = 0.001,
  verbose = FALSE
)
nromallung_gdT <- FindNeighbors(nromallung_gdT, dims = 1:gdT_preprocessing$dims_to_use, verbose = FALSE)
nromallung_gdT <- FindClusters(nromallung_gdT, resolution = gdT_preprocessing$resolution, verbose = FALSE)

nromallung_gdT$bc <- colnames(nromallung_gdT)
vd1_cells <- WhichCells(nromallung_gdT, expression = TRDV1 > 0.5)
vd2_cells <- WhichCells(nromallung_gdT, expression = TRDV2 > 0.5)
vd3_cells <- WhichCells(nromallung_gdT, expression = TRDV3 > 0.5)

nromallung_gdT$TRDVs <- case_when(
  nromallung_gdT$bc %in% vd1_cells & nromallung_gdT$bc %in% vd2_cells ~ "DP",
  nromallung_gdT$bc %in% vd1_cells ~ "Vd1",
  nromallung_gdT$bc %in% vd2_cells ~ "Vd2",
  nromallung_gdT$bc %in% vd3_cells ~ "Vd3",
  TRUE ~ NA_character_
)

nromallung_gdT$gdTCR <- if_else(nromallung_gdT$seurat_clusters == "7", "Vd2", "Vd1/Vd3")

nromallung_gdT$age_stage <- case_when(
  nromallung_gdT$donor_age %in% c("20 years", "20-24 years", "23 years", "25 years", "25-30 years", "27 years", "33 years") ~ "young adult",
  nromallung_gdT$donor_age %in% c("35-39 years", "48 years", "40-44 years", "50-54 years") ~ "mid age",
  nromallung_gdT$donor_age %in% c("64 years", "68 years", "75 years", "73 years", "70-74 years") ~ "old",
  TRUE ~ NA_character_
)

# -----------------------------------------------------------------------------
# Section 4: Score published signature modules
# -----------------------------------------------------------------------------
message("Scoring published signature modules ...")
stopifnot(file.exists(config$signature_xlsx))

signature_table <- openxlsx::read.xlsx(config$signature_xlsx)
colnames(signature_table) <- str_remove(colnames(signature_table), "\.\\(.+\\)")
signature_list <- as.list(signature_table) %>%
  map(~ .x[!is.na(.x)])
signature_list$Tissue.resident <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "RBPJ", "GPR25")

for (signature_name in names(signature_list)) {
  nromallung_gdT <- AddModuleScore(
    nromallung_gdT,
    features = list(signature_list[[signature_name]]),
    name = signature_name,
    assay = "RNA"
  )
}
colnames(nromallung_gdT@meta.data) <- str_replace(colnames(nromallung_gdT@meta.data), "(?<=\\w)1$", "")

if (!is.null(config$outputs$normal_lung_gdt_rds)) {
  saveRDS(nromallung_gdT, config$outputs$normal_lung_gdt_rds)
}

# -----------------------------------------------------------------------------
# Section 5: Prepare the GDTlung reference
# -----------------------------------------------------------------------------
message("Preparing GDTlung reference ...")
GDTlung_s_ref <- readRDS(config$reference$input_rds)

assays_to_drop <- intersect(c("GM", "CITE", "HTO", "AUC"), names(GDTlung_s_ref@assays))
if (length(assays_to_drop) > 0) {
  for (assay_name in assays_to_drop) {
    GDTlung_s_ref[[assay_name]] <- NULL
  }
}

GDTlung_s_ref <- RunUMAP(
  GDTlung_s_ref,
  dims = config$reference$dims,
  seed.use = config$reference$seed,
  umap.method = "uwot",
  return.model = TRUE,
  assay = "integrated",
  reduction = "pca",
  min.dist = config$reference$min_dist,
  reduction.key = "umap_"
)

ref_embeddings <- Embeddings(GDTlung_s_ref, "umap")
ref_embeddings[, "umap_2"] <- -ref_embeddings[, "umap_2"]
GDTlung_s_ref@reductions$umap@cell.embeddings <- ref_embeddings
GDTlung_s_ref@reductions$umap@misc$model$embedding[, 2] <- -GDTlung_s_ref@reductions$umap@misc$model$embedding[, 2]
GDTlung_s_ref$Umap_1 <- ref_embeddings[, "umap_1"]
GDTlung_s_ref$Umap_2 <- ref_embeddings[, "umap_2"]

if (!is.null(config$reference$output_rds)) {
  saveRDS(GDTlung_s_ref, config$reference$output_rds)
}

# -----------------------------------------------------------------------------
# Section 6: Map public datasets to the reference
# -----------------------------------------------------------------------------
message("Mapping nromallung_gdT to the reference ...")

normal_anchors <- FindTransferAnchors(
  reference = GDTlung_s_ref,
  query = nromallung_gdT,
  dims = 1:30,
  reference.reduction = "pca"
)
normal_predictions <- TransferData(
  anchorset = normal_anchors,
  refdata = GDTlung_s_ref$pheno,
  dims = 1:30
)

nromallung_gdT <- AddMetaData(nromallung_gdT, metadata = normal_predictions)
nromallung_gdT <- MapQuery(
  anchorset = normal_anchors,
  reference = GDTlung_s_ref,
  query = nromallung_gdT,
  refdata = list(celltype = "pheno"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
ref_umap_embeddings <- Embeddings(nromallung_gdT, "ref.umap")
ref_umap_embeddings[, "refUMAP_2"] <- -ref_umap_embeddings[, "refUMAP_2"]
nromallung_gdT@reductions$ref.umap@cell.embeddings <- ref_umap_embeddings
nromallung_gdT$Umap_1 <- ref_umap_embeddings[, "refUMAP_1"]
nromallung_gdT$Umap_2 <- ref_umap_embeddings[, "refUMAP_2"]
nromallung_gdT$pheno <- nromallung_gdT$predicted.id

message("Mapping gdt_integrated dataset to the reference ...")
gdt_integrated <- readRDS(config$query_reference_rds)

gdt_anchors <- FindTransferAnchors(
  reference = GDTlung_s_ref,
  query = gdt_integrated,
  dims = 1:30,
  reference.reduction = "pca"
)
gdt_predictions <- TransferData(
  anchorset = gdt_anchors,
  refdata = GDTlung_s_ref$pheno,
  dims = 1:30
)

gdt_integrated <- AddMetaData(gdt_integrated, metadata = gdt_predictions)
gdt_integrated <- MapQuery(
  anchorset = gdt_anchors,
  reference = GDTlung_s_ref,
  query = gdt_integrated,
  refdata = list(celltype = "pheno"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
ref_umap_integrated <- Embeddings(gdt_integrated, "ref.umap")
ref_umap_integrated[, "refUMAP_2"] <- -ref_umap_integrated[, "refUMAP_2"]
gdt_integrated@reductions$ref.umap@cell.embeddings <- ref_umap_integrated
gdt_integrated$Umap_1 <- ref_umap_integrated[, "refUMAP_1"]
gdt_integrated$Umap_2 <- ref_umap_integrated[, "refUMAP_2"]

gdTlung.anchors <- list(nromallung = normal_anchors, gdt_integrated = gdt_anchors)

# -----------------------------------------------------------------------------
# Section 7: Summaries and comparisons
# -----------------------------------------------------------------------------
message("Summarizing gdTRM proportions ...")

cluster_summary <- nromallung_gdT@meta.data %>%
  group_by(donor_age, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percent = n / sum(n) * 100, total = sum(n)) %>%
  ungroup()
trm_cluster_summary <- cluster_summary %>%
  filter(seurat_clusters == 8, total > 25)

prediction_summary <- nromallung_gdT@meta.data %>%
  group_by(donor_age, predicted.id) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percent = n / sum(n) * 100, total = sum(n)) %>%
  ungroup()
trm_prediction_summary <- prediction_summary %>%
  filter(predicted.id == "Lung_TRM")

healthy_trm <- trm_prediction_summary %>%
  mutate(disease = "healthy") %>%
  select(disease, percent)

copd_percentages <- c(65.51899, 45.04182, 11.11111, 70.89305, 60.33479, 29.16667, 10.37118)

copd_trm <- tibble(
  disease = "COPD",
  percent = copd_percentages
)
trm_disease_comparison <- bind_rows(healthy_trm, copd_trm)

normal_lung_mapping_results <- list(
  normal_lung = nromallung,
  normal_lung_gdT = nromallung_gdT,
  gdt_reference = GDTlung_s_ref,
  gdt_integrated = gdt_integrated,
  anchors = gdTlung.anchors,
  summaries = list(
    cluster = trm_cluster_summary,
    predicted = trm_prediction_summary,
    disease = trm_disease_comparison
  )
)

message("normal_lung_mapping_pipeline.R completed. Access 'normal_lung_mapping_results' for outputs.")
