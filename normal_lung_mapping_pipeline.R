#!/usr/bin/env Rscript

# Normal lung γδT mapping workflow
# ---------------------------------
# This script reorganizes the ad-hoc notebook-style code that was previously
# embedded in IL32_PBMC_coculture.R into a reusable set of helpers.  The goal is
# to make it easier to (1) load the public normal lung dataset, (2) score γδT
# markers, (3) preprocess the inferred γδT population, (4) prepare the
# GDTlung reference, and (5) map/query public datasets to that reference for
# downstream analyses.

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
})

source("funcs.r")

# ---------------------------------------------------------------------------
# Configuration -----------------------------------------------------------------
# ---------------------------------------------------------------------------

default_normal_lung_config <- list(
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
  ),
  gdT_preprocessing = list(
    min_cells = 100,
    nfeatures = 1500,
    dims_to_use = 20,
    resolution = 1.5,
    cd4_cutoff = 0.2,
    trd_threshold = 0.2,
    trab_threshold = 0
  ),
  copd_percentages = c(65.51899, 45.04182, 11.11111, 70.89305, 60.33479, 29.16667, 10.37118)
)

# ---------------------------------------------------------------------------
# Helper functions ---------------------------------------------------------
# ---------------------------------------------------------------------------

load_public_lung_dataset <- function(paths, min_cells = 100) {
  stopifnot(file.exists(paths$h5ad), file.exists(paths$mtx))

  mtx_data <- data.table::fread(
    paths$mtx,
    skip = 3,
    col.names = c("i_row_idx", "j_col_idx", "count")
  )
  dims <- as.numeric(strsplit(readLines(paths$mtx, n = 3)[3], " ")[[1]])

  ad <- read_h5ad(paths$h5ad)
  genes <- ad$var_names
  barcodes <- ad$obs_names
  cell_metadata <- ad$obs

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
    i = mtx_data$i_row_idx,
    j = mtx_data$j_col_idx,
    x = mtx_data$count,
    dims = c(dims[1], dims[2]),
    dimnames = list(barcodes, genes)
  )

  raw_counts_matrix <- raw_counts_matrix[, gene_translation$ENSEMBL]
  colnames(raw_counts_matrix) <- gene_translation$SYMBOL

  CreateSeuratObject(
    counts = t(raw_counts_matrix),
    assay = "RNA",
    min.cells = min_cells,
    meta.data = cell_metadata
  )
}

score_gdt_modules <- function(obj, marker_sets) {
  trab_genes <- grep(marker_sets$TRAB_pattern, rownames(obj), value = TRUE)

  obj %>%
    AddModuleScore(features = list(marker_sets$CD3), name = "CD3_score") %>%
    AddModuleScore(features = list(marker_sets$TRD), name = "TRD_score") %>%
    AddModuleScore(features = list(trab_genes), name = "TRAB_score")
}

infer_gdt_population <- function(obj, trab_threshold = 0, trd_threshold = 0.2) {
  obj$gdTcells_infered <- if_else(
    obj$TRAB_score1 <= trab_threshold & obj$TRD_score1 >= trd_threshold,
    "gdT",
    "nongdT"
  )
  obj
}

subset_gdt_population <- function(obj, cd4_cutoff = 0.2) {
  subset(obj, gdTcells_infered == "gdT" & CD4 < cd4_cutoff)
}

get_variable_genes <- function(obj) {
  VariableFeatures(obj) %>%
    str_subset("^(MT|RP[SL]|HIST|^AC|^AL|^AF|-AS1$|^AP|^TRA|^TRB|^IG|LINC|LOC|^MIR|$DT)", negate = TRUE) %>%
    union(c(
      "CD4", "CD8A", "CD8B", "ZNF683", "ITGAE", "CXCR6", "AREG", "CSF1", "CSF2", "GZMB", "GZMA", "RORC"
    ))
}

preprocess_gdt_population <- function(
    obj,
    nfeatures = 1500,
    dims_to_use = 20,
    resolution = 1.5,
    regressors = c("nFeature_RNA", "institute", "donor_age", "self_reported_ethnicity")
) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  var_genes <- get_variable_genes(obj)

  obj <- ScaleData(
    obj,
    vars.to.regress = regressors,
    verbose = FALSE
  )

  obj <- RunPCA(obj, features = var_genes, npcs = 100, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:dims_to_use, reduction.key = "UMAP_", min.dist = 0.001, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dims_to_use, verbose = FALSE)
  FindClusters(obj, resolution = resolution, verbose = FALSE)
}

assign_trdv_usage <- function(obj, threshold = 0.5) {
  obj$bc <- colnames(obj)

  vd1_cells <- WhichCells(obj, expression = TRDV1 > threshold)
  vd2_cells <- WhichCells(obj, expression = TRDV2 > threshold)
  vd3_cells <- WhichCells(obj, expression = TRDV3 > threshold)

  obj$TRDVs <- case_when(
    obj$bc %in% vd1_cells & obj$bc %in% vd2_cells ~ "DP",
    obj$bc %in% vd1_cells ~ "Vd1",
    obj$bc %in% vd2_cells ~ "Vd2",
    obj$bc %in% vd3_cells ~ "Vd3",
    TRUE ~ NA_character_
  )

  obj$gdTCR <- if_else(obj$seurat_clusters == "7", "Vd2", "Vd1/Vd3")
  obj
}

prepare_gdt_reference <- function(
    ref_obj,
    dims = 1:35,
    min_dist = 0.2,
    seed = 1955,
    assay = "integrated",
    reduction = "pca"
) {
  assays_to_drop <- intersect(c("GM", "CITE", "HTO", "AUC"), names(ref_obj@assays))
  if (length(assays_to_drop) > 0) {
    for (assay_name in assays_to_drop) {
      ref_obj[[assay_name]] <- NULL
    }
  }

  ref_obj <- RunUMAP(
    ref_obj,
    dims = dims,
    seed.use = seed,
    umap.method = "uwot",
    return.model = TRUE,
    assay = assay,
    reduction = reduction,
    min.dist = min_dist,
    reduction.key = "umap_"
  )

  ref_embeddings <- Embeddings(ref_obj, "umap")
  ref_embeddings[, "umap_2"] <- -ref_embeddings[, "umap_2"]
  ref_obj@reductions$umap@cell.embeddings <- ref_embeddings
  ref_obj@reductions$umap@misc$model$embedding[, 2] <- -ref_obj@reductions$umap@misc$model$embedding[, 2]

  ref_obj$Umap_1 <- ref_embeddings[, "umap_1"]
  ref_obj$Umap_2 <- ref_embeddings[, "umap_2"]
  ref_obj
}

map_query_to_reference <- function(
    query_obj,
    reference_obj,
    dims = 1:30,
    label = "pheno",
    reference_reduction = "pca",
    reduction_model = "umap"
) {
  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query = query_obj,
    dims = dims,
    reference.reduction = reference_reduction
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_obj[[label]],
    dims = dims
  )

  query_obj <- AddMetaData(query_obj, metadata = predictions)
  query_obj <- MapQuery(
    anchorset = anchors,
    reference = reference_obj,
    query = query_obj,
    refdata = list(celltype = label),
    reference.reduction = reference_reduction,
    reduction.model = reduction_model
  )

  ref_embeddings <- Embeddings(query_obj, "ref.umap")
  ref_embeddings[, "refUMAP_2"] <- -ref_embeddings[, "refUMAP_2"]
  query_obj@reductions$ref.umap@cell.embeddings <- ref_embeddings
  query_obj$Umap_1 <- ref_embeddings[, "refUMAP_1"]
  query_obj$Umap_2 <- ref_embeddings[, "refUMAP_2"]

  list(query = query_obj, anchors = anchors, predictions = predictions)
}

score_signature_modules <- function(obj, signature_list, assay = "RNA") {
  for (sig_name in names(signature_list)) {
    obj <- AddModuleScore(
      obj,
      features = list(signature_list[[sig_name]]),
      name = sig_name,
      assay = assay
    )
  }

  colnames(obj@meta.data) <- str_replace(colnames(obj@meta.data), "(?<=\\w)1$", "")
  obj
}

compute_percent_by_group <- function(meta, grouping, category) {
  meta %>%
    group_by(across(all_of(c(grouping, category)))) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(percent = n / sum(n) * 100, total = sum(n)) %>%
    ungroup()
}

annotate_age_stage <- function(meta) {
  if (!"age_stage" %in% colnames(meta)) {
    meta$age_stage <- NA_character_
  }

  meta %>%
    mutate(
      age_stage = case_when(
        donor_age %in% c("20 years", "20-24 years", "23 years", "25 years", "25-30 years", "27 years", "33 years") ~ "young adult",
        donor_age %in% c("35-39 years", "48 years", "40-44 years", "50-54 years") ~ "mid age",
        donor_age %in% c("64 years", "68 years", "75 years", "73 years", "70-74 years") ~ "old",
        TRUE ~ age_stage
      )
    )
}

load_signature_modules <- function(xlsx_path) {
  stopifnot(file.exists(xlsx_path))

  modules <- openxlsx::read.xlsx(xlsx_path) %>%
    `colnames<-`(str_remove(colnames(.), '.\\\(.+\\\)')) %>%
    as.list() %>%
    map(~ na.exclude(.x) %>% as.vector())

  modules$Tissue.resident <- c("ITGAE", "ITGA1", "ZNF683", "CXCR6", "RBPJ", "GPR25")
  modules
}

summarize_gdTRM_proportions <- function(normal_meta, copd_percentages) {
  trm_percent <- compute_percent_by_group(normal_meta, "donor_age", "predicted.id") %>%
    filter(predicted.id == "Lung_TRM") %>%
    mutate(disease = "healthy") %>%
    select(disease, percent)

  copd_df <- tibble(
    disease = "COPD",
    percent = copd_percentages
  )

  bind_rows(trm_percent, copd_df)
}

# ---------------------------------------------------------------------------
# Main orchestration ------------------------------------------------------
# ---------------------------------------------------------------------------

run_normal_lung_mapping <- function(config = default_normal_lung_config) {
  marker_sets <- config$marker_sets
  signatures <- load_signature_modules(config$signature_xlsx)

  message("Loading public normal lung dataset ...")
  nromallung <- load_public_lung_dataset(
    config$public_files,
    min_cells = config$gdT_preprocessing$min_cells
  ) %>%
    score_gdt_modules(marker_sets) %>%
    infer_gdt_population(
      trab_threshold = config$gdT_preprocessing$trab_threshold,
      trd_threshold = config$gdT_preprocessing$trd_threshold
    )

  if (!is.null(config$outputs$normal_lung_rds)) {
    saveRDS(nromallung, config$outputs$normal_lung_rds)
  }

  message("Preprocessing inferred γδT population ...")
  nromallung_gdT <- nromallung %>%
    subset_gdt_population(config$gdT_preprocessing$cd4_cutoff) %>%
    preprocess_gdt_population(
      nfeatures = config$gdT_preprocessing$nfeatures,
      dims_to_use = config$gdT_preprocessing$dims_to_use,
      resolution = config$gdT_preprocessing$resolution
    ) %>%
    assign_trdv_usage()

  nromallung_gdT <- score_signature_modules(nromallung_gdT, signatures)
  nromallung_gdT@meta.data <- annotate_age_stage(nromallung_gdT@meta.data)

  if (!is.null(config$outputs$normal_lung_gdt_rds)) {
    saveRDS(nromallung_gdT, config$outputs$normal_lung_gdt_rds)
  }

  message("Preparing GDT lung reference ...")
  gdt_reference <- readRDS(config$reference$input_rds) %>%
    prepare_gdt_reference(
      dims = config$reference$dims,
      min_dist = config$reference$min_dist,
      seed = config$reference$seed
    )

  if (!is.null(config$reference$output_rds)) {
    saveRDS(gdt_reference, config$reference$output_rds)
  }

  message("Mapping public gdT population to reference ...")
  nromallung_mapping <- map_query_to_reference(
    query_obj = nromallung_gdT,
    reference_obj = gdt_reference
  )
  nromallung_gdT <- nromallung_mapping$query

  gdt_integrated <- readRDS(config$query_reference_rds)
  gdt_mapping <- map_query_to_reference(
    query_obj = gdt_integrated,
    reference_obj = gdt_reference
  )
  gdt_integrated <- gdt_mapping$query

  percent_by_cluster <- compute_percent_by_group(nromallung_gdT@meta.data, "donor_age", "seurat_clusters")
  trm_cluster_summary <- percent_by_cluster %>% filter(seurat_clusters == 8, total > 25)
  trm_prediction_summary <- compute_percent_by_group(nromallung_gdT@meta.data, "donor_age", "predicted.id") %>%
    filter(predicted.id == "Lung_TRM")
  trm_disease_comparison <- summarize_gdTRM_proportions(
    nromallung_gdT@meta.data,
    config$copd_percentages
  )

  list(
    normal_lung = nromallung,
    normal_lung_gdT = nromallung_gdT,
    gdt_reference = gdt_reference,
    gdt_integrated = gdt_integrated,
    mappings = list(
      normal_lung = nromallung_mapping,
      gdt_integrated = gdt_mapping
    ),
    summaries = list(
      cluster = trm_cluster_summary,
      predicted = trm_prediction_summary,
      disease = trm_disease_comparison
    )
  )
}

message(
  "normal_lung_mapping_pipeline.R loaded. Call run_normal_lung_mapping() with a\n",
  "custom config (or default_normal_lung_config) to execute the workflow."
)
