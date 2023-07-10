# the analysis script for gammadelta T cells from lung and luLN
# author: Likai Tan

library(BiocManager)
library(devtools)
library(flowCore)
library(gtools)
library(flowWorkspace)
library(flowStats)
library(metR)
library(FlowSOM)
library(filesstrings)
library(reticulate)
# library(cytofkit2)
library(PeacoQC)
library(ConsensusClusterPlus)
# library(Rtsne)
library(umap)
library(matrixStats)

library(kableExtra)
library(tidyverse)
library(forcats)
library(reshape2)
library(pheatmap)
library(splitstackshape)
library(RColorBrewer)
library(Seurat)
# library('mlr')

source('/home/big/tanlikai/script/rscripts//funcs.r')


# readprocessed dataframe -------------------------------------------------

gdT_corrected <- readRDS('gdTcell_processed_dataframe.rds')




# ##metadata tabl ---------------------------------------------------------



setwd('/home/big/tanlikai/Lung/FACSdata/')


# setwd('/home/big/tanlikai/Lung/PBMCs_Lung_panel_test/lung_LN/')
dir.create('figs')
fcsfiles <- list.files(pattern = c('gdT', 'fcs') )
file_name <- list.files(pattern = c('gdT', 'fcs') )

# meta data table 
md <-  data.frame(
                  #file 
                  file_name = file_name, 
                  #id
                  donor= str_extract(file_name,'p\\d\\d') ,
                  #group
                  # stimulation =  str_extract(file_name,'stimulated|unstimulated')   ,
                  tissue =   str_extract(file_name,'LN|Lung')  ) %>%  
  mutate(ID = paste0(donor, '_',tissue))



md


###input trimmed FCS file
FCS_raw_gdT <- read.flowSet(
  files = md$file_name,
  transformation = FALSE,
  truncate_max_range = FALSE)



for (i in file_name) {
  print(i)
  dim(FCS_raw_gdT@frames[[i]]) %>%  print()
}

# add barcode to every cell
for (i in file_name) {
  print(i)
  FCS_raw_gdT@frames[[i]]
}



# setwd('..')
# colQuantiles(FCS_humangdT@frames[['export_DMHH190325 IELps_trialUnmixedSamplesKO_Live cells.fcs']]@exprs,
#              probs = c(0.01, 0.2, 0.5, 0.75, 0.99))


##panels
panel <- FCS_raw_gdT@frames[[file_name[1]]]@parameters@data[, 1:2 ]
nrow(panel)

# panel_fcs <- parameters(FCS_raw_gdT[[1]]) %>% Biobase::pData()
panel$desc <- gsub("-","",panel$desc)
Channel <- colnames(FCS_raw_gdT)

panel[7:32,] %>% kbl(row.names = F) %>% kable_styling(full_width = F, position="left")


# name	desc
# FJComp-APC-A	Vg9
# FJComp-APC-Fire 810-A	CD3
# FJComp-APC-Vio 770-A	Vd2
# FJComp-Alexa Fluor 700-A	CD127
# FJComp-BUV395-A	CD45RA
# FJComp-BUV496-A	CD16
# FJComp-BUV563-A	CD4
# FJComp-BUV661-A	CD49a
# FJComp-BUV737-A	KLRG1
# FJComp-BV421-A	GATA3
# FJComp-BV480-A	CD103
# FJComp-BV570-A	CD45RO
# FJComp-BV650-A	CD357
# FJComp-BV711-A	CCR6
# FJComp-BV750-A	CD26
# FJComp-BV785-A	CD8
# FJComp-FITC-A	Eomes
# FJComp-PE-A	AREG
# FJComp-PE-Cy5-A	TCRab
# FJComp-PE-Cy7-A	CTLA4
# FJComp-PE-Dazzle594-A	GMCSF
# FJComp-PE-Fire 700-A	CD25
# FJComp-PerCP-Cy5.5-A	GzmA
# FJComp-PerCP-Vio700-A	TCRgd
# FJComp-VioGreen-A	Vd1
# FJComp-Zombie NIR-A	Live


# # Lineage markers
all_markers <- panel$desc[c(7:32)] %>% as.character()
lineage_markers <- setdiff(all_markers, c('Live','CD3', 'CD45RO', 'AREG', 'Gzma', 'GMCSF','TCRab','TCRgd'))
lineage_markers_2 <-c('Vg9', 'Vd2', 'Vd1' ,'CD45RA', 'CD16', 'CD4', 'CD8', 'CD49a','KLRG1', 'GATA3','CCR6',
                      'CD103', 'CD45RO', 'CD357', 'CD26', 'Eomes', 'CTLA4', 'CD25', 'CD127')
# # Functional markers  we don't use it here but should be useful for furture experiment design
functional_markers <-c('Vd1','Vd2','CD3',"TCRgd", 'TCRab', 'AREG', 'Gzma', 'GMCSF', 'TCRgd', 'TCRab')

#abs the expr matrix

# for (i in file_name) {
#   FCS_humangdT@frames[[i]]@exprs <- abs(FCS_humangdT@frames[[i]]@exprs)
# }


# ##biexpontential transform of data --------------------------------------


##############check rawdata distribution
rawdata <- FCS_raw_gdT@frames[[file_name[1]]]@exprs[,7:32] %>% asinh()
rawdata
colnames(rawdata) <- all_markers

channales <- all_markers

raw <- list()
# raw[[11]]
colnm  <- all_markers
for (i in 1:length(channales)) {
  raw[[i]] <- ggplot(rawdata %>% as.data.frame(), aes_string(x = channales[i])) +
    geom_density()+ggtitle(all_markers[i] )
}
plot_grid(plotlist = raw, ncol = 6)

figsave(plot_grid(plotlist = raw, ncol = 6), "rawsignals_Lung_gdT.pdf", w = 300, h = 300)

# save_plot("raw_Lung.pdf", plot_grid(plotlist = raw, ncol = 6), base_height = 20, base_width = 20)

# figsave(plot_grid(plotlist = raw, ncol = 6), 'rawsignals_lung.pdf', 300, 300)


# data transformation 
fcs_gdT_transform_pre_qc <- fsApply(FCS_raw_gdT, function(x) {transform(x, estimateLogicle(x, c(Channel[7:32])))})

# save the transformed data
write.flowSet(fcs_gdT_transform_pre_qc, outdir = file.path(getwd(), "Transformed_FCS_files") , filename = paste0("",FCS_raw_gdT@phenoData@data$name)) 
#  the cell number too low, QCsckip---
fcs_gdT_transform <- read.flowSet(files = md$file_name, path=file.path(getwd(), "Transformed_FCS_files"), transformation = FALSE, truncate_max_range = FALSE)
fcs_gdT_transform <- fcs_gdT_transform[sampleNames(fcs_gdT_transform), colnames(fcs_gdT_transform)[7:32]]



# automatic QC,  the cell number too low, QCsckip------------------------------------------------


# 
# # Automated quality control ,
# # 
# dir.create('PeacoQCresults')
# 
# 
#   
# for(i in 1:length(sampleNames(fcs_gdT_transform_pre_qc))){
#   ff <-fcs_gdT_transform_pre_qc[[i]]
#   channels=Channel[7:32]
#   peacoqc_res <- PeacoQC(ff, Channel[7:32], determine_good_cells = "all",
#                               save_fcs = TRUE, plot=TRUE, output_directory = "PeacoQCresults")
# } 
# 
# 
# fcs_gdT_transform <- read.flowSet(pattern = 'gdT',  path=file.path(getwd(), "PeacoQCresults","PeacoQC_results","fcs_files"), transformation = FALSE, truncate_max_range = F)
# 
# # 
# 
# 
# 
# # QC result
# sample_ids_raw <- rep(file_name, fsApply(fcs_gdT_transform_pre_qc, nrow))
# 
# sample_ids_trans <- rep(file_name, fsApply(fcs_gdT_transform , nrow))
# cell_table_raw <- table(sample_ids_raw)
# cell_table_trans <- table(sample_ids_trans)
# 
# 
# qc_table <- data.frame(cell_table_raw,cell_table_trans)
# qc_table$sample_ids_trans <- NULL
# colnames(qc_table) <- c("Sample_ID","Pre-QC cell count","Post-QC cell count")
# qc_table$Removed <- qc_table$`Pre-QC cell count`-qc_table$`Post-QC cell count`
# qc_table$Removed_perc <- round(((qc_table$Removed/qc_table$`Pre-QC cell count`)*100), digits=2)
# colnames(qc_table) <- c("Sample_ID","Pre-QC cell count","Post-QC cell count","Removed (n)","Removed (%)")
# qc_table  %>% kbl(row.names = T) %>% kable_styling(full_width = F, position="left")
# 
# 
fcs_gdT_transform <- fcs_gdT_transform[sampleNames(fcs_gdT_transform), colnames(fcs_gdT_transform)[7:32]]



### Generate sample IDs corresponding to each cell in the `expr` matrix

sample_ids <- rep(md$file_name, fsApply(fcs_gdT_transform , nrow))

# donor <- rep(md$donor, fsApply(FCS_humangdT, nrow))

# stimulation <- rep(md$stimulation, fsApply(FCS_humangdT, nrow))




# Batch effection correction ----------------------------------------------




library(cyCombine)


write.csv(md, 'Transformed_FCS_files/metadata_gdT_lung_LN.csv')
md

uncorrected  <- prepare_data(data_dir = 'Transformed_FCS_files/',
                          
                            markers = all_markers,
                            transform = FALSE,
                            pattern = c(  "fcs"),      
                            metadata = 'Transformed_FCS_files/metadata_gdT_lung_LN.csv',
                            # Can also be .csv file or data.frame object
                            filename_col  = "file_name",
                            batch_ids = "donor",
                            condition = "tissue",
                            sample_ids = 'file_name',
                                
                            down_sample =F,
                            # sample_size = 500000,
                            # seed = 473,
                            cofactor = 5)



head(uncorrected)
uncorrected$condition

# uncorrected$sample_id
dim(uncorrected)




uncorrected<- uncorrected[,!is.na(colnames(uncorrected))]
as.data.frame(uncorrected)
# colnames(dr_umap)
# uncorrected <- dr_umap[,c( 3:26,29:32)]
dir.create('figs/batch_gd')

all_markers
uncorrected %>%
  detect_batch_effect(markers = all_markers,
                      batch_col = 'batch',
                      out_dir = 'figs/batch_gd', 
                      seed = 434,
                      name = 'gdTcells')

colnames(uncorrected)


map(all_markers, ~  ggplot(uncorrected, aes_string(x = ., color = 'batch',
                                                 group = 'batch')) +
      geom_density() ) %>% PG(nrow = 5) %T>% figsave('batch_gd/batch_uncorrected_histgram.pdf', 600, 400)

uncorrected$condition
corrected <- uncorrected %>%
  batch_correct(markers = all_markers,
                   # out_dir = 'figs/batch',
              xdim = 10, ydim = 10, covar = 'condition',
                norm_method = "scale", # "rank" is recommended when combining data with heavy batch effects
                rlen = 20, seed = 22)

head(corrected)
start_time_8 <- Sys.time()
labels <- corrected %>%
  cyCombine::create_som(rlen = 20,
                        xdim = 10,
                        ydim = 10,
                        seed = 22,
                        markers = gsub('_', '', all_markers))
end_time_8 <- Sys.time()

# Add labels
corrected <- corrected %>%
  dplyr::mutate(som = labels)




map(all_markers, ~  ggplot(corrected, aes_string(x = ., color = 'batch',
                                          group = 'batch')) +
      geom_density() ) %>% PG(nrow = 5) %T>% figsave('batch_gd/batch_corrected_histgram.pdf', 600, 400)


# Set column for evaluation of EMD (per-cluster)
celltype_col <- "som"
colnames(corrected)

saveRDS(corrected, 'gdt_corrected_FACS_data.rds')

# Transfer labels to uncorrected data
uncorrected <- corrected %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(uncorrected, by = "id")


uncorrected$som <- corrected$som

uncorrected$id <- corrected$id
# Evaluation using EMD
emd_val <- uncorrected %>%
  cyCombine::evaluate_emd(corrected,
                          binSize = 0.1,
                          markers = gsub('_', '', all_markers),
                          cell_col = celltype_col)

# Show plots
cowplot::plot_grid(emd_val$violin, emd_val$scatterplot)
plot_density(uncorrected, corrected, ncol = 6)  %T>% figsave('batch_gd/histogram_corrected_vs_uncorrected.pdf', 400, 400)

plot1 <- plot_dimred(uncorrected, name = 'Uncorrected', type = 'umap', markers = lineage_markers)
plot2 <- plot_dimred(corrected, name = 'Corrected 8x8', type = 'umap', markers = lineage_markers)




plot3 <- plot_dimred(corrected, name = 'corrected', type = 'umap', plot = 'som',
                     markers = lineage_markers_2)
(rasterise(plot3,dpi = 300) +facet_wrap(~Batch))

cowplot::plot_grid(  rasterise(plot1,dpi = 300) ,  rasterise(plot2,dpi = 300))
cowplot::plot_grid(  rasterise(plot1,dpi = 300) ,  rasterise(plot3,dpi = 300)) %T>% 
  figsave('batch_gd/Uncorrected_vs_batchcorrected.pdf', 400, 200)
plot1$data

(rasterise(plot1,dpi = 300) +facet_wrap(~Batch)) %T>% 
  figsave('batch_gd/Uncorrected_facet4.pdf', 400, 200)
(rasterise(plot2,dpi = 300) +facet_wrap(~Batch))%T>% 
  figsave('batch_gd/corrected_facet4.pdf', 400, 200)
(rasterise(plot3,dpi = 300) +facet_wrap(~Batch))%T>% 
  figsave('batch_gd/corrected_facet4_fewermarker.pdf', 400, 200)

colnames(uncorrected)
colnames(corrected)

plot3$data

head(corrected)

ncol(uncorrected)
head(uncorrected)[,3:27]


rng <- colQuantiles(as.matrix(corrected[3:28]), probs = c(0.01, 0.99))

##normalize the expression matrix by :
#expression value - minimal expression value across all cells / data range
expr01 <- t((t(as.matrix(corrected[3:28])) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1




gdT_corrected <- cbind(plot3$data[,1:2], corrected, 
                          rename_with(uncorrected[,3:28] , ~paste0(., '_uncorrected')),
                          rename_with(as.data.frame(expr01) , ~paste0(., '_scaled'))
                          ) %>%
  mutate(sample = paste0(batch, '_',condition)) %>% 
  as.data.frame()

head(gdT_corrected)

# gdT_corrected$sample <- paste0(gdT_corrected$batch, '_', gdT_corrected$condition)

Feature_rast(gdT_corrected,'condition' ,d1 = 'UMAP1', d2 = 'UMAP2')


Feature_rast(gdT_corrected, 'condition', d1 = 'UMAP1', d2 = 'UMAP2', facets  = 'batch', do.label = F) %T>% 
  figsave('gdT_four_donor_tissue_UMAP.pdf', 200, 200)

Feature_rast(gdT_corrected, 'sample', d1 = 'UMAP1', d2 = 'UMAP2',  do.label = T) %T>% 
  figsave('gdT_allsamples.pdf', 200, 150)

Feature_rast(gdT_corrected, all_markers, d1 = 'UMAP1', d2 = 'UMAP2', ncol = 4, sz = 0.3) %T>% 
  figsave('gdT_allmarkers.pdf', 300, 400)


Feature_rast(gdT_corrected, paste0(all_markers, '_scaled'), d1 = 'UMAP1', d2 = 'UMAP2', ncol = 4, sz = 0.3) %T>% 
  figsave('gdT_allmarkers_scaled.pdf', 300, 400)



# SOM clustering ----------------------------------------------------------

fsom <- ReadInput(fcs_gdT_transform , transform = FALSE, scale = F)
fsom
# fsom_ns <- ReadInput(FCS_humangdT_biexp)

# replace the data in fsom with batch corrected data. 
fsom$data %>% colnames() 
str(fsom)

typeof(fsom$data)

fsom$scale

dim(corrected)
colnames(corrected)
panel
abcolors <- fsom$data %>% colnames() %>%  as.vector()
abcolors

cdata <- corrected[, 3:28] %>%  `colnames<-`(abcolors) %>% as.matrix()

fsom$data <- cdata 

# build up som 
# color to use 

lineage_colors <- panel %>%  dplyr::filter(desc %in% lineage_markers_2) %>%  pull(name) %>% as.vector()
as.data.frame(cdata) %>%  colnames()


set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_colors)
head(som)
# Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)
nmc <- 15
codes <- som$map$codes
t(codes) %>% head
dir.create('som_corrected_gd')

dir.create('som_corrected_gd/consensus_plots')

plot_outdir <- "som_corrected_gd/consensus_plots"


mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 1000,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)
head(mc)

code_clustering1 <- mc[[15]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]




# plot_outdir <- "FCS/exp1125/consensus_plots"
# codes_ns <- som_ns$map$codes
# 
# mc <- ConsensusClusterPlus(t(codes_ns), maxK = nmc, reps = 1000,
#                            pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
#                            clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
#                            distance = "euclidean", seed = 1234)
# 
# 
# ## Get cluster ids for each cell
# # choose 10 clusters according to the plot above
# 
# 
# code_clustering2 <- mc[[15]]$consensusClass
# cell_clustering2 <- code_clustering2[som_ns$map$mapping[,1]]


length(cell_clustering1)


# marker heatmap ----------------------------------------------------------


###build the function for heatmap
plot_clustering_heatmap_wrapper <- function(expr, expr01,
#alternatively
                                            cell_clustering, color_clusters, cluster_merging = NULL){

  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(funs(median))

  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))

  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  #the heating value of heat map
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  #mark the percentage of each cluster
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)

  # Row annotation for the heatmap
  # I think cluster names can be changed here
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)

  ###if want to change color look here?
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  ####here to give new cluster names & merge clusters
  #here in comment is how-to
  # test <- data.frame(old_cluster =  c(1:10), new_cluster = c("A","A","A","B", "A","C","C","C","D", "E"))
  # test$new_cluster <- factor(test$new_cluster)
  # annotation_row$cluster_mergering <- test$new_cluster
  # annotation_row


  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  # Colors for the heatmap
  color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)

  pheatmap::pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 8,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)

}
heatmap_orig_hcluster <- plot_clustering_heatmap_wrapper(expr = corrected[,lineage_markers_2],
                                           expr01 = expr01[,lineage_markers_2],
                                           cell_clustering = cell_clustering1,
                                           color_clusters = ggplotColours(20))%>% ggplotify::as.ggplot()

heatmap_orig_hcluster
figsave(heatmap_orig_hcluster, 'gdT_heatmap.pdf', 200, 150)

gdT_corrected$som_cluster <- as.factor(cell_clustering1)
gdT_corrected$tissue <-  gdT_corrected$condition
table(gdT_corrected$som_cluster, gdT_corrected$tissue)

table(tissue)
Feature_rast(gdT_corrected, 'som_cluster',d1 = 'UMAP1', d2 = 'UMAP2', facets = 'batch')
Feature_rast(gdT_corrected, c('som_cluster', 'batch', 'tissue', 'sample'),d1 = 'UMAP1', d2 = 'UMAP2', ncol  =2, sz = 0.3, do.label = F) %T>% figsave('gdT_UMAP_cl_batch_tissue_sample.pdf', 300, 200) 
Feature_rast(gdT_corrected, c('som_cluster'),d1 = 'UMAP1', d2 = 'UMAP2', ncol  =2, sz = 0.8, do.label = T) %T>% figsave('gdT_UMAP_som_cluster.pdf', 200, 200) 

saveRDS(gdT_corrected, 'gdTcell_processed_dataframe.rds')


# 
# heatmap_orig_merge <- plot_clustering_heatmap_wrapper(expr = expr01,
#                                                          expr01 = expr01,
#                                                          cell_clustering = cell_clustering_merging,
#                                                          color_clusters =clm ) %>% ggplotify::as.ggplot()
# 
# 



#alternatively, with merging clusters
#set the dataframe for cluser merging



## Find and skip duplicates
dups <- which(!duplicated(corrected[, lineage_markers]))
dups
## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)
inds

## How many cells to downsample per-sample
set.seed(1234)
tsne_ncells <- pmin(table(sample_ids), 10000)
## Get subsampled indices
set.seed(8964)

tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, all_markers]
tsne_expr %>% head

expr

library(umap)

## umap embedding of 150 items in 2 dimensions
## object components: layout, data, knn, config


# Umap --------------------------------------------------------------------
library(reticulate)

#Run tSNE
#position info stored in the tsne_out$Y
set.seed(1234)
# tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE )
tsne_expr

umap_out <-umap(tsne_expr, method = 'umap-learn')
umap_out2 <-umap(tsne_expr, method = 'naive')
umap_out3 <- uwot::umap(tsne_expr, scale = "Z", n_neighbors = 10, ret_model = T, min_dist = 0.1)
uwot::save_uwot(umap_out3, paste0(dir_FACS,'umapmodel'))

plot(umap_out3$embedding)


umap_out$layout %>% head
umap_out2$layout %>% head

umap_out3[,1]
saveRDS(umap_out, 'umapraw.rds')

dim(umap_out3[,1])
umap_out3$embedding

# Umap data.frame for ggplot ----------------------------------------------


umap_out$layout %>% head
dr_umap <- data.frame(UMAP_1 = umap_out3$embedding[,1] , UMAP_2 = umap_out3$embedding[,2] ,
                      
                 expr[tsne_inds, all_markers])

#put scaled data in 
dr_umap <- dr_umap %>% dplyr::select_at(vars(-ends_with('_scaled')))


dr_umap <- data.frame(dr_umap,
                      expr01[tsne_inds, all_markers] %>% 
                        `colnames<-`(paste0(all_markers, '_scaled')))
head(dr_umap)

dr_umap <- data.frame(dr_umap, UMAP_1 = umap_out3$embedding[,1], UMAP_2 = umap_out3$embedding[,2])


##add TCR
hist(dr_umap$Vd1_scaled, breaks = 200)
hist(dr_umap$Vd2_scaled, breaks = 200, xlim = c(0,0.5))


ggplot(dr_umap, aes(x = Vd1, y = Vd2, color =  TCRgd)) + geom_point()

head(dr_umap)
###gate on TCR 
dr_umap <- dr_umap %>% mutate(TCRD = case_when(
    Vd1_scaled >= 0.2 & Vd2_scaled < 0.2  ~ "TRDV1",
    Vd2_scaled >= 0.2& Vd1_scaled < 0.35 ~ "TRDV2"
)) %>% mutate(TCRD = replace(TCRD, is.na(TCRD), 'Other TRDV'))
dr_umap$TCRD <- factor( dr_umap$TCRD ,c('TRDV1', 'TRDV2', 'Other TRDV'))
ggplot(dr_umap, aes(x = Vd1, y = Vd2, color =  TCRD)) + geom_point()

##gate on CD26 and CD161 
ggplot(dr_umap, aes(x = CD26, y = CD161)) + geom_point(color = alpha('grey', 0.5)) +
  theme_minimal()+geom_hline(yintercept = 7)+geom_vline(xintercept = 6.5)

ggplot(dr_umap, aes(x = CD26, y = CCR6)) + geom_point(color = alpha('grey', 0.5))+
  theme_minimal()+geom_hline(yintercept = 6.5)


dr_umap <- dr_umap %>% mutate(gateCD26 = case_when(
  CD26 >= 6.5 ~ 'CD26pos',
  CD26 < 6.5 ~ 'CD26neg'),
gateCD161 = case_when(CD161 >= 7.1 ~'CD161pos',
                      CD161 < 7.1 ~'CD161neg'),

)

dr_umap <- dr_umap %>% mutate(gateCCR6 = case_when(
  CCR6 >= 6.5 ~ 'CCR6pos',
  CCR6 < 6.5 ~ 'CCR6neg')
  
)


saveRDS(dr_umap, dr_umapRDS)

#We can color the cells by cluster. Ideally, cells of the same color should be close to each other (see Figure 5).


dr_umap$sample_id <- sample_ids[tsne_inds]
mm <- match(dr_umap$sample_id, md$file_name)
dr_umap$donor <- md$donor[mm]
dr_umap$protocol <- md$process[mm]
dr_umap$sti <- md$stimulation[mm]
dr_umap$cell_clustering <- factor(cell_clustering1[tsne_inds], levels = 1:23)


Feature_rast(dr_umap, 'cell_clustering')

head(dr_umap)

Feature_rast(dr_umap, all_markers, ncol = 6, sz = 0.1) %T>%  figsave('allmarkers.pdf', 300, 200)
Feature_rast(dr_umap, paste0(all_markers, '_scaled'), ncol = 6,sz = 0.1) %T>%  figsave('allmarkers_scaled.pdf', 300, 200)

paste0(all_markers, '_scaled')


Feature_rast(dr_umap,c('cell_clustering', 'donor', 'sti', 'protocol'), sz = 0.1, ncol = 2, do.label = F, noaxis = F)%T>%  
  figsave('UMAP_clustering.pdf', 200, 200)


plot2 <- plot_dimred(corrected, name = 'Corrected 8x8', type = 'umap')


dr_umap$UMAP_1 <- plot2$data$UMAP1
dr_umap$UMAP_2 <- plot2$data$UMAP2
colnames(corrected)
dr_umap_cr <- corrected[,3:33]
dr_umap_cr$UMAP_1 <- plot2$data$UMAP1
dr_umap_cr$UMAP_2 <- plot2$data$UMAP2
dim(dr_umap)

Feature_rast(as.data.frame(dr_umap_cr), all_markers, ncol = 4, sz = 0.1) %T>%  figsave('allmarkers_adjusted.pdf', 300, 200)

ggplot(dr_umap_cr, aes(x = sample_id, y = Vd1, fill = sample_id ))+ geom_jitter_rast(size = 0.3)


Feature_rast(as.data.frame(dr_umap_cr), 'CD4', facets = 'sti')

Feature_rast(as.data.frame(dr_umap_cr), c('sti', 'sample_id' , 'donor'), sz = 0.2)

# Feature_rast(dr_umap, all_markers, ncol = 6, sz = 0.1) %T>%  figsave('allmarkers_adjusted.pdf', 300, 200)


# feature map -------------------------------------------------------------
#proten to gene 
ptg <- c('CCR6', "CD26 (DPP4)", "CD161 (KLRB1)", "CD127 (IL7R)", 'CCR4', "CCR7", "CD94 (KLRD1)", "CD16 (FCGR3A)", 
         "PD1 (PDCD1)", "TCR-VD1", "TCR-VD2")


# umap_allmarkers <- map(c(lineage_markers, functional_markers) %>% as.list(), function(x) {
#   ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
#                               color = x)) +
#     geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_classic() +
# theme(legend.key.width = unit(1, 'mm'))+mytheme+
#    grd
# }
# ) %>% set_names(c(lineage_markers, functional_markers) %>% as.list())

umap_allmarkers$PD1+ facet_wrap(~sample_id)
umap_allmarkers$CCR4 + facet_wrap(~sample_id)

umap_allmarkers$CCR4

# umap_allmarkers_2 <- map(c(lineage_markers, functional_markers) %>% as.list(), function(x) {
#   ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
#                               color = x)) +
#     geom_point_rast(size = 0.5, raster.dpi = 300) +
#     theme_classic() +theme(legend.key.width = unit(1, 'mm'))+mytheme+
#     scale_color_gradient2(low = alpha('lightgrey', 0.3), high = 'red', mid = 'purple', 
#                           midpoint =quantile(dr_umap[[x]],0.8 ) )
# }
# ) %>% set_names(c(lineage_markers, functional_markers) %>% as.list())
# umap_allmarkers_2$CD16

umap_allmarkers_tc <- map2(c(lineage_markers, functional_markers) %>% as.list(), ptg, function(x,y) {
  ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                              color =paste0(x, '_scaled') )) +
    theme_classic() +
   ggtitle(y)+
    geom_point_rast(size = 0.2, raster.dpi = 300) +
    theme(legend.key.width = unit(1, 'mm'))+mytheme+
     NoAxes()+ylim(-8, 6.1)+
    guides(color = guide_colorbar(title.position = 'left', title ="scaled expression" ,
                                  title.theme = element_text(angle = 90, size = 8)))+
    scale_color_gradient2(low = alpha('lightgrey', 0.3), high = 'red', mid = 'purple',
                          midpoint = 0.65    )
}
) %>% set_names(c(lineage_markers, functional_markers) %>% as.list())
umap_allmarkers_tc$CCR4

PG(umap_allmarkers_tc)

quantile(dr_umap$CCR6_scaled, c(0.1,0.9))

hist(dr_umap$CCR6_scaled)

hist(dr_umap$CCR6_scaled, breaks = 100)

mean(dr_umap$PD1)
max(dr_umap$PD1)

quantile(dr_umap$PD1)


umap_allmarkers_tc

umap_allmarkers$CCR4

ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                            color = "CCR4")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +theme(legend.key.width = unit(1, 'mm'))+mytheme+
  scale_color_gradient2(low = alpha('lightgrey', 0.3), high = 'red', mid = 'purple', midpoint = 6.5)


umap_allmarkers$Vd1


umap_markers

umap_markers <- PG(umap_allmarkers, ncol = 3 )
umap_markers
umap_markers_2 <- PG(umap_allmarkers_2, ncol = 3 )
umap_markers_scaled <- PG(umap_allmarkers_tc, ncol = 3 )


ggsave2(filename = "FACS_umap_markers.pdf", umap_markers, width = 270, height = 200, units = 'mm')
ggsave2(filename = "FACS_umap_markers_2.pdf", umap_markers_2, width = 270, height = 200, units = 'mm')
ggsave2(filename = "FACS_umap_markers_scaled.pdf", umap_markers_scaled, width = 240, height = 200, units = 'mm')
# TCRD and group --------------------------------------------------------------------

TRDFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                            color = "TCRD")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +mytheme+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  guides(color = guide_legend( override.aes = list(size = 1.5))) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))

groupFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                            color = "condition")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +theme(legend.key.width = unit(1, 'mm'))+mytheme+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  guides(color = guide_legend( title = 'group' ,override.aes = list(size = 1.5))) +
  scale_color_manual(values = c("#ff6633", "#00cc66"))

donorFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                                         color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +theme(legend.key.width = unit(1, 'mm'))+mytheme+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  guides(color = guide_legend( override.aes = list(size = 1.5))) 

PG(list(TRDFACS, groupFACS))






# #We can color the cells bycluster. --------------------------------------


#We can color the cells by cluster. Ideally, cells of the same color should be close to each other (see Figure 5).

dr_umap$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)
dr_umap$cell_clustering2 <- factor(cell_clustering2[tsne_inds], levels = 1:nmc)
dr_umap$cell_clustering3 <- factor(cell_clustering3[tsne_inds], levels = 1:nmc)
dr_umap$SOM <- factor(SOM_clustering[tsne_inds], levels = 1:100)

dr_umap$cell_hcluster <- factor(cell_hcluster[tsne_inds], levels = 1:20)
dr_umap$cell_hcluster_2 <- factor(cell_hcluster_2[tsne_inds], levels = 1:20)
dr_umap$cell_hcluster_3 <- factor(cell_hcluster_3[tsne_inds], levels = 1:20)
dr_umap$cell_hcluster_4 <- factor(cell_hcluster_4[tsne_inds], levels = 1:20)


saveRDS(dr_umap, file = 'FCS/exp1125/dr_umap.rds')
saveRDS(ggdf, file = 'ggdf.rds')



ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = SOM)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  scale_color_manual(values = sample(ggplotColours(100)))+
  gglp('b')+
  theme_classic()+
  # facet_wrap(~cell_clustering3)+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

set_sample()

## Plot t-SNE colored by clusters
ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cell_clustering3)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) + ylim(-8,6.1)+
  scale_color_manual(values = set_sample(ggplotColours(19)))+
  
  theme_classic()+
  # facet_wrap(~cell_clustering3)+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))






ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = TCRD)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +

  theme_classic()+
  # facet_wrap(~cell_clustering3)+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))



ggplot(dr_umap,  aes(x = Vd1, y = Vd2, color = cell_clustering2)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic()+facet_wrap(~cell_clustering2)+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))


Umap_orig_hc <- ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cell_hcluster)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +
  scale_color_manual(values = sample(ggplotColours(50)))+
  geom_text(data = (dr_umap %>% group_by(cell_hcluster) %>% 
                      dplyr::select(UMAP1,  UMAP2) %>% 
                      summarize_all(mean) %>% dplyr::rename(center = "cell_hcluster")), 
            aes(x = UMAP1, y = UMAP2, label = center),color = 'black')+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
Umap_orig_hc

Umap_orig_som <- ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cell_clustering1)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +
  scale_color_manual(values = sample(ggplotColours(18)))+
  geom_text(data = (dr_umap %>% group_by(cell_clustering1) %>% 
                      dplyr::select(UMAP1,  UMAP2) %>% 
                      summarize_all(mean) %>% dplyr::rename(center = "cell_clustering1")), 
            aes(x = UMAP1, y = UMAP2, label = center),color = 'black')+

  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
Umap_orig_som



som2.cent = dr_umap %>% group_by(cell_clustering2) %>% dplyr::select(UMAP1, 
                                                                       UMAP2) %>% 
  summarize_all(mean) %>% dplyr::rename(center = "cell_clustering2")


Umap_orig_som2 <- ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cell_clustering2)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +
  scale_color_manual(values = sample(ggplotColours(15)))+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))+
  geom_text(data = som2.cent, aes(x = UMAP1, y = UMAP2, label = center),color = 'black')
clustersom2_donor <-    Umap_orig_som2+facet_wrap(~sample_id)
Umap_orig_som2
ggsave2(clustersom2_donor, filename = 'clustersom2_donor.pdf', width = 200, height = 200, units = 'mm')
ggsave2(Umap_orig_som2, filename = 'Umap_orig_som2.pdf', width = 200, height = 200, units = 'mm')

Umap_orig_som2 + facet_wrap(~cell_clustering2)
##cluster adjustment based on som2

# ##cluster adjustment based on som2 --------------------------------------


# dr_umap$cluster_adj <- NULL
# dr_umap <- dr_umap %>% mutate(cluster_merge = 
#                                 case_when(cell_clustering2 == c("7", "9") ~ "Vd2 CD94+ IL7R+ developing CTL",
#                                           cell_clustering2 == c('11', '14') ~ "Vd2 CD94+ CD16+  CTL",
#                                           cell_clustering2 == c("1", "2")~ "Vd2 CCR6+ CD26+ gd17",
#                                           cell_clustering2 == "3" ~ "Vd2 CCR4+IL7R+  developing",
#                                           cell_clustering2 == c("12", "13") ~ "Vd1 CD16+ CTL", 
#                                           cell_clustering2 == c("8", "15") ~"unidentified",
#                                           cell_clustering2 == c("4", "5")  ~ "mixed CCR7+ PD1+/- naive",
#                                           cell_clustering2 == c( "10") ~ "Vd1 CD16- PD1+  CTL",
#                                           cell_clustering2 == "6" ~ "Vd2 CCR4- developing"))

cluster_merging <- data.frame(old_cluster =  c(1:15), new_cluster = c("Vd2 CCR6+ CD26+ gd17", #1
                                                                      "Vd2 CCR6+ CD26+ gd17", #2
                                                                      "Vd2 CCR4+IL7R+  developing", #3
                                                                      "mixed CCR7+ PD1+/- naive", #4
                                                                      "mixed CCR7+ PD1+/- naive",#5
                                                                      "Vd2 CCR4- developing", #6
                                                                      "Vd2 CD94+ CD16- developing CTL", #7
                                                                      "unidentified",#8
                                                                      "Vd2 CD94+ CD16- developing CTL", #9
                                                                      "Vd1 CD16- PD1+  CTL",#10
                                                                      "Vd2 CD94+ CD16- developing CTL",#11
                                                                      "Vd1 CD16+ CTL", #12
                                                                      "Vd1 CD16+ CTL", #13
                                                                      "Vd2 CD94+ CD16+  CTL",#14
                                                                      "unidentified"#15
                                                                      ))




cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
mm <- match(cell_clustering2, cluster_merging$old_cluster)

cell_clustering_merging <- factor(cluster_merging$new_cluster[cell_clustering2],
                                  levels = as.vector(unique(cluster_merging$new_cluster[cell_clustering2]))[c(1,5,3,4,6,7,9,2,8)])
dr_umap$cluster_merge <- factor(cell_clustering_merging[tsne_inds])




dr_umap <- dr_umap %>% mutate(cluster_merge = 
                     case_when(
                                             cluster_merge ==   "Vd2 CCR4+IL7R+  developing" ~ "Vd2 CCR4+ NKT",
                                             cluster_merge ==  'Vd2 CCR6+ CD26+ gd17' ~ "Vd2 NKT-17", 
                                             cluster_merge == as.vector(unique(dr_umap$cluster_merge)[4]) ~ "Vd1 Exhauster CTL",
                               cluster_merge == as.vector(unique(dr_umap$cluster_merge)[5]) ~ 'Vd2 CCR7+ NKT-1naive',
                               cluster_merge == as.vector(unique(dr_umap$cluster_merge)[6])~ "Vd1 Active CTL",
                               cluster_merge == as.vector(unique(dr_umap$cluster_merge)[7]) ~ "Vd2 IL7R+ NKT-1EF",
                               cluster_merge == as.vector(unique(dr_umap$cluster_merge)[8]) ~ "Vd2 CD16+ NKT-1EF",
                               !is.na(cluster_merge) ~ as.character(cluster_merge)))

dr_umap$cluster_merge <- factor(dr_umap$cluster_merge, 
                                levels = sort(unique(dr_umap$cluster_merge))[c(3:9,1,2)])


as.vector(unique(pull(test, cluster_merge)))
test$cluster_merge %>% unique() %>% as.vector() %>% sort

ggplot(dr_umap, aes(x = UMAP1, y = UMAP2)) + geom_point(aes(color = cluster_merge))
rm(test)

saveRDS(dr_umap, 'dr_umap.rds')


ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cluster_merge)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) 

##cluster merging based on cell_clustering 3
# ##cluster merging based on cell_clustering 3 ----------------------------
plotly::ggplotly((ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cell_clustering3))+
                    geom_point()+   scale_color_manual(values = set_sample(ggplotColours(19))) 
))


cluster_merging <- data.frame(old_cluster =  c( 6,13,
                                                15,
                                                14,
                                                12,
                                                7,9,10,
                                                1,2,
                                               3,4,
                                               5,
                                               8, 11,
                                               16
                                               ), 
                              new_clusterP = c(  rep( "CCR7hi Naive",2),
                                                 "CCR7int Naive",
                                                 'Vd1 CD16+ CTLs',
                                                 'Vd1 CD16- CTLs',
                                                 rep("Vd2 naive",3),
                                          rep('Vd2 Th-17 like',2),
                                           rep("Vd2 IL7R+ Th-1 like", 2),
                                            'Vd2 CD16+ Th-1 like',


                                        
                                           rep('Vd2 CCR4+',2),
                                          'unidentified'), 
                              new_clusterN = c('FCS1', "FCS1",
                                               "FCS2",  "FCS3",  "FCS4" ,
                                               "FCS5","FCS5","FCS5",
                                               "FCS6", "FCS6",
                                               "FCS7", "FCS7",
                                               "FCS8" ,
                                               "FCS9" ,"FCS9" ,
                                               "FCS10"
                                               )) %>% arrange(old_cluster)

mm <- match(cell_clustering3, cluster_merging$old_cluster)
lv <- sort(unique(cluster_merging$new_clusterP[cell_clustering3]))[c(1,2,4:10,3)]
cell_clustering_merging2 <- factor(cluster_merging$new_clusterP[cell_clustering3],
                                  levels = lv)
dr_umap$cluster_mergeP <- factor(cell_clustering_merging2[tsne_inds])

cell_clustering_mergingN <- factor(cluster_merging$new_clusterN[cell_clustering3],
                                   levels = paste0("FACS",1:10))

dr_umap$cluster_mergeN <- factor(cell_clustering_mergingN[tsne_inds])

dr_umap$cluster_mergeN <- gsub('FCS', 'FACS', dr_umap$cluster_mergeN) %>% factor(levels = paste0("FACS",1:10))


Umap_cl_adj <- ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cluster_mergeP)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +
scale_color_brewer(palette = 'Paired')+
    guides(color = guide_legend(override.aes = list(size = 1.5), ncol = 2))
Umap_cl_adj

Umap_cl_adj+facet_wrap(~sample_id)
Umap_cl_adj_group <- Umap_cl_adj + facet_wrap(~condition)
Umap_cl_adj_group
Umap_cl_adj + umap_markers_2
ggsave2(Umap_cl_adj, filename = 'Umap_cl_adj.pdf', width = 200, height = 200, units = 'mm')

heatmap_merged <- plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers],
                                                  expr01 = expr01[, lineage_markers],
                                                  cell_clustering = cell_clustering_merging2,
                                                  color_clusters = ggplotColours(10))

save_plot(paste0(dir_FACS,'heatmap_merged.pdf')   , heatmap_merged, base_height = 5, base_width = 8)






plotlyUmap<-  ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cluster_mergeP)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_manual(values = sample(ggplotColours(50)))+
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))  




plotly::ggplotly(plotlyUmap)
# save_plot('FACS_UMAP_CLUSTER15.pdf', Umap_orig, ncol = 1, base_height = 5, base_width = 5)


ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cell_clustering1)) +
  geom_point(size = 1.2) +scale_color_manual(values = set_sample(rainbow(13)))+facet_wrap(~ sample_id)
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))


#seprate by condition
tsne_condition <- tsne_merged + facet_wrap(~ condition)
save_plot('tsne_condition.png', tsne_condition, ncol = 2, base_height = 5, base_width = 5)

#seprate by samples
tsne_sample <- tsne_merged + facet_wrap(~ sample_id)
save_plot('tsne_sample.png', tsne_sample, ncol = 2, base_height = 5, base_width = 5)



#8 Differential analysis
library(lme4)
library(multcomp)
## Model formula without random effects
model.matrix( ~ condition, data = md)
## Create contrasts bewteen 2 conditions
contrast_names <- c("KO vs WT")
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
K
FDR_cutoff <- 0.05


###creat counts table & frequencies table
count_table <- table(cell_clustering1, condition)
counts <- as.data.frame.matrix(count_table)
Freq_table <- t(t(count_table)/colSums(count_table))*100
Freq <- as.data.frame.matrix(Freq_table)

##Compare cell composition of the 2 conditions
ggdf_2 <- melt(data.frame(cluster = rownames(Freq), Freq), id.vars = 'cluster',
             value.name = 'proportion', variable.name = 'sample_id')
ggdf_2
mm <- match(ggdf_2$sample_id, md$sample_id)
mm
ggdf_2$condition  <- md$condition[mm]


###integerating clustering into ggdf
library(reshape2)
dr_long <- dr
dr_long[,1:2 ] <- NULL
dr_long <- gather(dr_long, Marker, Expression, 1:17)
dr_long %>% head
nrow(dr_long)
nrow(ggdf)
nrow(dr)

head(ggdf)
ggdf_w <- dcast(ggdf,sample_id + condition ~ antigen, value.var = 'expression')
ggdf_w %>% head
conditons <- condition
condition <- NULL

# saveRDS(dr, file = 'dr.rds')
saveRDS(dr_umap, file = 'dr_umap.rds')

saveRDS(ggdf, file = 'ggdf.rds')
dr %>% group_by(cell_clustering1, condition) %>% summarise(median(GFP), mean(GFP), sd(GFP))



# ####fig7 single-cell data is valiadated by facs -------------------------

setwd("/home/large/likAIiiiiiiiiiii/Human_GDT_2019/")

dr_umap <- readRDS('FCS/exp1125/dr_umap.rds')
Umap_HumanGDT <- readRDS('Integrated/Umap_metadata.rds')

HumanGDT <- readRDS('Integrated/Humangdt_20190712.RDS')

# #Fig7A  -----------------------------------------------------------------


#Fig7A   summarize of gdT phenotypes
HumanGDT <- RenameIdents(object = HumanGDT, 
                        'c1' = 'Mixed naive',
                        'c3' = 'Mixed naive',
                        "c2"= 'Replicating naive',
                        'c5'= "CCR4+ NKT",
                        'c4'= "CCR4+ NKT",
                        'c6' = 'naive NKT',
                        'c7' = 'Acute activated',
                        'c8' = 'NKT-17',
                        'c9' = "INF induced",
                        'c10' = 'NKT-1',
                        'c11' = "Exhausted CTL",
                        "c12" = "Active CTL")



HumanGDT$pheno <- Idents(HumanGDT)
Idents(HumanGDT) <- HumanGDT$number_cl

Umap_HumanGDT$pheno <- NULL

Umap_HumanGDT$barcode
Umap_HumanGDT<- HumanGDT@meta.data %>% dplyr::select(bc_backup, pheno) %>% dplyr::rename( 'barcode'= 'bc_backup') %>% 
   right_join(Umap_HumanGDT, by = 'barcode')

DimPlot(HumanGDT, group.by = 'pheno')
saveRDS(HumanGDT, file = 'Integrated/Humangdt_20190712.RDS')
saveRDS(Umap_HumanGDT, file = 'Integrated/Umap_metadata.rds')

scpheno_c <- c("#F781BF", "#ff0066", "#4DAF4A", "#FFFF33", "#00ffff", "#984EA3" , "#FF7F00" ,"#006600",
                "#E41A1C","#377EB8" )

Fig7A_scpheno <- ggplot(Umap_HumanGDT, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(size = 0.2, raster.dpi = 300, aes(color = pheno))   + theme_classic() +
  theme(legend.position = "bottom", legend.key.height = unit(0.3, "cm"))+
  scale_color_manual(values = scpheno_c)+
  guides(color = guide_legend(nrow = 4, title = NULL ,override.aes = list(size = 1.5, alpha = 1))) +
  mytheme+notick
Fig7A_scpheno


Fig7A_scpheno_nolegend <- ggplot(Umap_HumanGDT, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(size = 0.2, raster.dpi = 300, aes(color = pheno))   + theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual(values = scpheno_c)+
  mytheme+notick


Fig7A_scpheno_nolegend





# Fig7B table -------------------------------------------------------------
library(gridExtra)
library(grid)
library(gtable)
phenoT <- read.csv('pheno.csv')

tt <- theme_minimal(core = list(bg_params = list(fill = alpha(scpheno_c,0.8)),
                                 fg_params=list(cex = 0.66)),
                     colhead = list(fg_params=list( cex=0.66) ))


pt <- grid.table(phenoT, theme = tt, rows = NULL) 

Fig7A_scpheno+
  gridExtra::tableGrob(phenoT, theme = tt, rows = NULL) 



PG(list(Fig7A_scpheno,   gridExtra::tableGrob(phenoT, theme = tt, rows = NULL) ))

# fig7C feature map -------------------------------------------------------

Fig7C_facsfeature <-
Reduce("+", umap_allmarkers_tc[1:9]) + plot_layout(guides = 'collect') 

fig7clegend <- get_legend(umap_allmarkers_tc$CCR6)

umap_allmarkers_tc_nolegend <- map(umap_allmarkers_tc[1:9], function(x) x + theme(legend.position = 'none'))

umap_allmarkers_tc_nolegend$legend <- fig7clegend

Fig7C_facsfeature <- PG(umap_allmarkers_tc_nolegend, nrow = 3, ncol = 4)

PG(umap_allmarkers_tc_nolegend)

# Fig7D&E  & TRD? ------------------------------------------------------


groupFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                         color = "condition")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +theme(legend.key.width = unit(1, 'mm'))+mytheme+
  theme(axis.ticks = element_blank(), axis.text = element_blank())+
  guides(color = guide_legend( title = 'group' ,override.aes = list(size = 1.5))) +
  scale_color_manual(values = c("#ff6633", "#00cc66"))


donorcolor <- c("#003300", "#33cc33","#99ff99","#00cc99", "#669900", "#339933", #CB, green
                "#993333", "#ff0000", "#993300", "#ff6699", "#ff0066","#cc0099" #PB, red
)

Fig7D_donorFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                         color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +
  theme(legend.key.size = unit(3, 'mm'))+mytheme+notick+
  guides(color = guide_legend(title = 'donor', override.aes = list(size = 1.5))) +
  scale_color_manual(values = donorcolor)

Fig7D_donorFACS_nolegend <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                               color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +
  theme(legend.position = 'none')+mytheme+notick+
  scale_color_manual(values = donorcolor)


Fig7E_TRDFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                       color = "TCRD")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +mytheme+notick+
  guides(color = guide_legend( title = 'TCRD', override.aes = list(size = 1.5))) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))
Fig7E_TRDFACS_legendlow <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                             color = "TCRD")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +mytheme+notick+
  theme(legend.position = 'bottom', legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))

Fig7E_TRDFACS_legendinner <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                                       color = "TCRD")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +mytheme+notick+
  theme(legend.position = c(0.2,0.15), legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))

Fig7E_TRDFACS_legendinner
# fig7F clustering  -------------------------------------------------------

c("#F781BF", "#ff0066", "#4DAF4A", "#FFFF33", "#00ffff", "#984EA3" , "#FF7F00" ,"#006600",
  "#E41A1C","#377EB8" )

library(RColorBrewer)
display.brewer.pal(n = 9, name = 'Set1')
levels(dr_umap$cluster_merge)
brewer.pal(n = 9, name = 'Set1')
cl_color <- brewer.pal(n = 9, name = 'Set1')[c(2,1,3,6,7,5,4,8,9)]


fig7F_Umap_cl_adj <- ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cluster_merge)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +mytheme+ notick+
  theme(legend.position = 'bottom', legend.key.size = unit(2, 'mm'))+
  scale_color_manual(values = cl_color )+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1.5), ncol= 2))
fig7F_Umap_cl_adj


ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cluster_merge)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +mytheme+ notick+
  theme(legend.position = 'bottom', legend.key.size = unit(2, 'mm'))+
  scale_color_manual(values = cl_color )+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1.5), ncol= 2))



ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cluster_merge)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +mytheme+ notick+
  theme(legend.position = 'right')+
  scale_color_brewer(palette = 'Set1')+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1.5), ncol = 1))




fig7F_Umap_cl_adj_R <- ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = cluster_merge)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +mytheme+ notick+
  theme(legend.position = 'right')+
  scale_color_brewer(palette = 'Set1')+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1.5), ncol = 1))


ggplot(dr_umap,  aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +
scale_color_manual(labels = c('Neonate', 'Adult'), values = c('blue', 'red'))+
  theme_classic() 

# fig7G donor composition of cluster and TCR ------------------------------

dr_cl_comp <- dr_umap %>% dplyr::group_by(condition, sample_id, cluster_merge) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n)*100) %>% 
  dplyr::rename('cluster' = 'cluster_merge','donor'= "sample_id" , 'group'='condition' )
  
levels(dr_cl_comp$cluster)

dr_cl_comp$cluster <- factor(dr_cl_comp$cluster, 
                             levels = unique(dr_cl_comp$cluster)[c(8, 3,4,7, 6,5, 1, 2,9)] )



FigS8_comp<-  ggplot(dr_cl_comp, aes(x = cluster, y = percent, fill = factor(group))) + 
  geom_boxplot(  alpha = 0.4 , color = 'black', size = 0.3,
                width = 0.5,
                outlier.alpha = 0)+ylab('% of donor')+xlab(NULL)+
  scale_fill_manual(values = c("#00ff99","#ff66ff" ))+
  geom_jitter(aes(color = factor(donor)),position = position_dodge(width = 0.6), size = 1) + 
  scale_color_manual(values = donorcolor) +  
  theme_minimal()+mytheme+
  theme(axis.text.x = element_text(angle = 320, hjust = 0), legend.key.size = unit(4, 'mm')
  )+  guides(color = guide_legend(title = "donor", override.aes = list(size = 1.5), ncol = 2),
             fill =  guide_legend(title = "group"))




FigS8_comp



# # assemble --------------------------------------------------------------
# Fig7AB <- PG(list(Fig7A_scpheno_nolegend, grid.table(phenoT, theme = tt, rows = NULL) ), 
#              labels = c("A", "B"), rw = c(1,2), nrow = 1)
# FIg7CDE <- PG(list(Fig7C_facsfeature, Fig7D_donorFACS, Fig7E_TRDFACS), nrow = 1, rw = c(1.5, 1, 1),
#               labels = c("C", "D", "E"))
# Fig7FG <- PG(list(fig7F_Umap_cl_adj, Fig7G_comp), nrow = 1, rw = c(1, 2), labels = c('F', 'G'))


# Fig7 <- PG(list(Fig7AB, FIg7CDE, Fig7FG), ncol = 1, rh = c(1,1.3,1.2))
# ggsave2(Fig7, filename = 'manu_fig/fig7.test.2.pdf', height = 290, width = 200, units = 'mm')



##or 
Fig7A  <- PG(list(Fig7A_scpheno_nolegend,  gridExtra::tableGrob(phenoT, theme = tt, rows = NULL)), 
             labels = c("A", NA), rw = c(1,2.2), nrow = 1)
# Fig7BC <- PG(list(Fig7C_facsfeature, PG(list(fig7F_Umap_cl_adj_R, NA), ncol = 1, rh = c(1, 0.25))), 
#              labels = c('B', "C"), rw = c(0.85,1))
# Fig7DE <- PG(list(PG(list(Fig7E_TRDFACS_legendlow + NoAxes(), Fig7D_donorFACS_nolegend+ NoAxes()), labels = c('D', "E"),
#                      axis = 'tblr', align = 'hv'),
#                   Fig7G_comp), nrow = 1,
#              rw = c(1,1,2.4))
# Fig7 <- PG(list(Fig7A, Fig7BC, Fig7DE), ncol = 1,
#            rh = c(1, 1.3, 1.3))

#or
Fig7BCD <- PG(list(Fig7C_facsfeature, 
                   PG(list(fig7F_Umap_cl_adj, Fig7E_TRDFACS_legendlow), labels = c("C", "E"), axis = 'tblr', align = 'hv')
                   ), 
             labels = c('B'), rw = c(0.85,1))
Fig7BCD

Fig7EF <- PG(list(Fig7D_donorFACS_nolegend, Fig7G_comp), labels = c('E', NA),
             axis = 'bl', align = 'hv', rw = c(0.6,1.3))
Fig7EF
Fig7D_donorFACS_nolegend

Fig7 <- PG(list(Fig7A, Fig7BCD, Fig7EF), ncol = 1,
           rh = c(1, 1.3, 1.3))

ggsave2(Fig7, filename = 'manu_fig/fig7.191215.pdf', height = 2650, width = 190, units = 'mm')



# Fig7 new  ---------------------------------------------------------------

##A clustering 

Fig7A_clustering <- ggplot(dr_umap,  aes(x = UMAP_1, y = UMAP_2, color = cluster_mergeN)) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +
  theme_classic() +mytheme+ notick+ylim(-8,6.1)+
  theme(legend.position = 'right', legend.key.size = unit(5, 'mm'))+
  scale_color_brewer(palette = 'Paired')+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 1.5), ncol= 1))

Fig7A_clustering

##B features 

facsfeature <-
  Reduce("+", umap_allmarkers_tc[1:9]) + plot_layout(guides = 'collect') 

featurelegend <- get_legend(umap_allmarkers_tc$CCR6)

umap_allmarkers_tc_nolegend <- lapply(umap_allmarkers_tc[1:9], 
                                   function(x) {x + theme(legend.position = 'none')})
umap_allmarkers_tc_nolegend$CCR6
umap_allmarkers_tc_nolegend$legend <- featurelegend

facsfeature <- PG(umap_allmarkers_tc_nolegend, nrow = 2, ncol = 5)
facsfeature
PG(umap_allmarkers_tc_nolegend)


##C TCR 
TRDFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                                             color = "TCRD")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +
  mytheme+notick+ylim(-8,6.1)+
  guides(color = guide_legend( title = 'TCRD', override.aes = list(size = 1.5))) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))
#D donor
donorcolor <- c("#003300", "#33cc33","#99ff99","#00cc99", "#669900", "#339933", #CB, green
                "#993333", "#ff0000", "#993300", "#ff6699", "#ff0066","#cc0099") #PB, red


donorcolor2 <- c("#003300", "#33cc33","#99ff99","#00cc99", "#669900", "#339933", #CB, green
                "#660066", "#cc00cc", "#ff66ff", "#6600cc", "#cc99ff","#cc0099") #PB, purple

donorFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP_1', y = "UMAP_2",
                                         color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +ylim(-8,6.1)+
  theme(legend.key.size = unit(5, 'mm'), legend.position = 'right')+mytheme+notick+
  guides(color = guide_legend(title = 'donor', override.aes = list(size = 1.5), ncol = 2)) +
  scale_color_manual(values = donorcolor2)
                
donorFACS
#summerise table
library(gridExtra)
library(grid)
library(gtable)
phenoT <- read.csv('pheno_new.csv') 


phenoT$Phenotype<- gsub('Innate-naïv', 'Innate-naïve', phenoT$Phenotype)

tt <- cowplot::theme_minimal(core = list(fg_params=list(cex = 0.66)),
                    colhead = list(fg_params=list( cex=0.66) ))



theme_minimal_grid()

pt <- grid.table(phenoT, theme = tt, rows = NULL) 



##assemble

#Fig7new 
Fig7NEWAB <- list(Fig7A_clustering, facsfeature) %>% PG(labels = "auto", ncol =2, rw = c(1,1.4))
Fig7NEWCD <- list(TRDFACS, donorFACS,NA) %>% PG(labels = c('c','d',NA), rw = c(1,1.3,0.05), 
                                             ncol = 3, align = 'hv')

Fig7NEW <- list(Fig7NEWAB, Fig7NEWCD, gridExtra::tableGrob(phenoT, theme = tt, rows = NULL)) %>%
  PG(labels = c(NA,NA, "e"), nrow = 3, rh = c(1,1,1.8))+
  draw_figure_label('Figure 8', position = 'top.right', size = 10, fontface = 'plain')

Fig7NEW

figsave(Fig7NEW, 'fig8_new_20200629.pdf', h = 260, w = 190)
ggsave2(Fig7NEW, filename = 'manu_fig/fig8_new_20200609.pdf', height = 260, width = 190, units = 'mm')

layout <- "
AABBB
AABBB
CCDD#
CCDD#

"




Fig7new <- (Fig7A_clustering + facsfeature +TRDFACS + donorFACS +
  plot_layout(design  = layout)) /
  gridExtra::tableGrob(phenoT, theme = tt, rows = NULL) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10, face = 'bold'))

ggsave2(Fig7new, filename = 'manu_fig/fig7_new_20200522.pdf', height = 260, width = 190, units = 'mm')

# FigS7 -------------------------------------------------------------------
####surface deg from part1 DEGs 

FigS7_heat_surface <- DoHeatmap(subset(HumanGDT, downsample = 800), size = gs(8), 
                          angle = 45, group.by = 'pheno', group.colors = scpheno_c,
                          features = top10surface_list, raster = T) +
  theme(text = element_text(size = 8), axis.text = element_text(size = 8, face = 'italic'),
        legend.key.width = unit(2,'mm'))+guides(color = FALSE)+
  hmp
FigS7_heat_surface

ggsave2(FigS7_heat_surface, filename = 'manu_fig/Supfig7_surfaceheat_20200502.pdf', width = 190, height = 210, units = 'mm')



# FigS8  ------------------------------------------------------------------

donorcolor <- c("#003300", "#33cc33","#99ff99","#00cc99", "#669900", "#339933", #CB, green
                "#993333", "#ff0000", "#993300", "#ff6699", "#ff0066","#cc0099" #PB, red
)

dr_cl_comp <- dr_umap %>% dplyr::group_by(condition, sample_id, cluster_mergeN) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n)*100) %>% 
  dplyr::rename('cluster' = 'cluster_mergeN','donor'= "sample_id" , 'group'='condition' ) %>%
  mutate(cluster = factor(cluster, levels(dr_umap$cluster_mergeN) ))
dr_cl_comp %>% dcast(donor~cluster, value.var = 'percent') %>% write.xlsx('rawdata_FIg10A.xlsx')

# CD26CD161GATE 
dr_cl_comp2 <- dr_umap %>% dplyr::group_by(condition, sample_id, gateCD161, gateCD26, TCRD) %>%
  summarise(n = n()) %>%
  mutate(percent = n/sum(n)*100) 
dr_umap %>% count(condition, sample_id, gateCD26, gateCD161)

dr_umap %>% dplyr::count( sample_id,gateCD161, gateCD26,TCRD, name = 'no') %>%
  group_by(sample_id) %>% summarize(VD2 = n())
  
C26etc <- dr_umap %>% group_by( sample_id) %>% 
  summarise(VD2 = median(Vd2), VD1 = median(Vd1), CD26 = median(CD26), n= n())


C26etc <- dr_umap %>% group_by( sample_id) %>% 
  summarise(VD2pos = length(TCRD[TCRD == 'TRDV2'])/n()*100,
            CD26pos = length(gateCD26[gateCD26 == 'CD26pos'])/n()*100,
           CD161pos = length(gateCD161[gateCD161 == 'CD161pos'])/n()*100,
           CCR6pos = length(gateCCR6[gateCCR6 == 'CCR6pos'])/n()*100)


ggplot(C26etc, aes(x = CD26pos, y = VD2pos, color = sample_id))+geom_point_rast()+
  scale_color_manual(values = donorcolor)


ggplot(C26etc, aes(x = CCR6pos, y = VD2pos, color = sample_id))+geom_point_rast()+
  scale_color_manual(values = donorcolor)



FigS8_donorFACS <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                               color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +
  theme(legend.key.size = unit(3, 'mm'), legend.position = 'bottom')+mytheme+notick+
  guides(color = guide_legend(title = 'donor', override.aes = list(size = 1.5), ncol = 2)) +
  scale_color_manual(values = donorcolor)
FigS8_donorFACS
FigS8_donorFACS_nolegend <- ggplot(dr_umap,  aes_string(x = 'UMAP1', y = "UMAP2",
                                                        color = "sample_id")) +
  geom_point_rast(size = 0.2, raster.dpi = 300) +  theme_classic() +
  theme(legend.position = 'none')+mytheme+notick+
  scale_color_manual(values = donorcolor)

FigS8_donorFACS_nolegend


FigS8_comp<-  ggplot(dr_cl_comp, aes(x = cluster, y = percent, fill = factor(group))) + 
  geom_boxplot(  alpha = 0.4 , color = 'black', size = 0.3,
                 width = 0.5,
                 outlier.alpha = 0)+ylab('% of donor')+xlab(NULL)+
  scale_fill_manual(values = c("#00ff99","#9933ff" ))+
  geom_jitter(aes(color = factor(donor)),position = position_dodge(width = 0.6), size = 1) + 
  scale_color_manual(values = donorcolor2) +  
  theme_minimal()+mytheme+
  theme(axis.text.x = element_text(angle = 320, hjust = 0), legend.key.size = unit(4, 'mm'),
        legend.position = c(0.8,0.9)
  )+  guides(color = F,
             fill =  guide_legend(title = "group", title.position = 'left'))


FigS8_comp



donor_compostion <-
  ggplot(dr_cl_comp, aes(x = cluster, y = percent, fill = factor(group))) + 
  geom_boxplot(  alpha = 0.4 , color = 'black', size = 0.3,
                 width = 0.5,
                 outlier.alpha = 0)+ylab('% of donor')+xlab(NULL)+
  scale_fill_manual(values = c("#00ff99","#ff66ff" ))+
  geom_jitter(aes(color = factor(donor)),position = position_dodge(width = 0.6), size = 1) + 
  scale_color_manual(values = donorcolor2) +  
  theme_minimal()+mytheme+
  theme(axis.text.x = element_text(angle = 320, hjust = 0), 
        legend.key.size = unit(4, 'mm'),
        axis.line = element_blank()
  )+  guides(color = guide_legend(override.aes = list(size = 2), title = 'donor'),
             fill = guide_legend(title = 'group'))


donor_compostion

FigS8A
FigS8B_TRDFACS_legendinner <- ggplot(dr_umap,  aes_string(x = 'Vd1', y = "Vd2",
                                                         color = "TCRD")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = c(0.65,0.8), legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))

FigS8B_TRDFACS_legendinner


FigS8B_TRDFACS <- ggplot(dr_umap,  aes_string(x = 'Vd1', y = "Vd2",
                                                          color = "TCRD")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
gglp('b',5)+
    guides(color = guide_legend( title = NULL, override.aes = list(size = 2), nrow =1)) +
  scale_color_manual(values = c("#FFFF33", "#377EB8","#E41A1C"))

FigS8B_TRDFACS

ggplot(dr_umap,  aes_string(x = 'Vd1', y = "Vd2",
                            color = "cluster_merge")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = c(0.65,0.8), legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = cl_color) + facet_wrap(~cluster_merge)



CCR6vCD26 <- ggplot(dr_umap,  aes_string(x = 'CCR6', y = "CD26",
                                         color = "cluster_mergeN")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ 
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')
CCR6vCD26

CCR4vCD127 <- ggplot(dr_umap,  aes_string(x = 'CD127', y = "CCR4",
                            color = "cluster_mergeN")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')
CCR4vCD127

CD16vCD127 <- ggplot(dr_umap,  aes_string(x = 'CD127', y = "CD16",
                            color = "cluster_mergeN")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ylim(5,10)+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')

CD16vCD127

quantile(dr_umap$CCR4, c(0, 0.001, 0.01, 0.95, 0.99, 1))


PD1vCD16 <- ggplot(dr_umap,  aes_string(x = 'PD1', y = "CD16",
                                        color = "cluster_mergeN")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ylim(5,10)+xlim(5.2, 8.2)+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')
PD1vCD16




CCR7vsCD127<-ggplot(dr_umap,  aes_string(x = 'CD127', y = "CCR7",
                                         color = "cluster_mergeN")) +
  geom_point_rast(size = 0.5, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+notick+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')

CCR7vsCD127







CD161vCD94 <- ggplot(dr_umap,  aes_string(x = 'CD161', y = "CD94",
                            color = "cluster_mergeN")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+xlim(5,8.5)+ylim(5,9.5)+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ notick+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_brewer(palette = 'Paired')
CD161vCD94


CD161vCD26 <- ggplot(dr_umap,  aes_string(x = 'CD161', y = "CD26",
                                          color = "cluster_mergeN")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'bottom', legend.key.size = unit(3, 'mm'))+ xlim(5,8.5)+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 2), nrow = 2)) +
  scale_color_brewer(palette = 'Paired')

CD161vCD26


lg <-get_legend(CD161vCD26)




lg



ggplot(dr_umap,  aes_string(x = 'cluster_merge', y = 'CCR7')) +
  geom_violin(aes(fill = alpha(cluster_merge, 0.1)))

FigS8D_FACS <- PG(list(PG(list(CD161vCD94,
                               CCR6vCD26,
                               CCR4vCD127,
                               CD16vCD127,
                               PD1vCD16,
                               CCR7vsCD127), nrow = 1),
                       lg
), nrow = 2, rh = c(1,0.3), labels = 'c')
FigS8D_FACS





FigS8AB <- PG(list(donor_compostion,FigS8B_TRDFACS), labels = c('a', "b"),
               nrow = 1  ,   axis = 'bl', align = 'hv', rw = c(1.8,1))

FigS8AB

 


Vln_allmarkers <- map(c(lineage_markers, functional_markers) %>% as.list(), function(x) {
  ggplot(dr_umap,  aes_string(x = 'cluster_mergeN', y = x, fill = "cluster_mergeN")) +
    geom_boxplot_jitter(  aes_string(outlier.color =  'cluster_mergeN'), size = 0.2, 
                   width = 0.5,raster.dpi = 300,outlier.jitter.width = 0.1,
                   outlier.alpha = 0.1)+xlab(NULL)+
    geom_violin(alpha = 0.4) +  theme_minimal() +
    scale_fill_brewer(palette = 'Paired')+
  theme(legend.position = 'none', axis.text.x = element_blank()
            )+mytheme
}
) %>% set_names(c(lineage_markers, functional_markers)) %>% as.list()


Vln_allmarkers <- map(c(lineage_markers, functional_markers) %>% as.list(), function(x) {
  ggplot(dr_umap,  aes_string(x = 'cluster_mergeN', y = x, fill = "cluster_mergeN")) +
    geom_point_rast(aes_string(color = "cluster_mergeN"), alpha = 0.1, raster.dpi = 300,
                    position=position_jitter(0.2))+xlab(NULL)+
    geom_boxplot_jitter(  size = 0.5, alpha = 0.4,outlier.color = 'transparent',
                          width = 0.6,raster.dpi = 300)+
    geom_violin( alpha = 0.4) +  theme_minimal() +
    scale_fill_brewer(palette = 'Paired')+
    scale_color_brewer(palette = 'Paired')+
    theme(legend.position = 'none', axis.text.x = element_blank()
    )+mytheme+
    theme(axis.line = element_blank())
}
) %>% set_names(c(lineage_markers, functional_markers)) %>% as.list()

Vln_allmarkers$lg <- get_legend(CD161vCD26+
                                  guides(color = guide_legend(title = NULL, ncol = 1, override.aes = list(size = 1))))

