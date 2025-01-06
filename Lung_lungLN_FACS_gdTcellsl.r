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
fcsfiles <- list.files(pattern = c('gd.fcs') )
file_name <- list.files(pattern = c('gd.fcs') )
fcsfiles
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


FCS_raw_gdT
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
# FJComp-BV421-A	 
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


FCS_raw_gdT@frames[[file_name[1]]]@exprs[,7:32]

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


gdT_corrected %>%  colnames()

Feature_rast(gdT_corrected, d1 = "CD103", d2 = "GATA3", facets = "tissue", g= "batch", noaxis = F, axis.number = T)

Feature_rast(gdT_corrected, d1 = "CD103", d2 = "CD25", facets = "tissue", g= "batch", noaxis = F, axis.number = T)

TotalT_corrected$batch
TotalT_corrected$tissue <-  TotalT_corrected$condition

Feature_rast(TotalT_corrected, d1 = "CD4", d2 = "GATA3", facets = "tissue", g= "CD8_scaled", noaxis = F,
              colorgrd = "grd2", sz = 0.2,
             axis.number = T)



Feature_rast(TotalT_corrected, d1 = "Eomes", d2 = "GATA3", facets = "tissue", g= "batch", noaxis = F,
             colorgrd = "grd2", sz = 0.2,
             axis.number = T)


ggplot(gdT_corrected, aes(x = tissue, y = ""))


colnames(gdT_corrected)

Feature_rast(gdT_corrected, d1 = 'AREG', d2 = "tissue", g= 'AREG' , noaxis = F)


ggplot(gdT_corrected, aes(x = sample, y = AREG_uncorrected))+geom_boxplot()
ggplot(TotalT_corrected, aes(x = condition , y = AREG_uncorrected))+geom_boxplot()


ggplot(gdT_corrected, aes(x = AREG, color = sample))+ 
  geom_histogram( fill="white", binwidth=0.05, position="identity",alpha=0.3)


TotalT_corrected




# # effect phenos  --------------------------------------------------------

Feature_rast(TotalT_corrected, d1 = "TCRgd_uncorrected", d2 = "CD103", g= 'batch', 
             facets = "tissue",
             noaxis = F, axis.number = T) 


Feature_rast(TotalT_corrected, d1 = "CD45RA", d2 = "CD27",g= 'batch', 
             facets = "tissue",
             noaxis = F, axis.number = T) 



Feature_rast(gdT_corrected, d1 = "Vd2", d2 = "Vg9", g= 'batch', 
             facets = "tissue",
             noaxis = F, axis.number = T) +
  geom_hline(yintercept = 2.5)+geom_vline(xintercept = 3.2)



Feature_rast(gdT_corrected, d1 = "CD103", d2 = "CD49a", g= 'batch', 
             facets = "tissue",
             noaxis = F, axis.number = T) +geom_vline(xintercept = 3.5)


Feature_rast(TotalT_corrected, d1 = "CD103", d2 = "CD49a", g= 'batch', 
             facets = "tissue", sz = 0.1,
             noaxis = F, axis.number = T) +geom_vline(xintercept = 3.2)


gdT_corrected


vroom::vroom_write(gdT_corrected, "gdT_corrected_FACS.csv")

vroom::vroom_write(TotalT_corrected, "TotalT_corrected_FACS.csv")



gdT_corrected  %<>% mutate(gdTtype= case_when(  Vd2 > 3.2 & Vg9 > 2.5 ~ "Vg9Vd2", TRUE ~ "NonVd2" ),
 TRM = case_when(CD103 > 3.2 ~ "TRM", TRUE ~ "Tcirc" )
  )
TotalT_corrected %<>% mutate(
                             TRM = case_when(CD103 > 3.2 ~ "TRM", TRUE ~ "Tcirc" )
)



Feature_rast(TotalT_corrected, d1 = "CD4", d2 = "CD8", g= 'TRM', 
             facets = "tissue",
             noaxis = F, axis.number = T) +
  geom_hline(yintercept = 2.5)+geom_vline(xintercept = 3.2)


Feature_rast(TotalT_corrected, d1 = "TCRab", d2 = "TCRgd", g= 'TRM', 
             facets = "tissue",
             noaxis = F, axis.number = T) 





Feature_rast(gdT_corrected, d1 = "CD103", d2 = "CD49a", g= 'TRM', 
             facets = "tissue",
             noaxis = F, axis.number = T) +geom_vline(xintercept = 3.2)



Feature_rast(as_tibble(gdT_corrected), d1 = "gdTtype", d2 = "CCR6", g= 'TRM', 
             facets = "tissue",
             noaxis = F, axis.number = T) 


table(gdT_corrected$tissue, gdT_corrected$TRM) 

colnames(TotalT_corrected)

Feature_rast(as_tibble(gdT_corrected), d1 = "TRM", d2 = "GATA3", g= 'batch', 
             facets = "tissue",
             noaxis = F, axis.number = T) 



  map(c("GATA3", "Eomes", "CD357", "AREG", "GMCSF", "GzmA") ,~
        ggplot(gdT_corrected, aes_string(x = .x, fill = "TRM"))+geom_density(alpha=0.4)+facet_wrap(~tissue)) %>% PG(ncol = 3, labels = 'gdT cells')

  
  

PG(list(       
ggplot(gdT_corrected, aes(x=CD26, fill = gdTtype)) + geom_density(alpha=0.4) + fill_m(), 

ggplot(gdT_corrected, aes(x=KLRG1, fill = gdTtype)) + geom_density(alpha=0.4) + fill_m(),
ggplot(gdT_corrected, aes(x=CD103, fill = gdTtype)) + geom_density(alpha=0.4) + fill_m(),
ggplot(gdT_corrected, aes(x=CD16, fill = gdTtype)) + geom_density(alpha=0.4) + fill_m()



), ncol = 3)


map(c("GATA3", "Eomes", "CD357", "AREG", "GMCSF", "GzmA", "CD16") ,~
      ggplot(TotalT_corrected, aes_string(x = .x, fill = "TRM"))+geom_density(alpha=0.4)+facet_wrap(~tissue)) %>% PG(ncol = 3, labels = "abT cells")

exprs(FCS_raw_gdT)



gdMFI <-  
  fsApply(FCS_raw_gdT , Biobase::exprs)[,7:32] %>% as_tibble()%>% `colnames<-`(paste0(all_markers, "_MFI"))


totalMFI <- 
  fsApply(FCS_raw_totalT , Biobase::exprs)[,7:32] %>% as_tibble()%>% `colnames<-`(paste0(all_markers, "_MFI"))




gdT_corrected  %<>% cbind(gdMFI)

# TotalT_corrected  %<>%  cbind(totalMFI) %T>% vroom::vroom_write("TotalT_corrected_FACS.csv")
# 
# vroom::vroom_write(gdT_corrected, "gdT_corrected_FACS.csv")

