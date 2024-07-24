# the analysis script for Total T  cells from lung and luLN
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




TotalT_corrected <- readRDS('totalTcell_processed_dataframe.rds')




# ##metadata tabl ---------------------------------------------------------



setwd('/home/big/tanlikai/Lung/FACSdata/')


# setwd('/home/big/tanlikai/Lung/PBMCs_Lung_panel_test/lung_LN/')
dir.create('figs')
fcsfiles <- list.files(pattern = c('TotalT', 'fcs') )
file_name <- list.files(pattern = c('TotalT', 'fcs') )

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

read.flow

###input trimmed FCS file
FCS_raw_totalT <- read.flowSet(
  files = md$file_name,
  transformation = FALSE,
  truncate_max_range = FALSE)



for (i in file_name) {
  print(i)
  dim(FCS_raw_totalT@frames[[i]]) %>%  print()
}

# add barcode to every cell
for (i in file_name) {
  print(i)
  FCS_raw_totalT@frames[[i]]
}



# setwd('..')
# colQuantiles(FCS_humangdT@frames[['export_DMHH190325 IELps_trialUnmixedSamplesKO_Live cells.fcs']]@exprs,
#              probs = c(0.01, 0.2, 0.5, 0.75, 0.99))


FCS_raw_totalT@frames[["TotalT_LN1_p254_T4047_Live_T_cells.fcs"]]@parameters@data



##panels
panel <- FCS_raw_totalT@frames[[file_name[1]]]@parameters@data[, 1:2 ]
nrow(panel)
panel
# panel_fcs <- parameters(FCS_raw_totalT[[1]]) %>% Biobase::pData()
panel$desc <- gsub("-","",panel$desc)
Channel <- colnames(FCS_raw_totalT)

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
lineage_markers <- setdiff(all_markers, c('Live','CD3', 'CD45RO', 'AREG', 'Gzma', 'GMCSF'))
lineage_markers_2 <-c('Vg9', 'Vd2', 'Vd1' ,'CD45RA', 'CD16', 'CD4', 'CD8', 'CD49a','KLRG1', 'GATA3',
                      'CD103', 'CD45RO', 'CD357', 'CD26', 'Eomes','TCRab', 'CTLA4', 'CD25', 'TCRgd')
# # Functional markers  we don't use it here but should be useful for furture experiment design
functional_markers <-c('Vd1','Vd2','CD3',"TCRgd", 'TCRab', 'AREG', 'Gzma', 'GMCSF')

#abs the expr matrix

# for (i in file_name) {
#   FCS_humangdT@frames[[i]]@exprs <- abs(FCS_humangdT@frames[[i]]@exprs)
# }


# ##biexpontential transform of data --------------------------------------

file_name[1]
##############check rawdata distribution
rawdata <- FCS_raw_totalT@frames[[file_name[1]]]@exprs[,7:32] %>% asinh()
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

figsave(plot_grid(plotlist = raw, ncol = 6), "rawsignals_Lung_totalT.pdf", w = 300, h = 300)

# save_plot("raw_Lung.pdf", plot_grid(plotlist = raw, ncol = 6), base_height = 20, base_width = 20)

# figsave(plot_grid(plotlist = raw, ncol = 6), 'rawsignals_lung.pdf', 300, 300)


# data transformation 
fcs_totalT_transform_pre_qc <- fsApply(FCS_raw_totalT, function(x) {transform(x, estimateLogicle(x, c(Channel[7:32])))})

# save the transformed data
write.flowSet(fcs_totalT_transform_pre_qc, outdir = file.path(getwd(), "Transformed_FCS_files") , filename = paste0("",FCS_raw_totalT@phenoData@data$name)) 

fcs_totalT_transform_pre_qc <- read.flowSet(files = md$file_name, path=file.path(getwd(), "Transformed_FCS_files"), transformation = FALSE, truncate_max_range = FALSE)



# automatic QC,  ------------------------------------------------



# Automated quality control ,
# 
dir.create('PeacoQCresults')


  
for(i in 1:length(sampleNames(fcs_totalT_transform_pre_qc))){
  ff <-fcs_totalT_transform_pre_qc[[i]]
  channels=Channel[7:32]
  peacoqc_res <- PeacoQC(ff, Channel[7:32], determine_good_cells = "all",
                              save_fcs = TRUE, plot=TRUE, output_directory = "PeacoQCresults")
} 


fcs_totalT_transform <- read.flowSet(pattern = 'TotalT',  path=file.path(getwd(), "PeacoQCresults","PeacoQC_results","fcs_files"), transformation = FALSE, truncate_max_range = F)

# 



# QC result
sample_ids_raw <- rep(file_name, fsApply(fcs_totalT_transform_pre_qc, nrow))

sample_ids_trans <- rep(file_name, fsApply(fcs_totalT_transform , nrow))
cell_table_raw <- table(sample_ids_raw)
cell_table_trans <- table(sample_ids_trans)


qc_table <- data.frame(cell_table_raw,cell_table_trans)
qc_table$sample_ids_trans <- NULL
colnames(qc_table) <- c("Sample_ID","Pre-QC cell count","Post-QC cell count")
qc_table$Removed <- qc_table$`Pre-QC cell count`-qc_table$`Post-QC cell count`
qc_table$Removed_perc <- round(((qc_table$Removed/qc_table$`Pre-QC cell count`)*100), digits=2)
colnames(qc_table) <- c("Sample_ID","Pre-QC cell count","Post-QC cell count","Removed (n)","Removed (%)")
qc_table  %>% kbl(row.names = T) %>% kable_styling(full_width = F, position="left")


fcs_totalT_transform <- fcs_totalT_transform[sampleNames(fcs_totalT_transform), colnames(fcs_totalT_transform)[7:32]]


# manually exponential transform, skip  -----------------------------------




#the a value will decide how much the negative values will be put together
testraw <-FCS_raw_totalT@frames[[file_name[1]]]
file_name[[1]]


biexp <- biexponentialTransform('biexp_transform',w = 5,
                                a= 200000)
biexp2 <- biexponentialTransform('biexp_transform',w =5, 
                                a= 100 , c = 300)


after.1<- c()
  after.1 <- transform(testraw, transformList(panel[c(7:32), 1], biexp))
colnames(after.1) <- gsub("-| ", "_", colnames(after.1))


# check if transfroming is enough
channales <- colnames(after.1@exprs)[c(7:32)] %>% as.character()
bie_tr <- c()

for (i in 1:length(channales)) {
  bie_tr[[i]] <- ggplot(after.1@exprs %>% as.data.frame(), aes_string(x = channales[i])) +
    geom_density()+ggtitle(all_markers[i] )
}
plot_grid(plotlist = bie_tr, ncol = 6)



##do the formal transform
##biexp
FCS_humangdT <- FCS_raw
for (i in file_name) {
  FCS_humangdT@frames[[i]] <- transform(FCS_raw@frames[[i]], 
                                             transformList(panel[c(7:32),1], biexp))
}

# for (i in file_name) {
#   FCS_humangdT@frames[[i]] <- transform(FCS_humangdT@frames[[i]], 
#                                         transformList(panel[c(7:18),1], biexp))
# }


FCS_humangdT_biexp <- fsApply(FCS_humangdT, function(x){
  colnames(x) <- colnames(FCS_humangdT)
  expr <- Biobase::exprs(x)
  # expr <- expr[,all_markers ]
  exprs(x) <- expr
  x
})





# data scale, skip ----------------------------------------------------------------------





expr <- fsApply(fcs_totalT_transform , Biobase::exprs)

colnames(expr)

# colQuantiles(fcs_raw@frames$CD3neg.fcs@exprs[, 21], probs = c(0.01, 0.25, 0.5, 0.75, 0.99))
## Extract expression
# expr <- fsApply(FCS_humangdT_biexp, Biobase::exprs)
expr <- expr%>% `colnames<-`(all_markers)
dim(expr)

colnames(expr)


###data normalization

# testdf <- data.frame(HA = c(1:100), GD = c(2:101), C = (-1:98), D = (20:119)) %>% as.matrix()
# colQuantiles(testdf, probs = c(0.01, 0.50))
#determin the datarange  of each marker
rng <- colQuantiles(expr, probs = c(0.01, 0.99))

##normalize the expression matrix by :
#expression value - minimal expression value across all cells / data range
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

dim(expr01)

##scale 


### Generate sample IDs corresponding to each cell in the `expr` matrix

sample_ids <- rep(md$file_name, fsApply(fcs_totalT_transform , nrow))

# donor <- rep(md$donor, fsApply(FCS_humangdT, nrow))

# stimulation <- rep(md$stimulation, fsApply(FCS_humangdT, nrow))

tissue <- rep(md$tissue, fsApply(fcs_totalT_transform , nrow))
donor <- rep(md$donor, fsApply(fcs_totalT_transform , nrow)) 
ID <-  rep(md$ID, fsApply(fcs_totalT_transform , nrow)) 

## Diagnostic plots
ggdf <- data.frame(sample_id = sample_ids, tissue, donor, ID, expr)
 ggdf %>% head
##transfor the ggdf to long table
ggdf <- melt(ggdf, id.var = c('sample_id','tissue',
                              'donor', 'ID'
                              ),
             value.name = "expression", variable.name = "antigen")
ggdf %>% head
# mm <- match(ggdf$sample_id, md$sample_id)
# ggdf$condition <- md$condition[mm]

dim(ggdf)
##downsample


sampFreq<-function(cdf,col,ns) {
  x<-as.factor(cdf[,col])
  freq_x<-table(x)
  prob_x<-freq_x/sum(freq_x)
  df_prob = prob_x[as.factor(cdf[,col])]
  nr=nrow(cdf)
  sLevels = levels(as.factor(cdf[,col]))
  nLevels = length(sLevels)
  rat = ns/nr
  rdata = NULL
  for (is in seq(1,nLevels)) {
    ldata <- cdf[cdf[,col]==sLevels[is],]
    ndata <- nrow(ldata)
    nsdata = max(ndata*rat,1)
    srows <- sample(seq(1,ndata),nsdata,replace=rat>1)
    sdata <- ldata[srows,]
    rdata <- rbind(rdata,sdata)
  }
  return(rdata)
}

#downsample cell number to ~300,000
ggdf_ds <- sampFreq(ggdf, col = 1, ns =  20000)
head(ggdf_ds)
nrow(ggdf_ds)


table(ggdf_ds$sample_id)

ggdf_ds

histog_allmarkers <-ggplot(ggdf_ds, aes(x = expression, color = ID,
                          group = ID)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  # theme_classic() + scale_color_manual(labels = c('Neonate', 'Adult'), values = c('blue', 'red'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 8),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) 

histog_allmarkers %T>%  figsave('hist_permeabilization_without_batch_correction_lung_totalT_foursample.pdf', 300, 300)




# boxplot <-ggplot(ggdf, aes(y = expression, color = sample_id,
#                           x = antigen))+ geom_boxplot()
# boxplot


# Get the median marker expression per sample
# expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
#   group_by(sample_id) %>%
#   summarize_all(funs(median))
# #transform to long-table
# expr_median_sample <- t(expr_median_sample_tbl[, -1])
# colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

##multi-dimensional scaling (MDS) plotSuch plots show similarities between samples measured in an unsupervised way and give a sense of how much differential expression can be detected before conducting any formal tests.
# expr_median_sample <- as.data.frame(expr_median_sample)
# expr_median_sample$Condition <- 0
#
# mds <- plotMDS(expr_median_sample, plot = FALSE)
# ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
#                    sample_id = colnames(expr_median_sample))
# mm <- match(ggdf$sample_id, md$sample_id)
# ggplot(ggdf, aes(x = MDS1, y = MDS2)) +
#   geom_point(size = 2, alpha = 0.8) +
#   geom_label_repel(aes(label = sample_id)) +
#   theme_bw()

# Batch effection correction ----------------------------------------------




library(cyCombine)


write.csv(md %>%  mutate(file_name = str_replace(file_name, 'cells', 'cells_QC')), 'PeacoQCresults/PeacoQC_results/fcs_files/metadata_TotalT_lung_LN.csv')
md

uncorrected  <- prepare_data(data_dir = '/home/big/tanlikai/Lung/FACSdata/PeacoQCresults/PeacoQC_results/fcs_files/',
                          
                            markers = all_markers,
                            transform = FALSE,
                            pattern = "fcs",      
                            metadata = 'PeacoQCresults/PeacoQC_results/fcs_files/metadata_TotalT_lung_LN.csv', # Can also be .csv file or data.frame object
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
dir.create('figs/batch')

all_markers
uncorrected %>%
  detect_batch_effect(markers = all_markers,
                      batch_col = 'batch',
                      out_dir = 'figs/batch', 
                      seed = 434,
                      name = 'Tcells')

colnames(uncorrected)


map(all_markers, ~  ggplot(uncorrected, aes_string(x = ., color = 'batch',
                                                 group = 'batch')) +
      geom_density() ) %>% PG(nrow = 5) %T>% figsave('batch/batch_uncorrected_histgram.pdf', 600, 400)

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
      geom_density() ) %>% PG(nrow = 5) %T>% figsave('batch/batch_corrected_histgram.pdf', 600, 400)


# Set column for evaluation of EMD (per-cluster)
celltype_col <- "som"
colnames(corrected)
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
plot_density(uncorrected, corrected, ncol = 6)  %T>% figsave('batch/histogram_corrected_vs_uncorrected.pdf', 400, 400)

plot1 <- plot_dimred(uncorrected, name = 'Uncorrected', type = 'umap', markers = lineage_markers)
plot2 <- plot_dimred(corrected, name = 'Corrected 8x8', type = 'umap', markers = lineage_markers)




plot3 <- plot_dimred(corrected, name = 'corrected', type = 'umap', plot = 'som',
                     markers = lineage_markers_2)





(rasterise(plot3,dpi = 300) +facet_wrap(~Batch))


cowplot::plot_grid(  rasterise(plot1,dpi = 300) ,  rasterise(plot2,dpi = 300)) %>% 
  figsave('batch/Uncorrected_vs_batchcorrected.pdf', 400, 200)
plot1$data

(rasterise(plot1,dpi = 300) +facet_wrap(~Batch)) %T>% 
  figsave('batch/Uncorrected_facet4.pdf', 400, 200)
(rasterise(plot2,dpi = 300) +facet_wrap(~Batch))%T>% 
  figsave('batch/corrected_facet4.pdf', 400, 200)
(rasterise(plot3,dpi = 300) +facet_wrap(~Batch))%T>% 
  figsave('batch/corrected_facet4_fewermarker.pdf', 400, 200)

colnames(uncorrected)
colnames(corrected)

plot3$data

head(corrected)

ncol(uncorrected)
head(uncorrected)[,3:27]


rng <- colQuantiles(as.matrix(corrected[3:28]), probs = c(0.005, 0.995))

##normalize the expression matrix by :
#expression value - minimal expression value across all cells / data range
expr01 <- t((t(as.matrix(corrected[3:28])) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1




TotalT_corrected <- cbind(plot3$data[,1:2], corrected, 
                          rename_with(uncorrected[,3:28] , ~paste0(., '_uncorrected')),
                          rename_with(as.data.frame(expr01) , ~paste0(., '_scaled'))
                          ) %>%
  mutate(sample = paste0(batch, '_','condition')) %>% 
  as.data.frame()

head(TotalT_corrected)

# TotalT_corrected$sample <- paste0(TotalT_corrected$batch, '_', TotalT_corrected$condition)

Feature_rast(TotalT_corrected, 'condition', d1 = 'UMAP1', d2 = 'UMAP2', facets  = 'batch', do.label = F) %T>% 
  figsave('four_donor_tissue_UMAP.pdf', 200, 200)

Feature_rast(TotalT_corrected, 'sample', d1 = 'UMAP1', d2 = 'UMAP2',  do.label = T) %T>% 
  figsave('allsamples.pdf', 200, 150)

Feature_rast(TotalT_corrected, all_markers, d1 = 'UMAP1', d2 = 'UMAP2', ncol = 4, sz = 0.1) %T>% 
  figsave('allmarkers.pdf', 300, 400)


Feature_rast(TotalT_corrected, paste0(all_markers, '_scaled'), d1 = 'UMAP1', d2 = 'UMAP2', ncol = 4, sz = 0.1) %T>% 
  figsave('allmarkers_scaled.pdf', 300, 400)

saveRDS(TotalT_corrected, 'totalTcell.rds')


# SOM clustering ----------------------------------------------------------

fsom <- ReadInput(fcs_totalT_transform , transform = FALSE, scale = F)
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

cdata <- TotalT_corrected[, all_markers] %>%  `colnames<-`(abcolors) %>% as.matrix()

fsom$data <- cdata 

# build up som 
# color to use 

lineage_colors <- panel %>%  dplyr::filter(desc %in% lineage_markers_2) %>%  pull(name) %>% as.vector()



set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_colors)
head(som)
# Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)
nmc <- 30
codes <- som$map$codes
t(codes) %>% head
dir.create('som_corrected')

dir.create('som_corrected/consensus_plots')

plot_outdir <- "som_corrected/consensus_plots"


mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 1000,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)
head(mc)

code_clustering1 <- mc[[17]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]



# SOM with out CD4CD8 -----------------------------------------------------

lineage_markers_3 <- c(setdiff(lineage_markers_2, c('CD4', 'CD8')), 'TCRab', 'TCRgd')

lineage_colors_w4w8 <- panel %>%  dplyr::filter(desc %in% lineage_markers_3) %>%  pull(name) %>% as.vector()

set.seed(1234)
som2 <- BuildSOM(fsom, colsToUse = lineage_colors_w4w8)


library(ConsensusClusterPlus)
nmc <- 30
codes <- som2$map$codes
t(codes) %>% head
dir.create('som_corrected_w4w8')

dir.create('som_corrected_w4w8/consensus_plots')

plot_outdir <- "som_corrected_w4w8/consensus_plots"


mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 1000,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)
head(mc)

code_clustering2 <- mc[[17]]$consensusClass
cell_clustering2 <- code_clustering2[som2$map$mapping[,1]]


# # plot_outdir <- "FCS/exp1125/consensus_plots"
# # codes_ns <- som_ns$map$codes
# # 
# # mc <- ConsensusClusterPlus(t(codes_ns), maxK = nmc, reps = 1000,
# #                            pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
# #                            clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
# #                            distance = "euclidean", seed = 1234)
# # 
# # 
# # ## Get cluster ids for each cell
# # # choose 10 clusters according to the plot above
# # 
# # 
# # code_clustering2 <- mc[[15]]$consensusClass
# # cell_clustering2 <- code_clustering2[som_ns$map$mapping[,1]]
# 
# 
# length(cell_clustering1)
# # !!skip ##clustering by hclust --------------------------------------------------
# #complete
# #scale expression
# nrow(codes)
# 
# 
# celldist <- dist(codes,  method = 'euclidean')
# cell_hc <- hclust(celldist, method = 'complete')
# 
# plot(cell_hc)
# abline(h = 3.5, col = "red")
# abline(h = 3.9, col = "red")
# 
# h1 <- cutree(cell_hc, h = 3.5) 
# h2 <- cutree(cell_hc, h = 3.9) 
# h3 <- cutree(cell_hc, k = 15) 
# 
# 
# som$map$mapping
# 
# cell_hcluster <- h1[som$map$mapping[,1]]
# cell_hcluster_2 <- h2[som$map$mapping[,1]]
# cell_hcluster_3 <- h3[som$map$mapping[,1]]
# 
# length(cell_hcluster)
# ###assign colors
# color_clusters <- set_sample(ggplotColours(30))
# 
# #average or median
# cell_hca <- hclust(celldist, method = 'average')
# plot(cell_hca)
# abline(h = 2.8, col = "red")
# 
# h4 <- cutree(cell_hca, k = 15) 
# cell_hcluster_4 <- h4[som$map$mapping[,1]]
# 


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
heatmap_orig_hcluster <- plot_clustering_heatmap_wrapper(expr = corrected[,3:27],
                                           expr01 = expr01[,1:25],
                                           cell_clustering = cell_clustering1,
                                           color_clusters = ggplotColours(20))%>% ggplotify::as.ggplot()


TotalT_corrected$som_cluster <- as.factor(cell_clustering1)
TotalT_corrected$som_cluster_w4w8 <- as.factor(cell_clustering2)


Feature_rast(TotalT_corrected, c('som_cluster','som_cluster_w4w8'),d1 = 'UMAP1', d2 = 'UMAP2', facets = 'batch', sz = 0.2)

Feature_rast(TotalT_corrected, c('som_cluster','som_cluster_w4w8'),d1 = 'cUMAP_1', d2 = 'cUMAP_2', facets = 'batch', sz = 0.2)




saveRDS(TotalT_corrected, 'totalTcell_processed_dataframe.rds')


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


# Umap_withoutCD4CD8 --------------------------------------------------------------------
library(reticulate)

#Run tSNE
#position info stored in the tsne_out$Y
set.seed(1234)
# tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE )
tsne_expr

all_markers
lineage_markers_3 <- c(setdiff(lineage_markers_2, c('CD4', 'CD8')), 'TCRab', 'TCRgd')




umap_out_woCD4CD8 <- uwot::umap(TotalT_corrected[,lineage_markers_3], scale = "Z", n_neighbors = 10, ret_model = T, min_dist = 0.1)


# uwot::save_uwot(umap_out3, paste0(dir_FACS,'umapmodel'))

plot(umap_out_woCD4CD8$embedding)


#introduce UMAP without CD4 CD8 into datafrom----------------------------------------------


# umap_out$layout %>% head
TotalT_corrected <- data.frame(cUMAP_1 = umap_out_woCD4CD8$embedding[,1] , 
                      cUMAP_2 = umap_out_woCD4CD8$embedding[,2] ,
                      TotalT_corrected)


TotalT_corrected %>%  colnames()

Feature_rast(TotalT_corrected,c('CD4_scaled', 'CD8_scaled', 'KLRG1_scaled',"CD45RA_scaled","CD16_scaled",
                                'GzmA_scaled', 'CD103_scaled','som_cluster_withoutCD4CD8'),
             d1 = 'cUMAP_1', d2 = 'cUMAP_2' )

T


# CD4 and CD8 gating ------------------------------------------------------
TotalT_corrected %>%  colnames()
TotalT_corrected$condition
Feature_rast(TotalT_corrected, 'condition', d1 = 'CD4', d2= 'CD8', noaxis = F, axis.number = T,facets = 'batch')



celldist <- dist(TotalT_corrected[,c('CD4', 'CD8')],  method = 'euclidean')
library(fastcluster)
hc <- fastcluster::hclust(celldist,method = 'complete')

plot(hc)
abline(h = 4.5, col = 'red')

h45 <- cutree(hc, h = 4.5) 
h45


TotalT_corrected$CD4CD8_cluster <- as.factor(h45)
Feature_rast(TotalT_corrected, c( 'CD4CD8_cluster'), d1 = 'CD4', d2= 'CD8', noaxis = F, axis.number = T)


CD4CD8_colors <- panel %>%  dplyr::filter(desc %in% c('CD4', 'CD8')) %>%  pull(name) %>% as.vector()

set.seed(1234)
CD4CD8som <- BuildSOM(fsom, colsToUse = CD4CD8_colors)
head(som)
# Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)
nmc <- 5
codes <- CD4CD8som$map$codes
t(codes) %>% head
dir.create('som_48')

dir.create('som_48/consensus_plots')

plot_outdir <- "som_48/consensus_plots"


mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 1000,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "pdf",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)
head(mc)

code_48 <- mc[[3]]$consensusClass
cell_48 <- code_48[CD4CD8som$map$mapping[,1]]

TotalT_corrected$CD4CD8_cluster <- as.factor(cell_48)
Feature_rast(TotalT_corrected, c( 'CD4CD8_cluster'), d1 = 'CD4', d2= 'CD8', noaxis = F, axis.number = T)



# plot(cell_hc)
# abline(h = 3.5, col = "red")
# abline(h = 3.9, col = "red")
# 
# h1 <- cutree(cell_hc, h = 3.5) 
# h2 <- cutree(cell_hc, h = 3.9) 
# h3 <- cutree(cell_hc, k = 15) 
# 
# 
# som$map$mapping
# 
# cell_hcluster <- h1[som$map$mapping[,1]]
library(ggiraph)
library(shiny)
library(plotly)

colnames(TotalT_corrected)

Feature_rast(TotalT_corrected, g= "TCRgd" ,d1 ="CD4_uncorrected", d2= "CD8_uncorrected")


TotalT_corrected$rowname <- rownames(TotalT_corrected)
TCRab_CD4CD8 <- ggplot(TotalT_corrected, aes(y = CD4, x = CD8, color = TCRgd, key = rowname,
                                    tooltip = rowname, data_id = rowname))+geom_point(size = 0.3)+
  scale_color_gradientn( 
    colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  theme_minimal()
ggplotly(TCRab_CD4CD8)


ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brush")
)
# CD8

server <- function(input, output, session) {
  frame2<<- NULL
  
  frame2 <<-data.frame()
  # frame2$number <<- row.names(TCRabgd)
  # nms <- TCRabgd$rowname
  
  output$plot <- renderPlotly({
    ggplotly(TCRab_CD4CD8) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRab_CD4CD8[frame2$key, 'CD4orCD8'] <<- 'CD8+'
  })
  
}
shinyApp(ui, server)
TCRab_CD4CD8 <- ggplot(TotalT_corrected, aes(y = CD4, x = CD8, color = TCRgd, key = rowname,
                                             tooltip = rowname, data_id = rowname))+geom_point(size = 0.3)+
  scale_color_gradientn( 
    colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))+
  theme_minimal()

ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("click"),
  verbatimTextOutput("brush")
)
# CD8
server <- function(input, output, session) {
  frame2<<- NULL
  
  frame2 <<-data.frame()
  # frame2$number <<- row.names(TCRabgd)
  # nms <- TCRabgd$rowname
  
  output$plot <- renderPlotly({
    ggplotly(TCRab_CD4CD8) %>% layout(dragmode = "lasso")
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      frame2 <<- frame2[is.null(frame2$pointNumber), ] # Optional line to remove the previous selections 
      frame2 <<- rbind(frame2, d)
      
    }
    # frame2$TCRs<<- 'TCRab+'
    TCRab_CD4CD8[frame2$key, 'CD4orCD8'] <<- 'CD8+'
  })
  
}
shinyApp(ui, server)





















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

FigS8E_Vln_all <- PG(Vln_allmarkers, nrow = 3, labels = 'd')

FigS8E_Vln_all


FigS8 <- PG(list(FigS8AB, FigS8D_FACS, FigS8E_Vln_all), nrow = 3, rh = c(1.2, 0.8, 2))+
  draw_figure_label("Figure S8", position = 'top.right')

ggsave2(FigS8, filename = 'manu_fig/SupFigS8_FACS_20200609.pdf', width = 190, height = 255, units = 'mm')

figsave(FigS8AB,'figs8ab.pdf', w = 190, h = 76.5)

FigS8

CD161vCCR6 <- ggplot(dr_umap,  aes_string(x = 'CD161', y = "CCR6",
                                          color = "cluster_merge")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ xlim(5,8.5)+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = cl_color ) 

CD161vCCR6 



 



CCR7vsCD127<-ggplot(dr_umap,  aes_string(x = 'CD127', y = "CCR7",
                                         color = "cluster_merge")) +
  geom_point_rast(size = 0.1, raster.dpi = 300)

ggplot(dr_umap %>% filter(cluster_merge == 'Vd1 Active CTL'), aes_string(x = 'Vd1', y = "Vd2", color = "TCRD") )+
  geom_point_rast(size = 0.8, raster.dpi = 300) 

dr_umap %>% filter(cluster_merge == 'Vd1 Active CTL') %>% dplyr::count(TCRD)
dr_umap %>% dplyr::count(TCRD)


ggplot(dr_umap,  aes_string(x = 'cluster_mergeN', y = "Vd2", fill = "cluster_mergeN")) +
  geom_point_rast(aes_string(color = "cluster_mergeN"), alpha = 0.1, raster.dpi = 300,
                  position=position_jitter(0.2))+
  geom_boxplot_jitter(  size = 0.5, alpha = 0.4,outlier.color = 'transparent',
                        width = 0.6,raster.dpi = 300)+
  geom_violin( alpha = 0.4) +  theme_minimal() +
  scale_fill_brewer(palette = 'Paired')+
  scale_color_brewer(palette = 'Paired')+
  theme(legend.position = 'none', axis.text.x = element_blank())



ggplot(dr_umap,  aes_string(x = 'CD161', y = "CD26",
                            color = "CD26")) +
  geom_point_rast(size = 0.1, raster.dpi = 300) +  theme_bw() +mytheme+notick+
  theme(legend.position = 'none', legend.key.size = unit(3, 'mm'))+ xlim(5,8.5)+
  guides(color = guide_legend( title = NULL, override.aes = list(size = 1.5), ncol = 1)) +
  scale_color_manual(values = cl_color ) 

dr_umap$condition
Feature_rast(dr_umap %>% filter(TCRD == 'TRDV2' & condition == 'PB'), 
             g = 'CD161', d1 = 'CD94', d2 = 'CD26', facet = 'condition',
             color_grd = 'grd', noaxis = F, axis.number = T) 

Feature_rast(dr_umap %>% filter(TCRD == 'TRDV2' & condition == 'PB'), 
             g = 'CCR6', d1 = 'CD161', d2 = 'CD26', facet = 'condition',
             color_grd = 'grd', noaxis = F, axis.number = T) 

Feature_rast(dr_umap %>% filter(TCRD == 'TRDV2' ), 
             g = 'CD26', d1 = 'CCR6', d2 = 'CD94', facet = 'condition',
             color_grd = 'grd', noaxis = F, axis.number = T) + facet_grid(~condition)


Feature_rast(dr_umap, 
             g = 'CD161', d1 = 'Vd1', d2 = 'Vd2', facet = 'condition',
             color_grd = 'grd', noaxis = F, axis.number = T) + facet_grid(~condition)



GDT_2020AUG$group
GDT_2020AUG@assays$RNA@counts

FetchData(subset(GDT_2020AUG, subset = paired == 'GV9 DV2'), vars = c('CCR6', 'KLRD1', "DPP4"), slot = 'scale.data')


Feature_rast(GDT_2020AUG, 'CCR6', facet = 'orig.ident') +facet_grid(~orig.ident)



Feature_rast(FetchData(subset(GDT_2020AUG, subset = paired == 'GV9 DV2' ), 
                       vars = c('CCR6', 'KLRD1', "DPP4"), slot = 'scale.data'), 'CCR6', d1 = 'KLRD1', d2 = 'DPP4',
             noaxis = F, color_grd = 'grd',)

#

# R3 required VD1 vs other VD ---------------------------------------------

C1_CB <- subset(GDT_2020AUG, subset = group == 'Newborn' & Cell_cluster == 'c1' &
                  v_gene_TRD %in% c('TRDV1', 'TRDV2','TRDV3', 'TRDV5', 'TRDV8'))
ClusterCompare(C1_CB, id1 = 'TRDV1', id2 = c('TRDV3', 'TRDV5','TRDV2','TRDV8'), rm = 'CD4',
               group.by = "v_gene_TRD", log2fc = 0.25)




