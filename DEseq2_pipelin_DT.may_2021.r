
# DESeq2 piple ------------------------------------------------------------


setwd('/home/big/tanlikai/Demo_projects/DT/CleanData/')

library(DESeq2)
library(ggplot2)
library(ggrastr)
library(shadowtext)
library(ggrepel)
library(purrr)
library(cowplot)
library(magrittr)
library(dplyr)
library(stringr)
library(parallel)
library(vroom)
library(tibble)
library(limma)
library(clusterProfiler)
library(openxlsx)

# install.packages('rJava')
figsave <- function (p, filename,  w =50, h = 60, device = cairo_pdf,
                     path = 'result/fig/',
                     scale = 1, units = 'mm', dpi = 300
) {
  ggsave2( plot = p, filename = filename,device = device,
           path = path, scale = scale, width = w, height = h,
           units = units, dpi = dpi, limitsize = TRUE)
}


# creat count table -------------------------------------------------------



allcsv <- list.files(path = 'result', pattern = '.csv') 


allcounts <- map(allcsv, ~ vroom(paste0('result/',.x)) %>% select(1,7,8) %>% 
                   `colnames<-`(c('GeneID', 'GeneSymbol', str_remove(.x, '.csv')))  ) %>% set_names(allcsv)

allcounts$BCMI_1.count.csv

allcounts %<>% reduce(full_join, by = c('GeneID', 'GeneSymbol'))

allcounts %>% colnames
# remove duplicated genes 
allcounts %<>% distinct(GeneSymbol, .keep_all = T)

allcounts %<>% column_to_rownames('GeneSymbol')

colnames(allcounts) <- str_remove_all(colnames(allcounts), '.count')

head(allcounts)
dim(allcounts)

rowmeans <- rowMeans(allcounts[,2:23]) %>% sort(decreasing = T)

rowsum <- rowSums(allcounts[,2:23]) %>% sort(decreasing = T)


rowmeans %>% quantile()
rowsum %>% quantile()
rowmeans[rowmeans>0]
# filter low read gene
after <- apply(allcounts[,2:17],1,function(x) length(x[x>10])>=3)   
qcn <- after[after == 1] %>% names


allcounts_f <- allcounts[qcn,2:17]

allcounts

coldata <- data_frame(sample= colnames(allcounts_f), group = str_extract(colnames(allcounts_f), '.+(?=_\\d+)')  %>% as.vector() %>% as.factor())
coldata
dds <- DESeqDataSetFromMatrix(countData = allcounts_f,
                              colData = coldata,
                              design= ~ group)

dds <- DESeq(dds)

dds@assays$mu

res <- results(dds, name="group_BCNMI_vs_BCMI")
res$Symbol <- rownames(res)
write.xlsx(res, 'BCNMI_vs_BCMI.xlsx')

res %>% as_tibble()%>% dplyr::filter(padj <0.05 & abs(log2FoldChange >=0.5)) %>% pull(Symbol)

plotMA(resLFC, ylim=c(-2,2))

summary(res)
plotCounts(dds, gene=which.min(res$padj), intgroup="group")

res
# GSEA --------------------------------------------------------------------
BCNMI_vs_BCMI<- clusterProfiler::bitr(res$Symbol, fromType="SYMBOL",
                                   toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% `colnames<-`(c('Symbol', 'ENTREZID')) %>%
  right_join(as_tibble(res), by = 'Symbol') 



BCNMI_vs_BCMI %<>% dplyr::mutate(significancy = case_when(padj >=0.05 | abs(log2FoldChange) < 1 ~ 'Unsignificant',
                                                          padj < 0.05 & log2FoldChange >=1 ~'BCNMI upreg', 
                                                          padj < 0.05 & log2FoldChange <= -1 ~'BCMI upreg')) %>% 
  dplyr::filter(!is.na(padj))
BCNMI_vs_BCMI$log10padj = -log10(BCNMI_vs_BCMI$padj)


table(BCNMI_vs_BCMI$significancy)


top20genes <- BCNMI_vs_BCMI  %>% filter(significancy != 'Unsignificant' & !grepl('^ENS|^RP|^LINC|-AS', Symbol)) %>%
  group_by(significancy) %>%  
  top_n(20, abs(log2FoldChange)) %>% pull(Symbol)
BCNMI_vs_BCMI$top20 <- NULL
BCNMI_vs_BCMI %<>% mutate(top20 = case_when(Symbol %in% top20genes ~ Symbol))

Vol_DT <- Feature_rast(BCNMI_vs_BCMI,g = "significancy",  
             d1 = "log2FoldChange", d2 = "log10padj", noaxis = F, axis.number = T, do.label = F)+
  color_m(c('#9B1B30','#2A4B7C','lightgrey'))+
  ggtitle('BCMI vs BCNMI')+ylab('-log10(adjusted p value)')+
  geom_hline(yintercept = -log10(0.05),linetype="dashed", size = 0.1)+
  geom_vline(xintercept = c(-1,1),linetype="dashed", size = 0.1)+
  geom_text_repel(aes(label=top20, fontface = 'italic'),show.legend = F, size = gs(8),
                  max.overlaps = 30,
                  segment.size = 0.3, segment.color = 'grey50',
                  bg.color = 'white', bg.r = 0.2,
                   na.rm = T)


figsave(Vol_DT, 'Volcano_DT.pdf', 120, 120)


BCNMI_vs_BCMI_heat <- BCNMI_vs_BCMI %>% filter(padj < 0.05 & !is.na(Symbol)) %>%  distinct(Symbol, .keep_all = T) %>%
  select(paste0('KO',1:3),paste0('WT',1:3),'Symbol') %>% `rownames<-`(NC_vs_Sup_top100$Symbol)

NC_vs_Sup_heat$Symbol <- NULL
pheatmap::pheatmap(NC_vs_Sup_heat,scale = 'row')



genelist <- BCNMI_vs_BCMI$log2FoldChange
names(genelist) <- BCNMI_vs_BCMI$ENTREZID
genelist <- sort(genelist, decreasing = T)
genelist

DT_GSEA_GO <- gseGO(geneList     = genelist,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    nPerm        = 1000000,
                    minGSSize    = 10,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)%>%
  clusterProfiler::simplify(by = 'p.adjust',  select_fun = min) %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

DT_GSEA_GO@result <- DT_GSEA_GO@result %>% arrange(desc(NES))
write.xlsx(DT_GSEA_GO, 'DT_GSEA_GO.xlsx')



library(msigdbr)

as.data.frame(msigdbr::msigdbr_collections())

C4_CGN <- msigdbr(category = 'C4') %>% dplyr::select(gs_name, entrez_gene)

c7 <- msigdbr(category = 'C7') %>% dplyr::select(gs_name, entrez_gene)
C6_onco_sig <- msigdbr(category = 'C6') %>% dplyr::select(gs_description, entrez_gene)
 
Hallmark <- msigdbr(category = 'H') %>% dplyr::select(gs_name, entrez_gene)

REACTOME <- msigdbr(category = 'C2', subcategory = 'CP:REACTOME') %>% dplyr::select(gs_name, entrez_gene)


C8 <- msigdbr(category = 'C8') %>% dplyr::select(gs_name, entrez_gene)

DT_GSEA_MsigDB_c7 <- GSEA(geneList = genelist, TERM2GENE=c7,  
     minGSSize    = 20,
     pvalueCutoff = 0.01, pAdjustMethod = "BH") %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
DT_GSEA_MsigDB_c7@result <- DT_GSEA_MsigDB_c7@result %>% arrange(desc(NES))
DT_GSEA_MsigDB_c7@result %>% write.xlsx('DT_GSEA_MsigDB_c7result.xlsx')
DT_GSEA_MsigDB_c7@result %>% view


DT_GSEA_MsigDB_Hallmark <- GSEA(geneList = genelist, TERM2GENE=Hallmark,  
                          minGSSize    = 20,
                          pvalueCutoff = 0.01, pAdjustMethod = "BH") %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
DT_GSEA_MsigDB_Hallmark@result %<>%  arrange(desc(NES))
DT_GSEA_MsigDB_Hallmark %>%  write.xlsx('DT_GSEA_MsigDB_Hllmarkresult.xlsx')

DT_GSEA_MsigDB_c6 <-  GSEA(geneList = genelist, TERM2GENE=C6_onco_sig,  
                           minGSSize    = 20,
                           pvalueCutoff = 0.01, pAdjustMethod = "BH") %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
DT_GSEA_MsigDB_c6@result %<>%  arrange(desc(NES))
DT_GSEA_MsigDB_c6 %>%  write.xlsx('DT_GSEA_MsigDB_C6_onco_sig.xlsx')

DT_GSEA_MsigDB_c4 <-  GSEA(geneList = genelist, TERM2GENE=C4_CGN,  
                           minGSSize    = 20,
                           pvalueCutoff = 0.01, pAdjustMethod = "BH") %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
DT_GSEA_MsigDB_c4@result %<>%  arrange(desc(NES))
DT_GSEA_MsigDB_c4 %>%  write.xlsx('DT_GSEA_MsigDB_C4_CGN.xlsx')



DT_GSEA_MsigDB_c8<-  GSEA(geneList = genelist, TERM2GENE=C8,  
                           minGSSize    = 20,
                           pvalueCutoff = 0.01, pAdjustMethod = "BH") %>%
  setReadable( OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
DT_GSEA_MsigDB_c8@result %<>%  arrange(desc(NES))
DT_GSEA_MsigDB_c8 %>%  write.xlsx('DT_GSEA_MsigDB_C8.xlsx')

DT_GSEA_MsigDB_c8$ID


#  hallmarker dot plot
DT_GSEA_MsigDB_Hallmark_dot <- DT_GSEA_MsigDB_Hallmark@result %>% 
  mutate(group = case_when(NES > 0 ~ 'BCNMI',  NES < 0 ~'BCMI')) %>% 
  mutate(Description = factor(Description, levels = rev(Description)))
DT_GSEA_MsigDB_Hallmark_dot$NES_abs <-abs( DT_GSEA_MsigDB_Hallmark_dot$NES)


DTDOT <- ggplot(DT_GSEA_MsigDB_Hallmark_dot, aes(x = group, y = Description))+
  geom_point(aes(size = NES_abs, color = -log10(p.adjust)))+
  scale_color_continuous(low = "blue", high = "red")+
  # ggtitle('Enriched Gene Ontology')+
  theme_minimal()+scale_y_discrete(position = "left")+
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8),
        legend.key.width = unit(3,'mm'),
        axis.text.y  = element_text(size=8),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        # legend.position = 'bottom', legend.key.height = unit(2, 'mm'),
        legend.title = element_text(size = 8), title = element_text(size = 8),
        legend.text = element_text(size = 8))+
  guides(size = guide_legend(title.position = 'top',title = 'absolute NES'),
         color = guide_colorbar(title.position = 'top'))


figsave(DTDOT, 'GSEA_enrich_Hallmarks.pdf', 150, 150)


# heatmap
library(pheatmap)
library(ComplexHeatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
      decreasing=TRUE)

select <- BCNMI_vs_BCMI %>% filter(significancy != 'Unsignificant' & abs(log2FoldChange) >3 ) %>% pull(Symbol)

ntd <- normTransform(dds)

color_group <- c('blue', 'red')

names(color_group) <- c('BCNMI',  'BCMI')

anc <- list(color_group = color_group)

df <- as.data.frame(colData(dds))[,c('sample', 'group')]
dt_heat<- pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df, annotation_colors = anc)


dt_heat

ggplotify::as.ggplot(dt_heat)
figsave(ggplotify::as.ggplot(dt_heat), 'heatmap.pdf', 150,270)

Heatmap(as.array(ntd))

as.numeric(assay(ntd))


assay(ntd)['CD3D',]
