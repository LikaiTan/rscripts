################################################################
# Title: RNAseq data analsis pipeline: STAR
# Date: 2019-11-09
# Ver.: 0.1 
################################################################

setwd('/home/big/tanlikai/Demo_projects/DT/CleanData/')

library(edgeR)
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

source('/home/big/tanlikai/script/rscripts/funcs.r')
file.edit('/home/big/tanlikai/script/rscripts/funcs.r')

figsave <- function (p, filename,  w =50, h = 60, device = cairo_pdf,
                     path = 'result/fig/',
                     scale = 1, units = 'mm', dpi = 300
) {
  ggsave2( plot = p, filename = filename,device = device,
           path = path, scale = scale, width = w, height = h,
           units = units, dpi = dpi, limitsize = TRUE)
}



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
after <- apply(allcounts[,2:23],1,function(x) length(x[x>10])>=4)   
qcn <- after[after == 1] %>% names


allcounts_f <- allcounts[qcn,2:17]
# filter MT
allcounts_f <- allcounts_f[grep("^MT-",rownames(allcounts_f),ignore.case=T,invert=T),]
GeneID <- as.data.frame(row.names = rownames(allcounts_f), allcounts_f[,1])
allcounts_f$GeneID <- NULL

########################################
########################################
#groupName<-sampleDetail$Culture.condition  

groupnames <- str_extract(colnames(allcounts_f), '.+(?=_\\d+)')  %>% as.vector() %>% as.factor()
groupnames


# apply(allcounts_f[2:23],1,function(x) x>10  


########################################################################
############# edgeR #############
library("edgeR")

# 构建比较矩阵, unpaired 类型
design <- model.matrix(~0+groupnames)
rownames(design)<-colnames(allcounts_f)
colnames(design)<-levels(factor(groupnames))
##############################normalization##########################################

y <- DGEList(counts=allcounts_f, group=groupnames)
y

keep <- rowSums(cpm(y) > 0.5 ) >=3
y <- y[keep,,keep.lib.sizes=FALSE]
# 这里的0.5(即阈值）等于 10/(最小的文库的 read count数 /1000000)，keep.lib.size=FALSE表示重新计算文库大小
nrow(y)  # 20938  





y <- calcNormFactors(y,method="TMM")
dim(y)
glimpse(y)

# after fitlering
log2counts <- log2(y$counts+1)
hist(log2counts)
boxplot(log2counts, col="gray", las=3)

plotMDS(log2counts)

# PCA and UMAP

scaled_y = scale(t(y$counts)) %>% t() %>% as.data.frame()
library(uwot)

scale(t(y$counts))%>% t() %>% dim()

umap_y <- umap(t(scaled_y)) %>% as.data.frame() %>% `colnames<-`(c('UMAP_1', 'UMAP_2')) %>% `rownames<-`(colnames(y$counts))

umap_y$sample <- rownames(umap_y)

umap_y$group <- str_extract(umap_y$sample, '.+(?=_\\d+)')  %>% as.vector() 

ggplot(umap_y, aes(x = UMAP_1,y = UMAP_2, color = group))+
  geom_text(aes(label=sample))



# 大部分RNAseq使用TMM, 但scRNA不可以
# 如果有很多是0 使用TMMwsp
# 如果一半以上的基因都是差异表达的，也不可以
# plotRLE(log2(cpm(y)+1), outline=FALSE, col=colors[x], las=2)
# plotPCA(log2(cpm(y)+1), col=colors[x],  xlim=c(-0.6,0.6) )                                 # 检查PCA
########################################################################

my.contrasts = makeContrasts(BCMI_BCNMI="BCMI-BCNMI",
                             # MKC_CA_MKC_NA="MKC_CA-MKC_NA",
                             levels=design)

# 估计离散值（Dispersion）。前面已经提到负二项分布（negative binomial，NB)需要均值和离散值两个参数。
# edgeR 对每个基因都估测一个经验贝叶斯稳健离散值（mpirical Bayes moderated dispersion），
# 还有一个公共离散值（common dispersion，所有基因的经验贝叶斯稳健离散值的均值）以及一个趋势离散值

y.Disp <- estimateDisp(y, design, robust = TRUE)

plotBCV(y.Disp)
plotMeanVar(y.Disp, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE, ylim=c(1e+1,1e+12))

# 使用quasi-likelihood (QL) 拟合NB模型
fit <- glmQLFit(y.Disp, design, robust=TRUE)
head(fit$coefficients)

# 当没有design matrix时，estimateDisp等于estimateCommonDisp和estimateTagwiseDisp
# 当有design matrix时，　estimateDisp等价于estimateGLMCommonDisp, estimateGLMTrendedDisp和estimateGLMTagwiseDisp
# y <- estimateGLMCommonDisp(y, design)
# y <- estimateGLMTagwiseDisp(y, design)
# fit <- glmFit(y, design)

id2<-log2(cpm(y)+1)

library(pheatmap)

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=5,scale="none",cluster_cols=T, show_rownames=F, treeheight_col=15, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file="heatmap_all_genes_raw.png",width=12,height=10)
########################################################################


# cpm(y)

# ========================================================================
# MF_noculture              ==============================================
# 如何提高 MF 方法            ==============================================
# ========================================================================
#################################################################

# 如果是glmFit则使用glmLRT, 但是如果是glmQLFit则使用glmQLFTest,这个更严格
lrt.BCMI_BCNMI<- glmQLFTest(fit, contrast=my.contrasts[,"BCMI_BCNMI"])
# lrt.MKC_CA_MKC_NA <- glmQLFTest(fit, contrast=my.contrasts[,"MKC_CA_MKC_NA"])

# 可选矫正FDR
# deg.edger <- lrt.MF_noculture$table[p.adjust(lrt.MF_noculture$table$PValue, method = "BH") < 0.1, ]
# dim(deg.edger)
design
summary(decideTestsDGE(lrt.BCMI_BCNMI, adjust.method="BH", p.value=0.1, lfc=1))  # 看一下大概有多少个上下条, 把BH换成none就是原始的p-value
# summary(decideTestsDGE(lrt.MKC_CA_MKC_NA,adjust.method="BH", p.value=0.1, lfc=1))

deg1 <- decideTestsDGE( lrt.MKC_CA_MKC_NA, adjust.method = 'fdr', p.value =0.1, lfc=1.5)                
deg1.names <-rownames(deg1)[deg1!=0]  

# 或者用常规的提取办法, 但没有lfc threshold
# n<-topTags(lrt.MF_noculture, adjust.method="BH", sort.by="logFC", p.value=0.1, n=Inf)

#################################################################
# final volcano plot           ##################################
#################################################################

library("ggplot2")
library("ggrepel")

# Set it globally, obligatory after Jan. 2021
options(ggrepel.max.overlaps = Inf)

v0<-lrt.MKC_CA_MKC_NA$table
# deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]       # FDR 0.1
deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 ,]                        # show pV better than 0.01 and pV 1.5
nrow(deg.edger0)
# 26 DEGs

degnames0<-rownames(deg.edger0)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames0,]$genelabels<-T                                                        # show only DEGs
v0$change<-as.factor(ifelse(v0$PValue<=0.01 & abs(v0$logFC)>=1.5, ifelse(v0$logFC>=1.5,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(PValue), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 4, segment.color='lightgrey', color='white',fontface='bold', force=0.03, aes(x =logFC, y =-log10(PValue), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle("Volcanoplot of MF vs nocluture") + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1.5, 1.5), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)
volcano0

#################################################################
# ggsave("volvanoplot_MF_vs_nocluture.png",width=10,height=10)
#################################################################

ge<-log2(cpm(y)+1)
u<-c(grep("_MF|_noculture",colnames(ge)))                   #选几列
id2<-ge[degnames0,u]
pheatmap(id2,fontsize=9,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=T, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

########################################################################
#write.table(cbind(deg.edger0,id2),"MF_vs_nonculture_result.csv",sep="\t")



#################################################################
# clusterprofiler @version 2020-01-09       #####################
#################################################################

####### overrepresentation analysis #######
library("clusterProfiler")
library("scales")
library("patchwork")
require(DOSE)

v0<-lrt.MF_noculture$table
deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5,]                         
degnames0<-rownames(deg.edger0)
go0<-enrichGO(degnames0,OrgDb='org.Hs.eg.db',ont='ALL',pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff=0.2,keyType='SYMBOL')
t<-theme(title=element_text(size=8, face='bold',hjust=0)) 
dotplot0<-dotplot(go0, showCategory=20, split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scale="free") + ggtitle("MF vs noculture")+ labs(tag="A")+  t
lfc<-v0[degnames0,]$logFC
names(lfc)<-degnames0
cnet0<-cnetplot(go0, categorySize="pvalue", foldChange=lfc, showCategory = 20)+labs(tag="B") + ggtitle("Network and associated genes")+ t 
dotplot0+cnet0+plot_layout(design = "ABBBBB")

# ggsave("PathwayRegulation_MF_vs_nocluture.png",width=15,height=10)
####### overrepresentation analysis #######
