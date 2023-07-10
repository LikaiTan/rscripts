library(Seurat)
library(ggplot2)
library(Nebulosa)
library(ggrastr)
library(cowplot)
library(purrr)
library(ggrepel)
library(dplyr)
library(magrittr)
library(tibble)
library(stringr)
library(Cairo)
# themes ------------------------------------------------------------------




mytheme <- theme(plot.title = element_text(size = 8 , face = 'plain'),
                 plot.subtitle = element_text( face = 'plain',size = 8),
                 text = element_text(size = 8, face = 'plain'),
                 legend.title = element_text(size = 8),
                 strip.text = element_text(size = 8),
                 legend.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 axis.line = element_line(size = 0.25),
                 axis.text = element_text(size = 8))
heattheme <-   theme(axis.text.y = element_text(size = 8, face = 'italic'),
                     legend.key.height  = unit(2, 'mm'),
                     legend.position = 'bottom',
                     legend.margin = margin(-7,30,0,0, "mm"),
                     plot.subtitle = element_text(size = 8))

heat_theme <- function(gp,size = 8, legend.position = 'bottom',
                       legend.margin = margin(-7,30,0,0, "mm"),
                       m = 'white', l = 'blue', h = 'red',
                       color = F, fill = NULL, dotsize = F){
  gp+
  theme(axis.text.y = element_text(size = size, face = 'italic'),
        legend.key.height  = unit(2, 'mm'),
        legend.position = legend.position,
        legend.margin = legend.margin,
        plot.subtitle = element_text(size = size,face = 'plain'),
       legend.title = element_text(size = size),
        strip.text = element_text(size = size),
        legend.text = element_text(size = size),
        axis.title = element_text(size = size),
        axis.line = element_line(size = 0.25),
        axis.text = element_text(size = size))+
    scale_fill_gradient2(mid = m, low =l, high = h)+
    guides(color = color, fill = fill, size = dotsize)

}

scale_color_gradient2(low = '#003399', mid = '#ffccff',  high = "#990000")

nox <- theme(axis.text.x = element_blank(),
             axis.ticks.x  = element_blank(),
             axis.title.x = element_blank())

notick <-   theme(axis.ticks = element_blank(), axis.text = element_blank())

hmp <- scale_fill_gradient2(mid = 'white', low = 'blue', high = 'red')



hmp2 <- scale_fill_gradientn( colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(50))

grd <- scale_color_gradientn( na.value = alpha('lightgrey', 0.3),
  colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))

# grd_3c <-     scale_color_gradient2(low = alpha('lightgrey', 0.3), high = 'red', mid = 'purple',
#                                     midpoint = median(FetchData(HumanGDT, x)[,1][FetchData(HumanGDT, x)[,1]>0]) )


#'compare one seurat cluster with another
#gglegendpostion

gglp <- function(p = 'n', s= 3) {
  position <-c('right', 'left', 'top', 'bottom', 'none') %>%
    set_names(c('r', 'l', 't', 'b', 'n'))
  theme(legend.position = position[[p]], legend.key.size = unit(s, 'mm') )
}

# ClusterCompare()



ClusterCompare <- function(ob, id1, id2,log2fc = 0.25,group.by = NULL, rm = "^MT|^RP", test = 'bimod', 
                           p_cutoff = 0.05, assay = 'RNA', do.plot = TRUE, group.colors = NULL, features = NULL,
                           min.pct = 0.1, genetoshow = 50, ds = 500) {
  DefaultAssay(ob) <- assay
    if (!is.null(group.by)) {
      Idents(ob) <- ob[[group.by]]
    }
  result <- c()
  result$table <-  FindMarkers(ob, ident.1 = id1, ident.2 = id2, only.pos = F, features = features,
                               logfc.threshold = log2fc, min.pct = min.pct,
                               test.use = test)%>%
    tibble::rownames_to_column('gene')  %>% dplyr::filter(p_val_adj <= p_cutoff) %>%
    dplyr:: arrange(desc(avg_log2FC))%>%
    dplyr::  mutate(pct.dff = pct.1 - pct.2)

    result$table <- result$table %>%
      dplyr:: filter(!grepl(rm, gene) )

  print(result$table)
  
  if (do.plot == TRUE) {
    result$plot <- DoHeatmap(subset(ob, idents = c(id1, id2), downsample = ds),
                             raster = T,size = gs(8),group.colors = group.colors,
                             features = result$table[c(1:(genetoshow/2),
                                                       (nrow(result$table)-(genetoshow/2-1)):nrow(result$table)
                             ),]$gene)+
      theme(text = element_text(size = 8), axis.text = element_text(size = 8),
            legend.key.width = unit(2,'mm'),
            axis.text.y = element_text(face = 'italic'))+mytheme+
      scale_fill_gradient2(low = 'blue', mid = 'white',  high = "red")+
      guides(color = FALSE)
  }
  

  # print(result$plot)
  return(result)
}



#'feature plot  rasterized
#'Feature_rast()
# Feature_rast is a function to draw scatter plot, the use of this function is similar to the DimPlot and
# FeaturePlot function from Seurat package. However, it's more like a combination of DimPlot and FeaturePlot,
# for it can draw plots with categorical color (like cluter) and  gradient color (like gene expression).
# This function can draw plot from either a Suerat object, or directly from a simple data frame.
# the plot is generated by ggrastr, so it produce rasterized plot. If doing a categorical plot, the labeling of
# cells are done by shadowtext package. Meanwhile , purrr, dplyr, and cowplot packages are also required.

# explain of functions.
# data: the object for visualization. Either a Suerat object, or a dataframe.
# g: the variables you want to visualize, can be a vector. by default it is the ident of a suerat project.
# facet other variables you want to ewrapped into the ggolot object , that you can use facet_grid or facet_wrap to split fig
# sz+ size of dot. dpi: resolution.
# mid.point, when drawing a figure with numeric variable and gradient color, this can assign where should be the middle point. 0.5 means 50%
# ncol: when you want to show a vector of variables,  how many columns you want. Note that we dont have nrow in this function
# mythe: mytheme, a pre-define them in above.
# titleface: the title face of the plot.
# colorset: for  categorical variables, which set of color you wanto to use. It can be 'um", the umap.colors defined above (23 colors), or 'gg', ggplot color, or a vector you define.
# color_grd: for gradient variables, which set of color you like, can choose from three color gradient or 'grd' color
# if choosing 'thresscolor', you can define the three color by following variables : l (low), h (high), m (middle)
# do.label: do you want to labels on your plot or not. label sizeL the size of lable.
# titlesize: the size of title
# othertheme: theme element for ggplot
# d1 and d2: two dimensions, by default is UMAP_1 and UMAP_2
# noaxis: if or not you want to axis
# axis number:  whether or not to lable the axis numbers
# sort: (if gradient) to sort the data frame from low to high
# labels: if have multiple variables to show, you can assign labels for each one.
Feature_rast <- function(data, g = 'ident',facets = NULL, other = NULL,  sz = 0.8,
                         dpi = 300, mid.point = 0.5, ncol = min(5, length(g)), facetcol = NULL,
                         mythe =T, titleface = 'italic',colorset = c('um','gg'), 
                         color_grd = c('A', "B","C", "D",  'threecolor'),
                         do.label = T, labelsize = 10, nrow = NULL, titlesize =8,othertheme = NULL,
                         d1 = "UMAP_1", d2 = 'UMAP_2',noaxis = T, axis.number = F, legendcol = NULL, legendrow=NULL, 
                         labels = NULL, sort =TRUE, assay = DefaultAssay(data),slot = 'data',
                         l = alpha('lightgrey', 0.3), h = 'red', m = 'purple' , navalue =alpha('lightgrey', 0.4) ) {
  if (class(data) == 'Seurat') {
    DefaultAssay(data) <- assay
    fd <- FetchData(data, c(d1, d2,
                            facets , g,other), slot = slot)
  } else {
    fd <- data
  }
  if (sum(as.numeric(grepl('-', g))) != 0) {
    g <- gsub('-', '_', g)
    colnames(fd) <- gsub('-', '_', colnames(fd) )
  } else {
    g<-g
    fd <- fd
  }

  ###one variable
  if (length(g) == 1) {
    (if (isTRUE(is.numeric(fd[[g]]) & isTRUE(sort))){
      fd<-   fd%>%ungroup() %>% arrange_at(.var = vars(contains(g)))

    })
    gp <- ggplot(fd, aes_string(x = d1, y = d2)) +
      geom_point_rast(aes_string(color = g), size = sz, stroke = 0,raster.dpi = dpi) +
      theme_classic()  + 
      ( if (!is.null(facets)) {
        facet_wrap(facets, ncol = facetcol)
      })+
      ( if (isTRUE(mythe)) {
        mytheme
      })+
      #color
      (if (isTRUE(is.numeric(fd[[g]]))){ 
        if (color_grd[1] %in% c("A", "B", "C", "D")  ) {
        scale_color_viridis(discrete = F, option = color_grd[1], na.value = 'transparent')
        
      } else {
        scale_color_gradient2(low = l, high = h, mid = m, na.value = navalue,
                              midpoint = median(fd[[g]][fd[[g]]>0])*mid.point*2)
        }
      } else {
        scale_color_manual(values = (if (colorset[1] == "gg"){
          ggplotColours(length(unique(fd[[g]])))
        } else if (colorset[1] == 'um') {
          umap.colors
        } else {colorset}   ) ,  na.value = navalue) })+
      (if(isTRUE(do.label) & !is.numeric(fd[[g]]))
      { shadowtext::geom_shadowtext(data = (fd %>%  group_by(get(g)) %>% dplyr::select(d1, d2)%>%
                                              dplyr::summarise_all(median) %>%
                                              dplyr::rename_all(.funs = ~c('center', d1,d2))),
                                    mapping = aes_string(label = 'center', x = d1, y =d2), color = 'black',
                                    bg.colour = 'white', size = gs(labelsize))  })+
      ggtitle(g) +
      theme(legend.title  = element_blank(), legend.key.height = unit(4, 'mm'),
            legend.key.width = unit(1,'mm'),
            plot.title = element_text(size = titlesize, face = titleface) )+
      othertheme+
      ( if(is.character(as.vector(fd[[g]]))) {
        guides(color = guide_legend(override.aes = list(size = 2), nrow = legendrow, ncol = legendcol)) })+
      if (isTRUE(noaxis)) {
        NoAxes()
      } else if (!isTRUE(axis.number)) {
        notick
      }


  } else {
    #multiple variables
    gp <- purrr::map(g, function(i) {
      (if (isTRUE(is.numeric(fd[[i]]) & isTRUE(sort) )){
        fd<-   fd%>%ungroup() %>% arrange_at(.var = vars(contains(i)))

      })
      ggplot(fd, aes_string(x = d1, y = d2)) +
        geom_point_rast(aes_string(color = i),size = sz,stroke = 0, raster.dpi = dpi) +
        theme_classic() +
        ( if (!is.null(facets)) {
          facet_wrap(facets, ncol = facetcol)
        })+
        ( if (isTRUE(mythe)) {
          mytheme
        })+
        #color
        (if (isTRUE(is.numeric(as.vector(fd[[i]])))){
          if (color_grd[1] %in% c("A", "B", "C", "D")  ) {
          scale_color_viridis(discrete = F, option = color_grd[1], na.value = 'transparent')
          
        } else {
          scale_color_gradient2(low = l, high = h, mid = m,na.value = navalue,
                                midpoint = median(fd[[i]][fd[[i]]>0])*mid.point*2) }
        } else {
          scale_color_manual(values = (if (colorset[1] == "gg"){
            ggplotColours(length(unique(fd[[i]])))
          } else if (colorset[1] == 'um') {
            umap.colors
          } else {colorset}   ), na.value = navalue)
        })+
        (if(isTRUE(do.label) & !is.numeric(fd[[i]]))
        { shadowtext::geom_shadowtext(data = (fd %>%  group_by(get(i)) %>% dplyr::select(d1, d2)%>%
                                                dplyr::summarise_all(median) %>%
                                                dplyr::rename_all(.funs = ~c('center', d1,d2))),
                                      mapping = aes_string(label = 'center', x = d1, y =d2), color = 'black',
                                      bg.colour = 'white', size = gs(labelsize))  })+
        ggtitle(i) +
        theme(legend.title  = element_blank(), legend.key.height = unit(4, 'mm'),
              legend.key.width = unit(1,'mm'),
              plot.title = element_text(size = titlesize, face = titleface) )+
        othertheme+
        ( if(is.character(as.vector(as.vector(fd[[i]])))) {
          guides(color = guide_legend(override.aes = list(size = 2))) })+
        if (isTRUE(noaxis)) {
          NoAxes()
        }  else if (!isTRUE(axis.number)) {
          notick
        }

    }) %>% PG(nrow = nrow, ncol = ncol,  align = 'v', labels = labels )

  }
  return(gp)
}





Feature_density <- function(data, feature = NULL,sz = 0.5,  pal = "viridis", reduction = 'umap',
                             ncol = min(5, length(feature)), joint =F, method = c("ks", "wkde"),
                            adjust = 1,shape = 16, nrow = NULL,  othertheme = NULL,
                            mythe =T, titleface = 'italic',titlesize =8,
                           noaxis = T, axis.number = F,
                            labels = NULL,  assay = DefaultAssay(data),slot = NULL ) {
  
  DefaultAssay(data) <- assay
  if (length(feature) ==1 ) {
   
 gp<- plot_density(object=data, features=feature, joint = joint, reduction = reduction,
                             pal =pal, slot = slot, size =sz, method = method)+
      ( if (isTRUE(mythe)) {
        mytheme
      })+     
   theme(legend.key.height = unit(4, 'mm'),
                     legend.key.width = unit(1,'mm'),
                     plot.title = 
           element_text(size = titlesize, face = titleface) )+
   othertheme+ 
   if (isTRUE(noaxis)) {
     NoAxes()
   } else if (!isTRUE(axis.number)) {
     notick
   }
  } else {
   gp<- plot_density(object=data, features=feature, joint = joint, combine = F,reduction = reduction,
                     pal =pal, slot = slot, size =sz, method = method)
   gp <- map(gp, ~ .x +
               ( if (isTRUE(mythe)) {
                 mytheme
               })+      
               theme(legend.key.height = unit(4, 'mm'),
                     legend.key.width = unit(1,'mm'),
                     plot.title = element_text(size = titlesize,
                                               face = titleface) )+
               othertheme+ 
               if (isTRUE(noaxis)) {
                 NoAxes()
               } else if (!isTRUE(axis.number)) {
                 notick
               })
   
   gp <- gp %>%  PG(nrow = nrow, ncol = ncol,  align = 'v', labels = labels )
 }
  
  return(gp)
} 









#'ggplot color scheme
#'ggplotColours()
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#'a function to unify the size of geom_text and element_text
#'
#'geom text/label size
#'gs()
gs  <- function(x) {
  x*5/14
}

gs(8)

#'improved plot_grid
#'PG()
PG <- function (x, align = c("none", "h", "v", "hv"),
                axis = c("none", "l", "r", "t", "b", "lr", "tb", "tblr"),
                nrow = NULL, ncol = NULL, rw = 1, rh = 1,
                labels = NULL, label_size = 9, label_fontfamily = NULL,
                label_fontface = "bold", label_colour = NULL, label_x = 0,
                label_y = 1, hjust = 0, vjust = 1.2, scale = 0.93,
                greedy = TRUE, cols = NULL, rows = NULL) {
  plot_grid(plotlist = x, align = align, axis = axis, nrow = nrow, ncol =ncol, rel_widths =rw,
            rel_heights = rh, labels = labels, label_size = label_size, label_fontfamily = label_fontfamily,
            label_fontface = label_fontface, label_colour = label_colour, label_x = label_x, label_y= label_y, hjust = hjust,
            cols = cols, rows = rows)
}


#'sample with fixed seed
#'set_samle()
set_sample <- function(x,  n = NULL, s = 629)  {
  set.seed(s)
  if (is.null(n)) {
    set.seed(s)
    sample(x)
  } else {
    set.seed(s)
    sample(x, n)
  }
}



#'function to save objects
#'save object as rdata
#'saverdata()

# # saverdata <- function(x, file = NULL) {
#   otc <- function(x) {
#     c <-  as.character(substitute(x))
#     if(length(c) > 1) {
#       return(c[-1])
#     } else {
#       return(c)
#     }
#   }
#
#   save(list = otc(x), file = paste0(file,Sys.Date(),'.rdata'))
# }


#'save list
#'otc()
# otc <- function(x) {
#   c <-  as.character(substitute(x))
#   if(length(c) > 1) {
#     return(c[-1])
#   } else {
#     return(c)
#   }
# }




#'function to produce multiple GSEA term plot
#'GSEA_multipplot()

GSEA_multipplot <- function(x, description_to_show, legendpvalue = F,
                            rel_h = c(1, 0.1, 0.4),
                            legend.position = 'bottom', c1,c2,
                            base_size = 8, title = paste('GSEA',c1,'vs',c2), col = "green") {
  require(enrichplot)
  GSEAall_1 <- gseaplot2(x, color = col,
                         geneSetID = description_to_show,
                         base_size = base_size, pvalue_table = legendpvalue, subplots = 1) +
    theme(legend.position = 'none', line = element_line(size = 0.5),axis.line.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.line.y = element_line(size = 0.25),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())+
    geom_hline(aes(yintercept=0), linetype="dashed", size = 0.3)+
    ggtitle(title)+
    theme(plot.title = element_text(size = base_size, hjust = 1, vjust = 0.5),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()
    )
  GSEAall_2 <- gseaplot2(x, color = col,
                         geneSetID = description_to_show,base_size = base_size, pvalue_table = F, subplots = 2)+
    theme(legend.position = 'none', line = element_line(size = 0.5),axis.line.x = element_blank(),
          axis.line.y = element_line(size = 0.25),
          panel.grid = element_blank(),
          plot.margin=margin(l=-0.8,unit="cm"),
          axis.ticks = element_blank(), axis.text.x = element_blank())
  GSEAall_3 <- gseaplot2(x, color = col,
                         geneSetID = description_to_show,base_size = base_size,
                         pvalue_table = F, subplots = 3)+
    theme(legend.position = 'none',           axis.line = element_line(size = 0.25),
          axis.text = element_text(size = 8),
          plot.margin=margin(l=-0.8,unit="cm"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8))+
    scale_x_continuous(breaks = c(1, length(x@geneList)), expand = c(0, 0),
                       labels=c(paste(c1,'high'),paste(c2,'high')))

  GESA_legend <-
    cowplot::get_legend(gseaplot2(x, color = col,
                                  geneSetID = description_to_show,base_size = base_size,
                                  pvalue_table = legendpvalue, subplots = 1)+
                          theme(legend.position = 'right', legend.text  = element_text(size = 8)))
  GSEA_all <- plot_grid(plotlist = list(GSEAall_1, GSEAall_2, GSEAall_3), ncol = 1, scale = 1,
                        align = "v", axis = 'y', rel_heights = rel_h)
  if (legend.position == 'bottom'){
    GSEA_all <- plot_grid(plotlist = list(GSEA_all, GESA_legend), ncol =1,
                          rel_heights = c(1, 0.12)) %>% ggplotify::as.ggplot()
  } else   {
    if (legend.position != 'no')
      GSEA_all <- plot_grid(plotlist = list(GSEA_all, GESA_legend), ncol =2) %>%
        ggplotify::as.ggplot()
  }
  return(GSEA_all)
}

#' function to produce entrezlist compareing two seurat cluster
#' entrezlist_generator()

entrezlist_generator <- function(x, id1, id2, OrgDB = c('org.Hs.eg.db'),rm = "^MT|^RP") {
  logfclist <- FindMarkers(object = x, ident.1 = id1, ident.2 = id2,
                           test.use = 'bimod',logfc.threshold = 0, 
                           only.pos = F,  min.pct = 0.1) %>%
    tibble::rownames_to_column('SYMBOL') %>% dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr:: filter(!grepl(rm, SYMBOL) )
  logfclist <- clusterProfiler::bitr(logfclist$SYMBOL, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb=OrgDB) %>%
    right_join(logfclist, by = 'SYMBOL')
  entrezidlist <- logfclist$avg_log2FC
  names(entrezidlist) <- as.character(logfclist$ENTREZID)
  return(sort(entrezidlist, decreasing = T))
}


figsave <- function (p, filename,  w =50, h = 60, device = cairo_pdf,
                     path = 'figs',
                    scale = 1, units = 'mm', dpi = 300
) {
  ggsave2( plot = p, filename = filename,device = device,
           path = path, scale = scale, width = w, height = h,
           units = units, dpi = dpi, limitsize = TRUE)
}



ViolinPlot <- function(data, g, sz = 0.5, dpi = 300,
                        group.by = NULL,
                        facet = NULL,
                       ncol = min(3, length(g)), split = NULL, 
                       colors = ggplotColours(cln), othertheme = NULL,
                       idents = NULL,alpha_point =0.8, alpha_fill = 0.4, jitter = T, box = F,
                       x.angle = 0, width = 0.25, Plotgrid = T, ylabtext ='\nexpression',size = 8,
                       assay = DefaultAssay(data),slot = 'data',
                       labels = NULL, labelsize =8, labelface='plain',
                       mythe =F, titleface = 'italic'){
  if (length(g) == 1) {
    fig <- VlnPlot(data, g, pt.size = 0, idents = idents, group.by = group.by,
                   split.by = facet, assay = assay, slot = slot)
    cln = fig$data$ident %>% unique() %>% length()
                  fig +
            (if(isTRUE(jitter)){
              geom_jitter_rast (size = sz,  raster.dpi = 300, stroke = 0,
                                width = width, aes(color = ident) )
            })     +
                    (if (isTRUE(box)){
   geom_boxplot( alpha = 0.5, size = 0.3,  width = 0.5, outlier.alpha = 0)
                    })+
      xlab(NULL)+ylab(paste(g,ylabtext))+theme_minimal()   +
                    ( if (isTRUE(mythe)) {
                      mytheme
                    })+
     scale_fill_manual(values = alpha(colors, alpha_fill))+
      scale_color_manual(values = alpha(colors, alpha_point))+
      theme(legend.position = 'none',
            axis.text.x = element_text(angle = x.angle, size = size),
            axis.text.y = element_text(size = size),
            plot.title = element_blank(),
            axis.line = element_blank(),
            axis.title.y = element_text(face = titleface,size = size))+
                    othertheme
  } else {
    gp <- lapply(g, function(i) {
      fig <- VlnPlot(data, i, pt.size = 0, idents = idents, group.by = group.by,split.by = facet,assay = assay, slot = slot)
    cln = fig$data$ident %>% unique() %>% length()
    fig +
      (if(isTRUE(jitter)){
        geom_jitter_rast (size = sz,  raster.dpi = 300,stroke = 0,
                          width = width, aes(color = ident) )
      })     +
      (if (isTRUE(box)){
        geom_boxplot( alpha = 0.5, size = 0.3,  width = 0.5, outlier.alpha = 0)
      })+   xlab(NULL)+ylab(paste(i,ylabtext))+theme_minimal()   +
      ( if (isTRUE(mythe)) {
        mytheme
      })+
      scale_fill_manual(values = alpha(colors, alpha_fill))+
      scale_color_manual(values = alpha(colors, alpha_point))+
      theme(legend.position = 'none',
            axis.text.x = element_text(angle = x.angle, size = size),
            axis.text.y = element_text(size = size),
            plot.title = element_blank(),
            axis.line = element_blank(),
            axis.title.y = element_text(face = titleface,size = size))+
      othertheme

    })
    if (isTRUE(Plotgrid) ) {
      gp <-PG(gp,ncol = ncol, align = 'v', labels = labels, label_fontface = labelface, label_size = labelsize)
    } else {
      gp <-setNames(gp, g)
    }

    return(gp)
  }
}

nr <- function(x, a,b) {
  x >=a & x <= b

}
nr(10.1,5,11)
nr(c(1,2,3,4),2,4)



#label
do.label <- function(data = NULL, label = "center", color = 'black',
                     bg.colour = 'white', size = 10  ){
  shadowtext::geom_shadowtext(data = data,  mapping = aes_string(label = label), color = color,
                                bg.colour = bg.colour, size = gs(size))

}



# umapcolors --------------------------------------------------------------

umap.colors <- c(
  "#0F95B9",
  "#B4DC49", 
  "#EC5ECE", 
  "#FBD64A" ,
  "#638B83", 
  "#6FD6E8",
  "#CF5046",
  "#1F405C" ,
  "#F2895E",
  "#A35A33",
  "#DE342F", 
  "#DB8A0F", 
  "#A33A43",
  "#7D6C86", 
  "#D0D0D0" ,
  "#E7C595", 
  "#A34F23", 
  "#0F95DA", 
  "#5E2870" ,
  "#59BF30",
  "#A6E9DB", 
  "#9D43BB", 
  "#FB6C46"
)



color_m <- function(color = umap.colors, al =1,
                    na = alpha('lightgrey',0.5), labels = waiver()) {
  scale_color_manual(values = alpha(color, al), na.value = na, labels = labels)
}
fill_m <- function(color = umap.colors, al =1, na = alpha('lightgrey',0.5),labels = waiver()) {
  scale_fill_manual(values = alpha(color, al), na.value = na,labels = labels)
}
?scale_fill_manual


# options(future.fork.enable = TRUE)
# options(future.globals.maxSize= 6012896000)
# future::plan(strategy = "multicore", workers = 40)


multicores <- function(core=20, mem = 100, strategy = 'multicore') {
  options(future.fork.enable = TRUE)
  options(future.globals.maxSize= mem*1024^3,future.seed=TRUE)
  future::plan(strategy = strategy, workers = core)
}

# multicores(mem = 200)



# Hifreq_define <- map(id,  function(x) for (
#   i in 1:nrow(filter(TRDfreq_id, orig.ident == x))){
#   if (sum(filter(TRDfreq_id, orig.ident == x )$n[1:i])/sum(filter(TRDfreq_id, orig.ident == x )$n) >= 0.75 ) {
#     return(i/nrow(filter(TRDfreq_id, orig.ident == x))*100)
#     break
#   }
#   
# }
# ) %>% set_names(id) %>% as.data.frame()






























