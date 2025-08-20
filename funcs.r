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
library(RColorBrewer)
# themes ------------------------------------------------------------------



# mytheme -----------------------------------------------------------------



mytheme <- theme(plot.title = element_text(size = 6 , face = 'plain'),
                 plot.subtitle = element_text( face = 'plain',size = 6),
                 text = element_text(size = 6, face = 'plain'),
                 legend.title = element_text(size = 6),
                 strip.text = element_text(size = 6),
                 legend.text = element_text(size = 6),
                 axis.title = element_text(size = 6),
                 axis.line = element_line(size = 0.25),
                 axis.text = element_text(size = 6))
heattheme <-   theme(axis.text.y = element_text(size = 6, face = 'italic'),
                     legend.key.height  = unit(2, 'mm'),
                     legend.position = 'bottom',
                     legend.margin = margin(-7,30,0,0, "mm"),
                     plot.subtitle = element_text(size = 6))

heat_theme <- function(gp,size = 6, legend.position = 'bottom',
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



ClusterCompare <- function(ob, id1, id2,log2fc = 0.25,group.by = NULL,
                           rm = "^MT|^RP", test = 'bimod', angle = 20,
                           p_cutoff = 0.05, assay = 'RNA', slot = "data",
                           do.plot = TRUE, group.colors = NULL, features = NULL,
                           min.pct = 0.1, genetoshow = 50, ds = 500) {
  DefaultAssay(ob) <- assay
    if (!is.null(group.by)) {
      ob <-  SetIdent(ob, value = group.by)
    }
  result <- c()
  result$table <-  FindMarkers(ob, ident.1 = id1, ident.2 = id2, only.pos = F, features = features,
                               logfc.threshold = log2fc, min.pct = min.pct, slot = slot,
                               test.use = test)%>%
    tibble::rownames_to_column('gene')  %>% dplyr::filter(p_val_adj <= p_cutoff) %>% dplyr:: arrange(desc(avg_log2FC )) %>% dplyr::  mutate(pct.dff = pct.1 - pct.2)

    result$table <- result$table %>%
      dplyr:: filter(!grepl(rm, gene) )

  print(result$table)
  
  if (do.plot == TRUE) {
    result$plot <- DoHeatmap(subset(ob, idents = c(id1, id2), downsample = ds),angle = angle,
                             raster = T,size = gs(8),group.colors = group.colors,
                             features = result$table[c(1:(genetoshow/2),
                                                       (nrow(result$table)-(genetoshow/2-1)):nrow(result$table)
                             ),]$gene)+
      theme(text = element_text(size = 6), axis.text = element_text(size = 6),
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
                         dpi = 300, mid.point = 0.5, ncol = min(5, length(g)), 
                         facetcol = NULL,
                         mythe =T,
                         titleface = 'italic',colorset = c('um','gg'),    colorgrd = "grd1",
               
                         do.label = T, labelsize = 10, nrow = NULL, titlesize =6,othertheme = NULL,
                         d1 = "UMAP_1", d2 = 'UMAP_2',noaxis = T, axis.number = F, legendcol = NULL, legendrow=NULL, 
                         labels = NULL, sort =TRUE, assay = DefaultAssay(data),slot = 'data',
                      
                         navalue ="transparent" ) {

  if (class(data)[1] == 'Seurat') {
    DefaultAssay(data) <- assay
    fd <- FetchData(data, c(d1, d2,
                            facets , g,other), layer = slot)
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

  # gradient color define 
  color_list1 <- c( alpha(c("#D4EDF7", "#347B99"), 0.5), "#4424D6", "#110934")  # Example colors
  color_list2 <- c(alpha(c("#F0F7D4", "#B2D732"),0.5), "#347B11", "#092834")  # Example colors
  
  # Determine the color gradient to use
  if (length(colorgrd) == 1 && colorgrd == 'grd1') {
    colors <- color_list1
  } else if ( length(colorgrd) == 1 && colorgrd == 'grd2') {
    colors <- color_list2
  } else if (length(colorgrd) > 1 && is.vector(colorgrd) ) {
    colors <- colorgrd
  } else {
    stop("Invalid colorgrd value. Use 1, 2, or a vector of color hex codes.")
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
        scale_colour_gradientn(
          colors = colors,
          space = "Lab", 
          na.value = navalue,
          guide = "colourbar",
          aesthetics = "colour" )
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
          scale_colour_gradientn(
            colors = colors,
            space = "Lab",
            na.value = navalue,
            guide = "colourbar",
            aesthetics = "colour" )
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
                            mythe =T, titleface = 'italic',titlesize =6,
                           noaxis = T, axis.number = F,
                           colorgrd =  "grd1",navalue= "transparent",
                           
                            labels = NULL,  assay = DefaultAssay(data),slot = NULL ) {
  
  # gradient color define 
  color_list1 <- c( alpha(c("#D4EDF7", "#347B99"), 0.5), "#4424D6", "#110934")  # Example colors
  color_list2 <- c(alpha(c("#F0F7D4", "#B2D732"),0.5), "#347B11", "#092834")  # Example colors
  # Determine the color gradient to use
  if (length(colorgrd) == 1 && colorgrd == 'grd1') {
    colors <- color_list1
  } else if ( length(colorgrd) == 1 && colorgrd == 'grd2') {
    colors <- color_list2
  } else if (length(colorgrd) > 1 && is.vector(colorgrd) ) {
    colors <- colorgrd
  } else {
    stop("Invalid colorgrd value. Use 1, 2, or a vector of color hex codes.")
  }
  
  
  DefaultAssay(data) <- assay
  if (length(feature) ==1 ) {
   
 gp<- plot_density(object=data, features=feature, joint = joint, reduction = reduction,
                             pal =pal, slot = slot, size =sz, method = method)+
      ( if (isTRUE(mythe)) {
        mytheme
      })+     
   scale_colour_gradientn(
     colors = colors,
     space = "Lab",
     na.value = navalue,
     guide = "colourbar",
     aesthetics = "colour" )+
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
               scale_colour_gradientn(
                 colors = colors,
                 space = "Lab",
                 na.value = navalue,
                 guide = "colourbar",
                 aesthetics = "colour" )+
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









#'ggplot color scheme generate ggplot colors, random or not
#'ggplotColours()
ggplotColours <- function(n = 6, h = c(0, 360) + 15, r = F, seed = 1){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
 cl <-  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
 if (isTRUE(r)) {
cl <- set_sample(cl, s = seed) }
 return(cl)
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
                labels = NULL, label_size = 8, label_fontfamily = NULL,
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
                            base_size = 6, title = paste('GSEA',c1,'vs',c2), 
                            plots = c(1,2,3),
                            col = "green") {
  require(enrichplot)
  GSEAall_1 <- gseaplot2(x, color = col,
                         geneSetID = description_to_show,
                         base_size = base_size, pvalue_table = legendpvalue, subplots = 1) +
    theme(legend.position = 'none', line = element_line(size = 0.5),axis.line.x = element_blank(),
          axis.text.y = element_text(size = 6),
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
          axis.text = element_text(size = 6),
          plot.margin=margin(l=-0.8,unit="cm"),
          panel.grid = element_blank(),
          axis.title = element_text(size = 6))+
    scale_x_continuous(breaks = c(1, length(x@geneList)), expand = c(0, 0),
                       labels=c(paste(c1,'high'),paste(c2,'high')))
  
    GSEA_all_list <-  list(GSEAall_1, GSEAall_2, GSEAall_3)
  GESA_legend <-
    cowplot::get_legend(gseaplot2(x, color = col,
                                  geneSetID = description_to_show,base_size = base_size,
                                  pvalue_table = legendpvalue, subplots = 1)+
                          theme(legend.position = 'right', legend.text  = element_text(size = 6)))
  GSEA_all <- plot_grid(plotlist = GSEA_all_list[plots], ncol = 1, scale = 1,
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

entrezlist_generator <- function(x, id1, id2, OrgDB = c('org.Hs.eg.db'),rm = "^MT|^RP", sort = "avg_log2FC") {
  logfclist <- FindMarkers(object = x, ident.1 = id1, ident.2 = id2,
                           test.use = 'bimod',logfc.threshold = 0, 
                           only.pos = F,  min.pct = 0.1) %>% mutate(sig = avg_log2FC*-log10(p_val_adj)) %>% 
    tibble::rownames_to_column('SYMBOL') %>% 
    # dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr:: filter(!grepl(rm, SYMBOL) )
  logfclist <- clusterProfiler::bitr(logfclist$SYMBOL, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb=OrgDB) %>%
    right_join(logfclist, by = 'SYMBOL')
  entrezidlist <- logfclist[[sort]]
  names(entrezidlist) <- as.character(logfclist$ENTREZID)
  return(sort(entrezidlist, decreasing = T))
}



Genelist_generator <- function(x, id1, id2, rm = "^MT|^RP",sort = "avg_log2FC", group.by = NULL) {
  
  if (!is.null(group.by)) {
   x <-  SetIdent(x, value = group.by)
  }
 
   logfclist <- FindMarkers(object = x, ident.1 = id1, ident.2 = id2, 
                           test.use = 'bimod',logfc.threshold = 0, 
                           only.pos = F,  min.pct = 0.1) %>%
    tibble::rownames_to_column('SYMBOL') %>% 
    mutate( # Replace p_val_adj = 0 with the smallest positive double to avoid -log10(0) = Inf
      p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
      sig = avg_log2FC * -log10(p_val_adj)) %>% 
    # dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr:: filter(!grepl(rm, SYMBOL) )
  genelist <- logfclist[[sort]]
  names(genelist) <- as.character(logfclist$SYMBOL)
  return(sort(genelist, decreasing = T))
}




figsave <- function (p, filename,  w =50, h = 60, 
                     path = 'figs',
                    scale = 1, units = 'mm', dpi = 300
) {
  ggsave2( plot = p, filename = filename,
           path = path, scale = scale, width = w, height = h,
           units = units, dpi = dpi, limitsize = TRUE)
}



ViolinPlot <- function(data, g, sz = 0.5, dpi = 300,
                        group.by = NULL,
                        facet = NULL,
                       ncol = min(3, length(g)), split = NULL, 
                       colors = ggplotColours(cln), othertheme = NULL,
                       idents = NULL,alpha_point =0.8, alpha_fill = 0.4, 
                       jitter = T, box = F,
                       x.angle = 0, width = 0.25, Plotgrid = T, ylabtext ='\nexpression',size = 6,
                       assay = DefaultAssay(data),slot = 'data',
                       labels = NULL, labelsize =6, labelface='plain',
                       mythe =T, titleface = 'italic'){
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


scale_m <- function( colorgrd = "grd1", 
                     
                     breaks=NULL,
                     limits=NULL,
                     na = "transparent"
  
) {
  
  # gradient color define 
  color_list1 <- c( alpha(c("#D4EDF7", "#347B99"), 0.5), "#4424D6", "#110934")  # Example colors
  color_list2 <- c(alpha(c("#F0F7D4", "#B2D732"),0.5), "#347B11", "#092834")  # Example colors
  
  # Determine the color gradient to use
  if (length(colorgrd) == 1 && colorgrd == 'grd1') {
    colors <- color_list1
  } else if ( length(colorgrd) == 1 && colorgrd == 'grd2') {
    colors <- color_list2
  } else if (length(colorgrd) > 1 && is.vector(colorgrd) ) {
    colors <- colorgrd
  } else {
    stop("Invalid colorgrd value. Use 1, 2, or a vector of color hex codes.")
  }  
  scale_colour_gradientn(
    colors = colors,
    space = "Lab", 
    na.value = na,
    guide = "colourbar",
    breaks = breaks,
    limits = limits,
    aesthetics = "colour" )
  
}




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






# DoMultiHeatmap, revised from https://github.com/satijalab/seurat/issues/2201 and https://github.com/elliefewings/DoMultiBarHeatmap/blob/main/R/domultiheatmap.func.R -------------------------------------------
suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is.null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is.null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is.null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is.null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.unique <- paste(group.use[[colname]], group.use[[group.by]], sep="+=$")
          print(length(group.use$x))
          print(length(label.unique))
          
          label.x.pos <- tapply(X = group.use$x, INDEX = label.unique,
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          label.x.pos$group <- label.x.pos$group %>% lapply( function(x) gsub("\\+\\=\\$.*","", x))
          
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}


##' Plot TCR Sharing Across Phenotypes as a Circular Dendrogram
#' 
#' @param tcr_data A data frame with columns: pheno, TCR, n, uniquename
#' @param colors Color palette for visualization (default = rev(brewer.pal(6, "Set1")))
#' @param mytheme Custom ggplot theme (default = theme_minimal())
#' @param title Plot title (default = "Paired TCR Sharing")
#' @return A ggraph plot object
#' @export
plot_tcr_sharing <- function(tcr_data, 
                             colors = rev(RColorBrewer::brewer.pal(6, "Set1")),
                             mytheme = ggplot2::theme_minimal(),
                             size_range = c(0.1, 5),
                             title = "TCR Sharing") {
  
  # Load required libraries (comment out if already loaded)
  # library(ggraph)
  # library(igraph)
  # library(tidyverse)
  # library(RColorBrewer)
  colnames(tcr_data) <-  c("pheno", "TCR", "n", "uniquename")
  # Ensure input data has required columns
  required_cols <- c("pheno", "TCR", "n", "uniquename")
  if (!all(required_cols %in% colnames(tcr_data))) {
    stop("Input data must contain columns: pheno, TCR, n, uniquename")
  }
  
  # Ensure uniquename is unique and convert to character for consistency
  tcr_data <- tcr_data %>%
    dplyr::mutate(uniquename = as.character(uniquename),
                  pheno = as.character(pheno)) %>%
    dplyr::distinct(uniquename, .keep_all = TRUE)
  
  # Create connections for shared TCRs
  connect <- tcr_data %>%
    dplyr::inner_join(tcr_data, by = "TCR", suffix = c(".from", ".to"), 
                      relationship = "many-to-many") %>%
    dplyr::filter(uniquename.from < uniquename.to) %>%
    dplyr::select(from = uniquename.from, to = uniquename.to)
  
  # Define edges
  unique_phenos <- unique(tcr_data$pheno)
  d1 <- data.frame(from = "origin", to = unique_phenos, stringsAsFactors = FALSE)
  d2 <- data.frame(from = tcr_data$pheno, to = tcr_data$uniquename, stringsAsFactors = FALSE)
  edges <- rbind(d1, d2)
  
  # Define vertices explicitly including all required names
  vertices <- data.frame(
    name = c("origin", unique_phenos, tcr_data$uniquename),
    value = c(NA, rep(NA, length(unique_phenos)), tcr_data$n),
    group = c(NA, rep("origin", length(unique_phenos)), as.character(tcr_data$pheno)),
    stringsAsFactors = FALSE
  )
  
  # Remove any duplicate vertices (just in case)
  vertices <- vertices %>% dplyr::distinct(name, .keep_all = TRUE)
  
  # Debug: Check for mismatches
  missing_in_vertices <- setdiff(tcr_data$uniquename, vertices$name)
  if (length(missing_in_vertices) > 0) {
    warning("These uniquenames are missing from vertices: ", 
            paste(missing_in_vertices, collapse = ", "))
  }
  
  # Create igraph object
  mygraph <- igraph::graph_from_data_frame(edges, vertices = vertices)
  
  # Generate plot
  plot <- ggraph(mygraph, layout = "dendrogram", circular = TRUE) +
    geom_node_point(aes(filter = leaf, x = x, y = y, colour = group, size = value, alpha = 0.2)) +
    scale_edge_colour_manual(values = colors, guide = FALSE) +
    geom_conn_bundle(data = get_con(from = match(connect$from, vertices$name),
                                    to = match(connect$to, vertices$name)), 
                     alpha = 0.5, width = 0.5, aes(colour = group)) +
    scale_colour_manual(values = colors) +
    scale_size_continuous(range = size_range) +
    mytheme +
    ggtitle(title) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(alpha = "none", group = guide_legend(order = 1)) +
    labs(size = "clonal size") +
    theme_void(base_size = 6)
  
  # Print and return plot
  print(plot)
  return(plot)
}
# Example usage with your pre-processed data:

