# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)


GDTlung_s@meta.data %<>%  mutate(pheno = str_replace(pheno, "memory", "circ") ) 

TRD_data_meta <- GDTlung_s@meta.data %>%  filter(cdr3_paired_freq > 1 & !(Cell_cluster %in% c("Vg9Vd2_Mix", "gd_Naive_LN3"))) %>%  filter((grepl("_LG", Cell_cluster))  & tissue == "Lung" | !(grepl("_LG", Cell_cluster) ) & tissue == "LLN")  %>% mutate(pheno = factor(pheno, levels= c("Lung_TRM", "Lung_circ", "LN_circ"))) 

TRD_data <- TRD_data_meta %>%  group_by(pheno) %>%  dplyr::count(cdr3_paired) %>% arrange(pheno, desc(n))%>% mutate(pheno = factor(pheno, levels= c("Lung_TRM", "Lung_circ", "LN_circ"))) %>%   filter(n>1 )

GDTlung_s$cdr3_paired_freq

TRD_data$uniquename <-paste(TRD_data$pheno, TRD_data$cdr3_paired)




# paired TCR between phenos
connect <- TRD_data %>%
  inner_join(TRD_data, by = "cdr3_paired", suffix = c(".from", ".to")) %>%
  # Ensure that we don't pair a row with itself
  filter(uniquename.from < uniquename.to) %>%
  select(from = uniquename.from, to = uniquename.to)
connect

# create a vertices data.frame. One line per object of our hierarchy
# create the edge (clones in each pheno) 

d1 <- data.frame(from="origin", to= c("Lung_TRM", "Lung_circ", "LN_circ"))
d2 <- data.frame(from=TRD_data$pheno, to=TRD_data$uniquename)
edges <- rbind(d1, d2)
nrow(vertices)

# vertices  <-  data.frame(
#   # name is each clone
#   name = unique(c(as.character(edges$from), as.character(edges$to))) , 
# # the value in vertices is the clonal size!!
#     value = c(NA,NA,NA,NA,TRD_data$n)
# ) 
# # Let's add a column with the group of each name. It will be useful later to color points 
# vertices$group  <-  edges$from[ match( vertices$name, edges$to ) ]
# vertices

vertices <- data.frame(name = c("origin", "Lung_TRM", "Lung_circ", "LN_circ", TRD_data$uniquename),
                       # the value in vertices is the clonal size!!
                       value = c(NA,NA,NA,NA,TRD_data$n),
                       group = c(NA,"origin","origin","origin", as.character(TRD_data$pheno )))

vertices




mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
glimpse(mygraph)
# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)

F2D <-  (ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 

  geom_node_point(aes(filter = leaf, x = x, y=y, colour=group, size=value, alpha=0.2)) +
  scale_edge_colour_manual(values= umap.colors, guide = F) +
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.5, width=0.5, aes(colour=group)) +
  scale_colour_manual(values= umap.colors) +
  scale_size_continuous( range = c(0.1,5) ) +
  mytheme+ggtitle("Paired TCR sharing")+
 
  theme(
 
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +guides(alpha = "none", group =  guide_legend(order = 1))+
    labs(value = "clonal size")+
    mytheme+ theme_void(base_size = 6) +
    NULL
) %T>% print()




test <-  get_legend(F2D)

test
test




 