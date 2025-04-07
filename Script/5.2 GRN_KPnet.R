rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(GENIE3)# GENIE3 Infers a gene regulatory network (in the form of a weighted adjacency matrix) 
library(segmented)
library(igraph)
library(circlize)  # for the colorRamp2 function
library(grid)
library(ComplexHeatmap)

# construct a directed GRN using the same input matrix previously used to generate the co-expression
# network to predict potential regulatory targets of the hub gene Kinases and PPases 
Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%dplyr::select(2:25)%>%t()%>%t() # make a matrix
module_df <- readRDS("../SFCH/Data/WGCNA/netwk13/netwk13_module_df.rds")
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary_df_color.rds") #2011
# pre-analysis for reproducibility of results
set.seed(123)
# run analysis
regulators <- summary_df%>%filter(super_class%in%c("Kinase","PPase")) # 175 kinases; 36 PPases
table(regulators$super_class)

weightMat_all<-GENIE3(Input, 
                      regulators = regulators$ppep.site,
                      targets = NULL,treeMethod = "RF",K = "sqrt",nCores = 12,
                      returnMatrix = TRUE,verbose = FALSE)
#saveRDS("../SFCH/Data/WGCNA/genie3/weightMat_all.rds",object=weightMat_all) #424321=211*2011

weightMat_all <- readRDS("../SFCH/Data/WGCNA/genie3/weightMat_all.rds")
# get regulatory links :424110
# remove fromNode1-toNode from the same protein ID
linkList_all <- getLinkList(weightMat_all)%>%mutate(regProtein=substr(regulatoryGene,1,18),tarProtein=substr(targetGene,1,18))%>%
  filter(regProtein!=tarProtein) # 423768





# generate network for cytoscape #######################################################################################
# calculate edge weight threshold
# plot breakpoint analysis: x is the vector to plot, y is the number of breakpoint to calculate
brkpnt_fun <- function(x, y){
  df <- data.frame(vec = sort(x), num = seq_along(sort(x)))  # generate data frame
  glm_ <- glm(vec ~ num, data = df)  # create linear model
  seg <- segmented(glm_, seg.Z =  ~ num, npsi = y)  # calculate breakpoints
  brkpnt <- seg[["psi"]][, 2] %>% round() # return the 2nd breakpoint
  print(paste0("Breakpoints are ", paste(brkpnt[1:y], collapse = ", ")))  # prints the breakpoints
  thresh <- df[as.numeric(brkpnt),][[1]]  # return the eigencentrality threshold
  ggplot() +  # generate plot
    geom_point(data = df, aes(x = num, y = vec)) +
    geom_vline(xintercept = as.numeric(brkpnt), color='red') +  # xintercept is the breakpoint (psi) +
    theme_classic()
}
brkpnt_fun(linkList_all$weight[seq(from = 1000, to = length(linkList_all$weight),length.out = nrow(linkList_all)/1000)], 3)
#"Breakpoints are 22, 331, 403"
sort(linkList_all$weight[seq(from = 1000, to = length(linkList_all$weight),length.out = nrow(linkList_all)/1000)]) %>% .[403] # 0.01561742
# filter dataframe with eigencentrality threshold
grn_edges_all <- linkList_all[linkList_all$weight > 0.01561742,]%>%
  mutate(regProtein_class=summary_df[as.character(regulatoryGene),"Protein_class"],
         tarProtein_class=summary_df[as.character(targetGene),"Protein_class"]) #21988
#saveRDS("../SFCH/Data/WGCNA/genie3/grn_edges_all.rds",object=grn_edges_all)
grn_edges_all <- readRDS("../SFCH/Data/WGCNA/genie3/grn_edges_all.rds")


# assign protein names to list of network nodes
## define function to generate dataframe
grn_nodes_all <- summary_df%>%filter(ppep.site%in%unique(c(grn_edges_all$regulatoryGene,grn_edges_all$targetGene)))%>%
  dplyr::select(-c(3,8:13,28:33,35:191)) #2011
# calculate GRN hubs ###################################################################################################
#grn_nodes_all <- readRDS("../SFCH/Data/WGCNA/genie3/grn_nodes_all.rds")

grn_nodes_all_eigen <- graph_from_data_frame(grn_edges_all)  # 2011 21988
grn_all_eigen <- eigen_centrality(grn_nodes_all_eigen)  # calculate eigenvector centrality
grn_all_eigen <- data.frame(egnvctr = unname(grn_all_eigen[["vector"]]),gene = names(grn_all_eigen[["vector"]]))
brkpnt_fun(sqrt(grn_all_eigen$egnvctr), 2)  # calculate eigenvector centrality threshold to define GRN hubs
#"Breakpoints are 1360, 1910"
grn_all_hubs <- grn_all_eigen[order(grn_all_eigen$egnvctr),][1910:2011,]
grn_all_hubs <- grn_all_hubs[order(grn_all_hubs$gene),] # sort dataframe by gene
rownames(grn_all_hubs) <- grn_all_hubs$gene 
# 102 hubs
grn_hubs <- data.frame(ppep.site=grn_all_hubs$gene,eigencentrality=grn_all_hubs$egnvctr,summary_df[grn_all_hubs$gene,c(1:2,4:6,14:27,34)])
#saveRDS("../SFCH/Data/WGCNA/genie3/grn_hubs.rds",object=grn_hubs)
table(grn_hubs$color)
# ORA:GRN_hub??? they are already kinase and PPase
# KP mainly in green and turquoise

grn_hubs <- readRDS("../SFCH/Data/WGCNA/genie3/grn_hubs.rds")%>%mutate(GRN_hub="Yes",.after=ppep.site)
colnames(grn_nodes_all)[1] <- "module"
grn_nodes_all<- grn_nodes_all%>%mutate(GRN_hub=grn_hubs[ppep.site,"GRN_hub"],.after=module)
#saveRDS("../SFCH/Data/WGCNA/genie3/grn_nodes_all.rds",object=grn_nodes_all)
table(grn_hubs$color)

grn_nodes_all<- readRDS("../SFCH/Data/WGCNA/genie3/grn_nodes_all.rds")


# output GRN network
write_delim(grn_edges_all, file = "../SFCH/Data/WGCNA/genie3/grn_network_edges.txt", delim = "\t")
write_delim(grn_nodes_all, file = "../SFCH/Data/WGCNA/genie3/grn_network_nodes.txt", delim = "\t")
