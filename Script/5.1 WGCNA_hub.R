rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(WGCNA) #“weighted gene correlation networks for analysis”
library(segmented)
library(doParallel)
library(cowplot)
library(ComplexHeatmap)
library(igraph) # remove loop and duplicated undirected edges

net_edges <- readRDS("../SFCH/Data/WGCNA/netwk13/net_edges.rds")# 150335
net_nodes <- readRDS("../SFCH/Data/WGCNA/netwk13/net_nodes.rds")%>%dplyr::select(-3)# 1994
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
module_df <- readRDS("../SFCH/Data/WGCNA/netwk13/netwk13_module_df.rds")

## extract subnetworks
### create function to extract subnetworks by module ID
xtrxt_net <- function(x){  # put in desired module in quotes
  graph_from_data_frame(net_edges[((net_edges$fromNode %in% module_df$gene_id[module_df$colors == x]) &
                                     (net_edges$toNode %in% module_df$gene_id[module_df$colors == x])),],
                        vertices = net_nodes[net_nodes$module== x,])
}
####  output co-expression hub genes ##################
## define function to output co-expression hub genes
OutputHubs <- function(x){
  index <- get(paste0("egnvctr_", x))$gene %in% summary_df$ppep.site
  gn <- get(paste0("egnvctr_", x))$gene
  egncntrlt <- get(paste0("egnvctr_", x))$egnvctr
  hubs <- merge(data.frame(ppep.site = gn, eigencentrality = egncntrlt,hub=x),summary_df,"ppep.site")
}
# thresholding plots
## define plotting function
thrsh_plt <- function(x){ 
  ggplot()+geom_point(data = mod1, aes(x = nmbr, y = egnvctr))+
    geom_vline(xintercept = as.numeric(brkpnt), color='red')+labs(title = x)+theme_classic()}
########## black module ############ 
### run subnetwork extraction
net_bk <- xtrxt_net("black")
# calculate eigencentrality
egnvctr_bk <- eigen_centrality(net_bk)
egnvctr_bk <- data.frame(egnvctr = unname(egnvctr_bk[["vector"]]),gene = names(egnvctr_bk[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_bk$egnvctr),nmbr = seq_along(1:length(egnvctr_bk$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 76(the 2nd brkpnt)
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.3457546
thrsh_plt("black module") 
egnvctr_bk <- egnvctr_bk[egnvctr_bk$egnvctr > thresh,]

########## turquoise module ############ 
net_trqs <- xtrxt_net("turquoise")
egnvctr_trqs <- eigen_centrality(net_trqs)
egnvctr_trqs <- data.frame(egnvctr = unname(egnvctr_trqs[["vector"]]),gene = names(egnvctr_trqs[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_trqs$egnvctr),nmbr = seq_along(1:length(egnvctr_trqs$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 448
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.5160207
thrsh_plt("turquoise module") 
egnvctr_trqs <- egnvctr_trqs[egnvctr_trqs$egnvctr > thresh,]
########## brown module ############ 
net_brw <- xtrxt_net("brown")
egnvctr_brw <- eigen_centrality(net_brw)
egnvctr_brw  <- data.frame(egnvctr = unname(egnvctr_brw[["vector"]]),gene = names(egnvctr_brw[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_brw$egnvctr),nmbr = seq_along(1:length(egnvctr_brw$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 215
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.6447883
thrsh_plt("brown module") 
egnvctr_brw <- egnvctr_brw[egnvctr_brw$egnvctr > thresh,]
########## yellow module ############ 
net_yllw <- xtrxt_net("yellow")
egnvctr_yllw <- eigen_centrality(net_yllw)
egnvctr_yllw <- data.frame(egnvctr = unname(egnvctr_yllw[["vector"]]),gene = names(egnvctr_yllw[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_yllw$egnvctr),nmbr = seq_along(1:length(egnvctr_yllw$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 158
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.4338487
thrsh_plt("yellow module") 
egnvctr_yllw <- egnvctr_yllw[egnvctr_yllw$egnvctr > thresh,]
########## red module ############ 
net_red <- xtrxt_net("red")
egnvctr_red <- eigen_centrality(net_red)
egnvctr_red <- data.frame(egnvctr = unname(egnvctr_red[["vector"]]),gene = names(egnvctr_red[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_red$egnvctr),nmbr = seq_along(1:length(egnvctr_red$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 179
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.4784781
thrsh_plt("red module") 
egnvctr_red <- egnvctr_red[egnvctr_red$egnvctr > thresh,]
########## blue module ############ 
net_blue <- xtrxt_net("blue")
egnvctr_blue <- eigen_centrality(net_blue)
egnvctr_blue <- data.frame(egnvctr = unname(egnvctr_blue[["vector"]]),gene = names(egnvctr_blue[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_blue$egnvctr),nmbr = seq_along(1:length(egnvctr_blue$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 269
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.5909727
thrsh_plt("blue module") 
egnvctr_blue <- egnvctr_blue[egnvctr_blue$egnvctr > thresh,]
########## green module ############ 
net_grn <- xtrxt_net("green")
egnvctr_grn <- eigen_centrality(net_grn)
egnvctr_grn <- data.frame(egnvctr = unname(egnvctr_grn[["vector"]]),gene = names(egnvctr_grn[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_grn$egnvctr),nmbr = seq_along(1:length(egnvctr_grn$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2)
brkpnt <- mod3[["psi"]][[2, 2]] %>% round() # 212
thresh <- mod1[as.numeric(brkpnt),][[1]]  # 0.6794623
thrsh_plt("green module") 
egnvctr_grn <- egnvctr_grn[egnvctr_grn$egnvctr > thresh,]
########## grey module ############ 
net_grey <- xtrxt_net("grey")
egnvctr_grey <- eigen_centrality(net_grey)
egnvctr_grey <- data.frame(egnvctr = unname(egnvctr_grey[["vector"]]),gene = names(egnvctr_grey[["vector"]]))
mod1 <- data.frame(egnvctr = sort(egnvctr_grey$egnvctr),nmbr = seq_along(1:length(egnvctr_grey$egnvctr)))
mod2 <- glm(egnvctr ~ nmbr, data = mod1)  
mod3 <- segmented(mod2, seg.Z =  ~ nmbr, npsi = 2) #Breakpoint estimate(s) outdistanced to allow finite estimates and st.errs
brkpnt <- mod3[["psi"]][[2, 2]]  # null
#thresh <- mod1[as.numeric(brkpnt),][[1]]  # 6.020654e-15
#thrsh_plt("grey module") 
egnvctr_grey <- egnvctr_grey[egnvctr_grey$egnvctr > 0.416,] #Pick top3

## execute output function
hubs_black <-OutputHubs("bk") # 41
hubs_trqs <- OutputHubs("trqs") # 100
hubs_brw <- OutputHubs("brw") # 20
hubs_yllw <- OutputHubs("yllw") # 74
hubs_red <- OutputHubs("red") # 15
hubs_blue <- OutputHubs("blue") # 33
hubs_green <- OutputHubs("grn") # 14
hubs_grey <- OutputHubs("grey") # 3
hub <- rbind(hubs_black,hubs_trqs,hubs_brw,hubs_yllw,hubs_red,hubs_blue,hubs_green,hubs_grey) # 300 hubs
table(hub$hub)
rownames(hub) <- hub$ppep.site
#saveRDS("../SFCH/Data/WGCNA/netwk13/hub.rds",object=hub)

#### WGCNA:hubs check #####
hub <- readRDS("../SFCH/Data/WGCNA/netwk13/hub.rds")
table(hub$hub)

###################### history #################################
####### exportNetworkToCytoscape ##########
# eigencentrality?? remove edges for substrate-substrate only
Node <-net_nodes%>%mutate(hub=hub[nodeName,"hub"],summary_df[nodeName,])%>%dplyr::select(-c(2,5,7,12:17,32:195))
Edge <-net_edges 

write.table(Node, file = "../SFCH/Result/WGCNA/Cytoscape/Node.txt", sep = "\t",row.names = F, col.names = TRUE,quote = F)
write.table(Edge, file = "../SFCH/Result/WGCNA/Cytoscape/Edge.txt", sep = "\t",row.names = F, col.names = TRUE,quote = F)


