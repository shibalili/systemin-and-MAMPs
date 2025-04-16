rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)

############## extract subnetwork ######################
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary4701_df_color.rds")
# !!! as.character(regulatoryGene) because they are factors!!!
grn_edges_all<- readRDS("../SFCH/Data/WGCNA/genie3/grn_edges_all.rds") 
class(grn_edges_all$regulatoryGene) ## factor!!! 

grn_nodes_all<- readRDS("../SFCH/Data/WGCNA/genie3/grn_nodes_all.rds")
class(grn_nodes_all$ppep.site)



############## PLL1a centered subnetwork ######################
# remove PLL1a-PLL1b edges
E_sub <- grn_edges_all%>%filter_at(vars(regulatoryGene,targetGene),any_vars(str_detect(., pattern = "Solyc08g077150.3.1")))%>%
  filter(tarProtein!="Solyc08g007000.3.1",regProtein!="Solyc08g007000.3.1") # 1068
class(E_sub$regulatoryGene) # !!! factor
E_node <- grn_nodes_all[as.character(unique(c(E_sub$regulatoryGene,E_sub$targetGene))),]%>%filter(super_class!=0) # 256
table(E_node$super_class)

# output GRN network
#write_delim(E_sub, file = "../SFCH/Data/WGCNA/genie3/grn_network_edges_PLL1a.txt", delim = "\t")
#write_delim(E_node, file = "../SFCH/Data/WGCNA/genie3/grn_network_nodes_PLL1a.txt", delim = "\t")

# ROS #
E_node_ros <- E_node%>%filter(super_class=="ROS")
table(E_node_ros$Protein.ID)
#knitr::kable(table(E_node_ros$Protein_class), format = "simple")
E_PLL1_ros <- E_sub%>%filter_at(vars(regProtein_class,tarProtein_class),any_vars(str_detect(., pattern = "ROS")))


# Kinase #
E_node_kins <- E_node%>%filter(super_class=="Kinase")
#knitr::kable(table(E_node_kins$Protein_class), format = "simple")
E_PLL1_pbl <- E_sub%>%filter_at(vars(regProtein_class,tarProtein_class),any_vars(str_detect(., pattern = "PBL")))
E_PLL1_BIK <- E_sub%>%filter_at(vars(regProtein_class,tarProtein_class),any_vars(str_detect(., pattern = "BIK")))

# Ca2+  #
E_node_ca <- E_node%>%filter(super_class=="Ca2+")
knitr::kable(table(E_node_ca$Protein_class), format = "simple")
E_PLL1_CPK <- E_sub%>%filter_at(vars(regProtein_class,tarProtein_class),any_vars(str_detect(., pattern = "CPK")))


###### add edge annotation to subnetwork of PLL1a ###########
library(openxlsx)
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df_Net.rds")

my_data <- read.table(pipe("pbpaste"),sep = "\t",header = TRUE)
my_data$regProtein <-summary_df[my_data$regulatorySite,"Protein_class"] 
my_data$tarProtein <-summary_df[my_data$targetSite,"Protein_class"] 
write.xlsx("../SFCH/Data/WGCNA/genie3/KPnet_edgePLL1a_1681_anno.xlsx",x =my_data,rowNames=F,overwrite = T)














