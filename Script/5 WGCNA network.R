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
registerDoParallel(cores=4) # Warning: executing %dopar% sequentially: no parallel backend registered 

## scaling for P_sites time profile (P develop after time 0; so scaling = blank with time 0)
Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%dplyr::select(2:25)%>%t()

### Pick a soft threshold power near the curve of the plot ###
plot_stps <- function(x){
  powers = c(1:20) # Choose a set of soft-thresholding powers
  sft = pickSoftThreshold(x, powerVector = powers, verbose = 5) # Call the network topology analysis function
  par(mfrow = c(1,2)) ; cex1 = 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1],sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
}
plot_stps(Input) # power 13(scale)
# blockwiseModules: network construction and module detection in a block-wise manner
blckws_mdls <- function(x){
  picked_power = 13
  temp_cor <- cor
  cor <- WGCNA::cor
  allowWGCNAThreads() # allows multithreading
  y <- blockwiseModules(x,# input
                        # == Adjacency Function ==
                        power = picked_power,  # soft-thresholding power for network construction
                        networkType = "signed",  # allows module membership for negatively-correlated genes
                        corType = "pearson",
                        # == Tree and Block Options ==
                        deepSplit = 0,  # module splitting sensitivity
                        pamRespectsDendro = F,
                        minModuleSize = 30,
                        maxBlockSize = 5000,
                        # == Module Adjustments ==
                        mergeCutHeight = 0.37,  # threshold correlation to merge modules
                        # == TOM == Archive the run results in TOM file (saves time)
                        saveTOMs = T,
                        saveTOMFileBase = "WGCNA", # pls rename TOM (WGCNA13-block.1.RData) to avoid overiding 
                        # == Output Options
                        numericLabels = T,
                        verbose = 3)
  cor <- temp_cor
  return(y)
}
netwk<- blckws_mdls(Input) 
#save(netwk, file = "../SFCH/Data/WGCNA/netwk13/netwk_13.Rdata")


## load data
load("../SFCH/Data/WGCNA/netwk13/netwk_13.Rdata")

# data frame of module assignments for all genes
module_df <- data.frame(gene_id = names(netwk$colors), colors = labels2colors(netwk$colors))
rownames(module_df) <- module_df$gene_id
#saveRDS("../SFCH/Data/WGCNA/netwk13/netwk13_module_df.rds",object=module_df) #update:27.03.2025

module_df <- readRDS("../SFCH/Data/WGCNA/netwk13/netwk13_module_df.rds")
sort(table(module_df$colors),decreasing = T)
# grey: not part of any module
# turquoise:the largest module, next comes blue, next brown, etc.
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
table(summary_df$super_class)

# plot dendrogram
# The heights relate to the distance metric which is  One minus Pearson correlation (1 - r)
# Samples at the very bottom are weakly correlated to all of those in the tree structure at the top 
mergedColors = labels2colors(netwk$colors)  # Convert labels to colors for plotting
table(mergedColors)
plot_dendro <- function(x){
  plotDendroAndColors(x$dendrograms[[1]],mergedColors[x$blockGenes[[1]]],#"module colors"
                      dendroLabels = FALSE,hang = 0.03,addGuide = TRUE,guideHang = 0.05,
                      theme(plot.title = element_blank()))}
plot_dendro(netwk)
ggsave("../SFCH/Figure/Diss_chap1/WGCNA_8modules_dendro.pdf",width = 6,height =4)


# To figure out which modules are associated with each trait (treatment x time) group
MEs0 <- moduleEigengenes(Input, module_df$colors)$eigengenes  # data frame of module eigengenes (1st principal component)
MEs0 <- orderMEs(MEs0)  # reorder modules so similar modules are next to each other
module_order = names(MEs0) %>% gsub("ME","", .)  # create vector of module order

#module_order <- c("turquoise","red","black","brown","blue","yellow","grey","green")
module_order <- c("yellow","red","black","blue","turquoise","brown","green","grey")


MEs0$sample = row.names(MEs0)  # add column of treatment ID
# mME: the long dataframe of MEs0 for [treatment * time point]
mME = MEs0%>%pivot_longer(-sample)%>%mutate(name = gsub("ME", "", name),name = factor(name, levels = module_order))
setup<-readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")%>%mutate(name=paste0(treatment,time.point))
mME$treat <- substr(mME$sample,1,1)
mME$treat <- factor(mME$treat,levels = c("S","F","C","H"))
mME$time <- substr(mME$sample,2,3)
mME$time <- factor(mME$time,levels = c("0","1","2","5","15","45"))
mME$sample <- factor(mME$sample,levels = unique(mME$sample))
ggplot(mME, aes(x=time, y=name, fill=value))+facet_grid(.~treat)+scale_y_discrete(limits=rev)+
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0,limit = c(-1,1)) +
  geom_tile() +theme_bw()+theme(axis.text.x = element_text(angle=90)) +labs(title = "", x = "", y = "", fill="Correlation")
#ggsave("../SFCH/Figure/Diss_chap1/WGCNA_8modules.pdf",width = 6,height = 2)
dev.off()



# Pull out list of P_sites in that module or all P_sites profiles in each module
# pick out a few modules of interest here
#modules_of_interest = c("black","blue",'green',"brown","red","yellow","turquoise","grey")
modules_of_interest <- labels2colors(netwk$colors)%>%unique()
submod <- module_df %>%subset(colors %in% modules_of_interest)# subset module of interest
df <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%dplyr::select(2:25) # t(Input)
subexpr <- df[submod$gene_id,]
submod_df <-  data.frame(subexpr)%>%mutate(gene_id = row.names(.))%>%pivot_longer(-gene_id)%>%
  mutate(module = module_df[gene_id,]$colors,treat=substr(name,1,1),time=substr(name,2,3))
submod_df$treat <- factor(submod_df$treat,levels = c("S","F","C","H"))
submod_df$time <- factor(submod_df$time,levels = c("0","1","2","5","15","45"))
submod_df$name <- factor(submod_df$name,levels = unique(submod_df$name))
submod_df$module <- factor(submod_df$module,levels = module_order)
## simple plot ##
ggplot(submod_df, aes(x=time, y=value, group=gene_id)) +scale_color_manual(values= module_order,labels=module_order)+
  geom_line(aes(color = module),alpha = 0.05) +theme_bw() +theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module),cols = vars(treat)) +labs(x = "treatment",y = "normalized expression")
dev.off()
## add mean line and SD
se <- function(x) sqrt(var(x)/length(x))#SE calculation
Percluster <- submod_df %>% group_by(module, treat, time) %>% summarize(mean=mean(value), SE=se(value))
ggplot(data = submod_df, aes(x = time, y = value,group=gene_id))+labs(y= "Normalized intensity", x = "Time (min)")+
  theme(legend.position = "none")+scale_color_manual(values= module_order,labels=module_order)+
  geom_line(aes(color=module),linewidth =0.6,alpha=0.02)+
  facet_grid(as.factor(module)~as.factor(treat))+
  geom_line(Percluster,mapping=aes(y =mean,group=interaction(module,treat)),alpha=0.5)+
  geom_errorbar(Percluster,mapping=aes(x = time,y =mean,ymin=mean-SE, ymax=mean+SE,group=interaction(module,treat)),width=0.2)+
  theme(axis.text = element_text(colour = "black",size=13),
        axis.title=element_text(size=16,colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.text= element_text(size =14,colour = "black"),
        panel.grid = element_blank(),strip.background = element_blank(),
        #legend.position="top",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_blank(),
        legend.text=element_blank())+
  geom_text(data=plyr::count(submod_df, vars =c("treat","module")),
            aes(x=1.8, y=1.8, label=paste0("N=",freq/6)),colour="black", inherit.aes=FALSE, parse=FALSE,size = 3.5)
#sort(table(module_df$colors),decreasing = T)

#### Generate and Export Networks #####
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary4701_df_color.rds")
load("../SFCH/Data/WGCNA/netwk13/netwk_13.Rdata") # load netwk
load("../SFCH/Data/WGCNA/netwk13/WGCNA13-block.1.RData") # load TOM

# data frame of module assignments for all genes
module_df <- data.frame(gene_id = names(netwk$colors), colors = labels2colors(netwk$colors),Protein_class=summary_df[names(netwk$colors),"Protein_class"])
table(module_df$colors)
rownames(module_df) <- module_df$gene_id
modules_of_interest =labels2colors(netwk$colors)%>%unique()
genes_of_interest = module_df %>% subset(colors %in% modules_of_interest)

#### The network was  filtered for low edge weight (edge weight < 0.154) using the R package ‘segmented’ before visualization with Cytoscape.
## calculate edge weight threshold
temp1 <- exportNetworkToCytoscape(TOM,
                                  edgeFile = NULL,
                                  nodeFile = NULL,
                                  weighted = TRUE,
                                  threshold = 0.01,
                                  nodeNames = rownames(genes_of_interest),
                                  altNodeNames = genes_of_interest$Protein_class,
                                  nodeAttr = genes_of_interest$colors,
                                  includeColNames = TRUE)

temp2 <- temp1$edgeData$weight[seq(1, length(temp1$edgeData$weight),100)]

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
brkpnt_fun(temp2, 2) #"Breakpoints are 4028, 5286"
sort(temp2) %>% .[4028] # edge weight cutoff calculated to 0.03560391


temp1 <- exportNetworkToCytoscape(TOM,edgeFile = NULL,nodeFile = NULL,weighted = TRUE,threshold = 0.03560391,
                                  nodeNames = rownames(genes_of_interest),altNodeNames = genes_of_interest$Protein_class,
                                  nodeAttr = genes_of_interest$colors,includeColNames = TRUE)

###### co-expression hub gene prediction ####################################################################################
# loop removed by (exportNetworkToCytoscape)
# remove P_sites from same protein
# add color
net_edges <- temp1$edgeData%>%mutate(module1=module_df[fromNode,"colors"],module2=module_df[toNode,"colors"],protein1=substr(fromNode,1,18),protein2=substr(toNode,1,18))%>%filter(protein1!=protein2)
net_nodes <- temp1$nodeData%>%filter(nodeName%in%unique(c(net_edges$fromNode,net_edges$toNode)))%>%mutate(module=module_df[nodeName,"colors"])
rownames(net_nodes) <- net_nodes$nodeName
#saveRDS("../SFCH/Data/WGCNA/netwk13/net_edges.rds",object = net_edges)
#saveRDS("../SFCH/Data/WGCNA/netwk13/net_nodes.rds",object = net_nodes)

