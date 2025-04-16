rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(factoextra)
library(ggplot2)
library(tidyverse)
library(stringr)
library(reshape2)# melt
library(Biobase)
library(Mfuzz)

###### determine the optimal K: 7 clusters #########
imp30_kmc_df<- readRDS("../SFCH/Data/clustering/Newimp30_kmc_df_PS2+MF1.rds")
df <- t(scale(t(na.omit(imp30_kmc_df[,2:7]))))%>%as.data.frame()
set.seed(123)
fviz_nbclust(x = df,FUNcluster = kmeans,method = "gap_stat",diss =get_dist(df,method ="pearson"),nboot =100,verbose = TRUE,k.max =10)+labs(subtitle = "Gap Statistic Method")
#ggsave(filename = "../SFCH/Figure/Diss_chap1/Gap Statistic Method_imp7_PS2MF1.pdf",width =5,height =4)
dev.off()
########## Mfuzz : soft clustering ########
#Assay data
imp30_kmc_df<- readRDS("../SFCH/Data/clustering/Newimp30_kmc_df_PS2+MF1_threshold.rds")
df <- imp30_kmc_df[,2:7]%>%as.matrix()
minimalSet <- ExpressionSet(assayData=df)
#feature data:The number of rows in featureData must match the number of rows in assayData. Row names of featureData must match row names of the matrix / matricies in assayData.
featureData <- data.frame(labelDescription=word(rownames(minimalSet@assayData[["exprs"]]),1,sep = fixed("_")),
                          treatment=word(rownames(minimalSet@assayData[["exprs"]]),2,sep = fixed("_")),
                          row.names=rownames(minimalSet@assayData[["exprs"]]))
table(featureData$treatment)
exampleSet <- ExpressionSet(assayData=df,featureData = new("AnnotatedDataFrame",data=featureData))
dfSet <- exampleSet%>%filter.NA()%>%standardise()
set.seed(123)
cl <- mfuzz(dfSet,c=7,m = 1.5,iter.max = 100)
cl$size
#saveRDS("../SFCH/Data/clustering/SFCH_kmc7_PS2+MF1_threshold_mfuzz_cl.rds",object = cl)

cl <- readRDS("../SFCH/Data/clustering/SFCH_kmc7_PS2+MF1_threshold_mfuzz_cl.rds")

####### arrange 7 clusters:P_drop according to time point ######
gg_Max <- data.frame(cl$membership,membership=apply(cl$membership,1,max),cluster=apply(cl$membership,1,which.max))%>%mutate(ppep.site=word(rownames(.),1,sep = fixed("_")),.before=1)
a<- as.data.frame(df)%>%mutate(membership=gg_Max$membership,cluster= as.character(cl$cluster),treatment=str_split_fixed(rownames(df),"_",2)[,2])%>%mutate(ppep.site=word(rownames(.),1,sep = fixed("_")),.before=1)
a$cluster <- factor(a$cluster, levels = c("1","2","3","4","5","6","7"))
a$treatment <- factor(a$treatment, levels = c("sys","flg22","chitin","H2O"))
a$cluster <- factor(a$cluster, levels =c("2","4","3","5","7","1","6"),labels = c("Cluster 1", "Cluster 2", "Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7"))
#saveRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg.rds",object = a)
Order_a <- a[order(a$cluster,a$treatment),]
#saveRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order.rds",object = Order_a)

####### dataframe prep: Order_a and distribution #######
Order_a <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order.rds")
distribution <- Order_a%>%group_by(cluster=as.character(word(cluster,2)))%>%dplyr::select(cluster,treatment,ppep.site)%>%
  pivot_wider(names_from=treatment, values_from=cluster)%>%mutate(Protein.ID=substr(ppep.site,1,18),.before=1)%>%as.data.frame()
rownames(distribution) <- distribution$ppep.site
#saveRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds",object=distribution)

######## plot: seven P profile arranged by the 1st drop ##########
a <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg.rds")
b <- melt(a%>%select(-membership),id.vars = c("cluster","treatment","ppep.site"),variable.name = c("time_point"), value.name = "norm_intensity")
se <- function(x) sqrt(var(x)/length(x))#SE calculation
Percluster <- b %>% group_by(cluster, treatment, time_point) %>% summarize(mean=mean(norm_intensity), SE=se(norm_intensity))
ggplot(data = b, aes(x = time_point, y = norm_intensity,group=ppep.site))+labs(y= "Normalized intensity", x = "Time (min)")+
  theme(legend.position = "none")+scale_color_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"),labels=c("Systemin","Flg22","Chitin","H2O"))+
  geom_line(aes(color=treatment),linewidth =0.6,alpha=0.02)+
  facet_grid(as.factor(treatment)~as.factor(cluster))+
  geom_line(Percluster,mapping=aes(y =mean,group=interaction(cluster,treatment)),alpha=0.5)+
  geom_errorbar(Percluster,mapping=aes(x = time_point,y =mean,ymin=mean-SE, ymax=mean+SE,group=interaction(cluster,treatment)),width=0.2)+
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
  geom_text(data=plyr::count(a, vars =c("treatment","cluster")),aes(x=1.8, y=1.8, label=paste0("N=",freq)),colour="black", inherit.aes=FALSE, parse=FALSE,size = 3.5)
#ggsave("../SFCH/Figure/Diss_chap1/PS2+MF1_threshold_Mfuzz_k7.pdf",width =9,height =6)
#ggsave("../SFCH/Figure/Diss_chap1/PS2+MF1_threshold_Mfuzz_k7_sys.pdf",width =9,height =6)
dev.off()

######## stacking plot: P_sites distribution ##########
#Rm water do not display water... #
a <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg.rds")%>%filter(treatment!="H2O")
## stacking plot ##
count_by_treat <- NULL
for (i in levels(a$cluster)){num <-data.frame(a%>%filter(cluster==i)%>%count(treatment),cluster=as.character(word(i,2)))
                             count_by_treat <- rbind(count_by_treat,num)}
count_by_treat <- merge(count_by_treat,count_by_treat%>%group_by(cluster)%>%summarise(sum=sum(n)))
## y percentage 
ggplot(count_by_treat, aes(fill=treatment, y=n/sum, x=cluster)) + geom_bar(position="stack", stat="identity",show.legend = T,alpha=0.8)+
  ggtitle(" ")+ylab("P_sites %")+guides(fill=guide_legend(title="Treatment"))+
  scale_fill_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"),labels=c("Systemin","Flg22","Chitin","H2O"))+
  geom_text(aes(label=n),size=4,color="white",position = position_stack(vjust = 0.5))+#theme_minimal()+
  theme(plot.title = element_text(colour = "black",size=12),
        axis.title=element_text(size=12,colour = "black"),
        axis.text = element_text(colour = "black",size=10),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_text(size=12),
        #legend.position="top",
        legend.text=element_text(size=10))

#### go to L79 to make a stacking plot ######
#ggsave("../SFCH/Figure/Diss_chap1/k7_imp30_stack_sys_mfuzz(no water).pdf",width =4.8,height =4)

#### go to L47 to plot the temporal profile ######
dev.off()

### remove shared P_sites among SFC ##########
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")%>%filter(sys!=flg22,sys!=chitin)
Order_a_rem <- Order_a[Order_a$ppep.site%in%rownames(distribution),]
int_norm <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
a <- Order_a_rem[Order_a_rem$ppep.site%in%rownames(int_norm),]
#Rm water:do not display water
a <- a%>%filter(treatment!="H2O")
#ggsave("../SFCH/Figure/Diss_chap1/k7_imp30_stack_sys_mfuzz(no water).pdf",width =4.8,height =4)





