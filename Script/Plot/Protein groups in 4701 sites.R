rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)
# Q: how many protein groups in this dataset?
# Q: Is specificity achieved at Protein groups level?
# input raw by presence only  #
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
count <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%mutate(ppep.site=rownames(.),.before = 1)%>%mutate(Protein.ID=substr(rownames(.),1,18),.before = 1)
sp_s <-count%>%filter(Count_systemin!=0)%>%mutate(attr="systemin",summary_df[ppep.site,c("super_class1","Protein_class")],.before = Count_Sum)
sp_f <- count%>%filter(Count_flg22!=0)%>%mutate(attr="flg22",summary_df[ppep.site,c("super_class1","Protein_class")],.before = Count_Sum)
sp_chi<- count%>%filter(Count_chitin!=0)%>%mutate(attr="chitin",summary_df[ppep.site,c("super_class1","Protein_class")],.before = Count_Sum)
sp_w <- count%>%filter(Count_H2O!=0)%>%mutate(attr="H2O",summary_df[ppep.site,c("super_class1","Protein_class")],.before = Count_Sum)
sp <- rbind(sp_s,sp_f,sp_chi,sp_w)%>%filter(!super_class1%in%c("0","PPP:PP4"))
table(sp$super_class1)
rownames(sp) <- NULL
sp$attr <- factor(sp$attr,levels = c("systemin","flg22","chitin","H2O"))
gr <- plyr::count(plyr::count(sp, vars =c("super_class1","Protein.ID","attr"))%>%select(1:3),vars =c("super_class1","attr"))

# order protein groups 
unique(gr$super_class1)[c(25:29)] # Kinase
unique(gr$super_class1)[c(21:24,2)] # PPase
unique(gr$super_class1)[c(30)] # ROS
unique(gr$super_class1)[c(10,7:8,6,9,12,11)]# Ca2+
unique(gr$super_class1)[c(19,18,17,20)] # MAPK
unique(gr$super_class1)[c(16,34:40)] # Transporter
unique(gr$super_class1)[c(32,31,33)]# TF
unique(gr$super_class1)[c(1,3,5,4,13:14)]# hormone
unique(gr$super_class1)[c(15)]# ETI

gr$super_class1 <- factor(gr$super_class ,levels = c(unique(gr$super_class)[c(25:29)],unique(gr$super_class)[c(21:24,2)],
                                                     unique(gr$super_class)[30],unique(gr$super_class)[c(10,7:8,6,9,12,11)],
                                                     unique(gr$super_class)[c(19,18,17,20)],unique(gr$super_class)[c(16,34:40)],
                                                     unique(gr$super_class)[c(32,31,33)],unique(gr$super_class)[c(1,3,5,4,13:14)],
                                                     unique(gr$super_class)[15]))



## stack bar plot ###
ggplot(gr, aes(fill=attr, y=super_class1, x=freq)) + geom_bar(position="stack", stat="identity",show.legend = T,alpha=0.8)+scale_y_discrete(limits=rev)+
  ggtitle(" ")+xlab("Protein numbers")+ylab("")+guides(fill=guide_legend(title="Treatment"))+facet_grid(~as.factor(attr))+
  scale_fill_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"),labels=c("Systemin","Flg22","Chitin","H2O"))+
  geom_text(aes(label=freq),size=3,color="white",position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(colour = "black",size=12),
        axis.title=element_text(size=12,colour = "black"),
        axis.text = element_text(colour = "black",size=10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        legend.background =  element_rect(fill = "transparent",color = NA),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10))
ggsave("../SFCH/Figure/Diss_chap1/Raw_Protein_groups.pdf",width =7,height =6)

dev.off()

