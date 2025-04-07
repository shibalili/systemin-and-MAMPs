rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)
library(ggpubr)

summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")%>%filter(Count_Sum>30)
##### same protein group,different proteins, each specific to 1 treatment #########
# protein group:PBLs 
pbl <- summary_df%>%filter(str_detect(super_class1,"PBL"))
length(unique(pbl$Protein.ID)) #15
plyr::count(pbl, vars =c("Protein_class","Protein.ID"))%>%arrange(desc(freq))
# PBL30 (sys_sp):Solyc03g032150.3.1 251; BIK1 Solyc06g005500.3.1 25(a complete profile)
# PBL11 (flg22_sp):Solyc10g084770.2.1 19; maybe PBL7 Solyc10g085990.2.1 374 is better
# PBL19 (chitin_sp):Solyc01g088690.3.1 18
# PBL3 (sfc_shared) :Solyc05g007140.3.1 67; PBL8 Solyc10g074710.2.1 46(a complete profile)
# SlPBS1 (sfc_sp) : Solyc05g024290.3.1 105 (a complete profile or box)




##### one site, different kinetics #########
# a transporter: SLAC (impute data) Solyc09g014610.3.1 51
slac <- summary_df%>%filter(str_detect(super_class1,"SLAC"))


##### one proteins, different sites, different kinetics #########
# PPase:PLL1a
pll <- summary_df%>%filter(str_detect(Protein_class,"PLL1a"))




setup_combine <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
int_combine <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%dplyr::select(7:138)
#int_combine <-readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")%>%dplyr::select(c(2:133))
#int_combine <- readRDS("../SFCH/Data/imputation/SFCH_imp30.rds")%>%select(3:134)#imp30 MissForest

df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
#ID <-"Solyc03g032150.3.1 251"
#ID <-"Solyc10g084770.2.1 19" #scale_y_continuous(limits=c(20, 23.5))+
#ID <-"Solyc10g085990.2.1 374" #scale_y_continuous(limits=c(20.5, 23.5))+
#ID <-"Solyc01g088690.3.1 18" #scale_y_continuous(limits=c(21.5, 24.5))+
#ID <-"Solyc05g007140.3.1 67" #scale_y_continuous(limits=c(21, 24.5))+
#ID <- "Solyc05g024290.3.1 105"
#ID <- "Solyc08g077150.3.1 89" #scale_y_continuous(limits=c(21.9, 23.3))+
#ID <- "Solyc08g077150.3.1 190" # seems to be sys_sp
ID <- "Solyc08g007000.3.1 251" # PLL1b but this site is conserved at PLL1a WT-WT-A/D #scale_y_continuous(limits=c(20.6, 24.2))+

ID2 <- unique(df[ID,"Protein_class"])
plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                        intensity = as.numeric(int_combine[ID,]),
                        treatment = setup_combine[,"treatment"])#%>%filter(treatment=="F")
group.colors <- c(S = "#E41A1C", F = "#377EB8", C ="#4DAF4A", H = "#E7B800")
ggboxplot(plot.data, x = "minutes", y = "intensity",color = "treatment",add = ("jitter"),shape = "treatment",bxp.errorbar=T,
          title = ID,subtitle =ID2,xlab = c("Time (min)"),ylab ="Nomalized intensity",na.rm=T)+scale_y_continuous(limits=c(20.6, 24.2))+
  facet_wrap(vars(treatment),nrow = 1)+
  scale_color_manual(values = group.colors)+scale_linetype_manual(values=c("dashed"))+
  scale_x_discrete(breaks=0:5, labels = c("0","1","2","5","15","45"))+
  stat_summary(aes(x = minutes, y = intensity, color = treatment, group = treatment),fun=mean, geom="line",na.rm = TRUE,show.legend = F,linewidth =0.2,alpha=0.8)+
  theme(axis.text = element_text(size=10,color = "black"),title = element_text(size=10,color = "black"),
        plot.title = element_text(size=10,color = "black",hjust = 0),
        panel.background = element_rect(fill = "transparent"),strip.text = element_text(size=16,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background =element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),
        legend.key = element_blank(),legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.text=element_text(size=10))+guides(size=F,shape=F,color=F)
ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/",ID,"_box.pdf",sep=""), width =5,height =2.5)



#ID <- "Solyc06g005500.3.1 25"
#ID <- "Solyc10g074710.2.1 46"
#ID <- "Solyc05g024290.3.1 105"

#ID <-"Solyc09g014610.3.1 51"  # SLAC

#ID <- "Solyc08g077150.3.1 50" # 50-54-59 flg22_sp # A/D-WT-WT
#ID <- "Solyc08g077150.3.1 181" # 163;172-181 sys_sp # WT-A/D-WT
#ID <- "Solyc08g077150.3.1 86;88" # chitin_sp:86;88 and 89 # not yet mutate
#ID <- "Solyc08g007000.3.1 251" # PLL1b but this site is conserved at PLL1a WT-WT-A/D

#ID <- "Solyc08g077150.3.1 163;172"



############ ROS ###############
####### how many RBOHs in the dataset #########
rboh <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df_Net.rds")%>%filter(super_class=="ROS")
table(rboh$Protein.ID,rboh$Protein_class)

# 
# Solyc03g117980.3.1 340
# Solyc03g117980.3.1 33
# Solyc01g099620.3.1 648
# Solyc01g099620.3.1 18
# Solyc01g099620.3.1 291/287
# Solyc06g068680.3.1 15

ID <- "Solyc06g068680.3.1 15"
ID <- "Solyc01g099620.3.1 18"
ID <- "Solyc03g117980.3.1 116"
ID <- "Solyc01g099620.3.1 648"
ID <- "Solyc03g117980.3.1 340"
ID <- "Solyc01g099620.3.1 291"
ID <- "Solyc01g099620.3.1 287"


####### PLL1a-PBL/BIK-RBOH ############
ID <- "Solyc04g011520.3.1 13"
ID <- "Solyc09g010850.3.1 388"
ID <- "Solyc10g084770.2.1 19"
ID <- "Solyc05g007140.3.1 51" # ??

ID <- "Solyc04g011520.3.1 20" 
ID <- "Solyc06g005500.3.1 13"


ID <- "Solyc10g074710.2.1 49"


####### PLL1a-CPK-RBOH  ##############
ID <- "Solyc04g009800.3.1 572"
ID <- "Solyc11g065660.2.1 524"
ID <- "Solyc10g078390.2.1 54"
ID <- "Solyc07g064610.3.1 491"
ID <- "Solyc03g082500.3.1 344;348;347;363;383;331"
ID <- "Solyc01g010780.3.1 58"
ID <- "Solyc01g112250.3.1 523"


######### CRK ###########
ID <- "Solyc05g005070.2.1 385"



########## complete data ############## 
ID <- "Solyc03g123920.3.1 127" #132
########## MAR #############
ID <- "Solyc02g081810.3.1 441" #122
ID <- "Solyc05g005620.3.1 448" #115
ID <- "Solyc04g051310.3.1 512" #89
########## MNAR ###############
ID <- "Solyc03g120980.3.1 148"
ID <- "Solyc02g082690.3.1 404" #scale_y_continuous(limits = c(19.5,24.2))+
ID <- "Solyc06g065710.3.1 367"


ID2 <- unique(df[ID,"Protein_class"])
plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                        intensity = as.numeric(int_combine[ID,]),
                        treatment = setup_combine[,"treatment"])#%>%filter(treatment=="F")
ggplot(plot.data,aes(x = minutes, y = intensity, color = treatment, group = treatment))+facet_wrap(vars(treatment),nrow = 1)+
  geom_point(na.rm = TRUE)+stat_summary(fun=mean, geom="line",na.rm = TRUE)+
  scale_x_continuous(breaks=0:5, labels = c("0","1","2","5","15","45"))+scale_y_continuous(limits = c(19.5,24.2))+
  labs(title =ID,subtitle = ID2)+ylab(label = "Nomalized intensity")+xlab("Time (min)")+guides(size=F,shape=F,color=F)+
  scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"))+
  theme(axis.title = element_text(size=12,color = "black"),axis.line = element_line(colour = "black"),
        axis.text = element_text(size=8,color = "black"),
        plot.title = element_text(size=12,color = "black",hjust = 0),
        plot.subtitle = element_text(color = "black",size=10,hjust = 0),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))
#ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/",ID,"facet.pdf",sep=""), width =5,height =2.5)
#ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/",ID,".pdf",sep=""), width =3,height =3)
ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/Solyc02g082690.3.1 404/",ID,"_imp_mix.pdf",sep=""), width =5,height =2.5)
ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/Solyc02g082690.3.1 404/",ID,"_Ori.pdf",sep=""), width =5,height =2.5)
ggsave(paste("../SFCH/Figure/Diss_chap1/P_sites/Solyc02g082690.3.1 404/",ID,"_MFonly.pdf",sep=""), width =5,height =2.5)

dev.off()



################### history #####################

plot_Ori_Psite_line <- function(x){
  setup_combine <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
  int_combine <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%dplyr::select(7:138)
  df <- readRDS("../SFCH/Data/annotation/summary4701/SFCH_summary4701_PSMF_threshold_Mfuzz_protein_class_df.rds")
  plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                          intensity = as.numeric(int_combine[x,]),
                          treatment = setup_combine[,"treatment"])
  ggplot(plot.data,aes(x = minutes, y = intensity, color = treatment, group = treatment))+facet_wrap(vars(treatment),nrow = 1)+
    geom_point(na.rm = TRUE)+stat_summary(fun=mean, geom="line",na.rm = TRUE)+
    scale_x_continuous(breaks=0:5, labels = c("0","1","2","5","15","45"))+#scale_y_continuous(limits = c(22,27))+
    labs(title =x,subtitle = unique(df[x,"Protein_class"]))+ylab(label = "Phosphorylation intensity")+xlab("Time (min)")+
    scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"))+
    theme(axis.title = element_text(size=12,color = "black"),axis.line = element_line(colour = "black"),
          axis.text = element_text(size=8,color = "black"),
          plot.title = element_text(size=12,color = "black",hjust = 0),
          plot.subtitle = element_text(color = "black",size=10,hjust = 0),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          panel.grid.major = element_blank(),
          legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
          legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))
}
plot_imp_Psite_line <- function(x){
  setup_combine <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
  int_combine <-readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")%>%dplyr::select(c(2:133))
  df <- readRDS("../SFCH/Data/annotation/summary4701/SFCH_summary4701_PSMF_threshold_Mfuzz_protein_class_df.rds")
  plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                          intensity = as.numeric(int_combine[x,]),
                          treatment = setup_combine[,"treatment"])
  ggplot(plot.data,aes(x = minutes, y = intensity, color = treatment, group = treatment))+facet_wrap(vars(treatment),nrow = 1)+
    geom_point(na.rm = TRUE,alpha=0.7)+stat_summary(fun=mean, geom="line",na.rm = TRUE,alpha=0.7)+
    scale_x_continuous(breaks=0:5, labels = c("0","1","2","5","15","45"))+#scale_y_continuous(limits = c(22,27))+
    labs(title =x,subtitle = unique(df[x,"Protein_class"]))+ylab(label = "Phosphorylation intensity")+xlab("Time (min)")+
    scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"))+
    theme(axis.title = element_text(size=12,color = "grey"),axis.line = element_line(colour = "grey"),
          axis.text = element_text(size=8,color = "grey"),
          plot.title = element_text(size=12,color = "grey",hjust = 0),
          plot.subtitle = element_text(color = "grey",size=10,hjust = 0),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"),panel.border = element_rect(colour = "grey", fill=NA, linewidth=1),
          panel.grid.major = element_blank(),
          legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
          legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=10))
}

plot_Ori_Psite_box <- function(x){
  setup_combine <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
  int_combine <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%dplyr::select(7:138)
  df <- readRDS("../SFCH/Data/annotation/summary4701/SFCH_summary4701_PSMF_threshold_Mfuzz_protein_class_df.rds")
  plot.data <- data.frame(minutes = as.numeric(as.character(setup_combine$time.IDs)), 
                          intensity = as.numeric(int_combine[x,]),
                          treatment = setup_combine[,"treatment"])
  group.colors <- c(S = "#E41A1C", F = "#377EB8", C ="#4DAF4A", H = "#E7B800")
ggboxplot(plot.data, x = "minutes", y = "intensity",color = "treatment",add = ("jitter"),shape = "treatment",bxp.errorbar=T,
          title = x,subtitle =unique(df[x,"Protein_class"]),xlab = c("Time (min)"),ylab =F,na.rm=T)+
  facet_wrap(vars(treatment),nrow = 1)+
  scale_color_manual(values = group.colors)+scale_linetype_manual(values=c("dashed"))+
  scale_x_discrete(breaks=0:5, labels = c("0","1","2","5","15","45"))+
  stat_summary(aes(x = minutes, y = intensity, color = treatment, group = treatment),fun=mean, geom="line",na.rm = TRUE,show.legend = F,linewidth =0.2,alpha=0.8)+
  theme(axis.text = element_text(size=10,color = "black"),title = element_text(size=10,color = "black"),
        plot.title = element_text(size=10,color = "black",hjust = 0),
        panel.background = element_rect(fill = "transparent"),strip.text = element_text(size=16,color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background =element_blank(),
        legend.background = element_rect(fill = "transparent",color = NA),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.text=element_text(size=10))+guides(size=F,shape=F,color=F)
}

summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary4701_df_color_pti_group.rds")

plot_Ori_Psite_line("Solyc03g032150.3.1 251")
plot_imp_Psite_line("Solyc10g074710.2.1 46")
plot_Ori_Psite_box("Solyc03g032150.3.1 251")

ggsave("../SFCH/Figure/Diss_chap1/P_sites/Solyc01g088690.3.1 18.png",width =5,height =2)
ggsave("../SFCH/Figure/Diss_chap1/P_sites/Solyc03g032150.3.1 251.png",width =5,height =2.5)


dev.off()

