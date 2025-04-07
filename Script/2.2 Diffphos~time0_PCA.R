rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)

## PCA input:elicitor-responsive P_sites ##
setup <- readRDS("../SFCH/Data/setup_ave.rds")
sys <- readRDS("../SFCH/Data/limma/Rem_water_systemin_sig_time0.rds")
length(unique(sys$P_site_ID)) #760
flg22<- readRDS("../SFCH/Data/limma/Rem_water_flg22_sig_time0.rds")
length(unique(flg22$P_site_ID)) #731
chitin<- readRDS("../SFCH/Data/limma/Rem_water_chitin_sig_time0.rds")  
length(unique(chitin$P_site_ID)) #616
a <- unique(c(sys$P_site_ID,flg22$P_site_ID,chitin$P_site_ID))
int <- readRDS("../SFCH/Data/Ave_int_SFCH.rds")%>%filter(rownames(int)%in%a)#1545
pca <- prcomp(x = t(na.omit(int)),center =TRUE,scale=TRUE)
summary(pca)
pred.pca <- predict(pca)
df_out <- as.data.frame(pca$x)
df_out$treatment <- setup[rownames(df_out),1]
df_out$time.point <- setup[rownames(df_out),2]
#df_out <- df_out%>%filter(treatment!="H")

ggplot(df_out,aes(x=PC1,y=PC2,color =treatment,shape=time.point,label=row.names(df_out)))+
  scale_colour_manual(name="Treatment", values= c("#D16103", "#0072B2", "#00AFBB", "#f8ff68"),
                      labels=c("Systemin","Flg22","Chitin","H2O"))+geom_point(show.legend = T,size=2)+ggtitle("PC1 and PC2")+
  scale_shape_manual(name="Time(min)",values=c(15,16,17,0,1,2))+geom_text(show.legend = F,size=3,vjust = 0, nudge_y = 0.8)+
  xlab(label = paste0("PC1:",round(summary(pca)[[6]][2,1]*100,2),"% variance"))+
  ylab(label = paste0("PC2:",round(summary(pca)[[6]][2,2]*100,2),"% variance"))+
  #scale_y_continuous(limits=c(-18, 15))+scale_x_continuous(limits=c(-21, 31))+
  guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  theme(title = element_text(size=12,color = "black"),axis.text = element_text(size=12,color = "black"),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.position="right",
        legend.background = element_rect(fill = "transparent",color = NA),legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),legend.text=element_text(size=12))
# change color palate in illustrator
#ggsave(filename = "../SFCH/Figure/Diss_chap1/SFCH_DEG_PC12_lightwater.pdf",width =5,height =4)
dev.off()

############## prep Ori:Ave 6bioreps for PCA ############## 
int <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
setup.new <- setup[setup[,1]=="S",]
s <- int[,rownames(setup.new)]#%>%select(c(1,3:7,9:13,15:18,19,21:25,27:31,33))#remove systemin rep2
S0 <- apply(s[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
S1 <- apply(s[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
S2 <- apply(s[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
S5 <- apply(s[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
S15 <-apply(s[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
S45 <-apply(s[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Ave_s<- cbind(S0,S1,S2,S5,S15,S45)
setup.new <- setup[setup[,1]=="F",]
f<- int[,rownames(setup.new)]
F0 <- apply(f[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
F1 <- apply(f[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
F2 <- apply(f[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
F5 <- apply(f[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
F15 <- apply(f[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
F45 <- apply(f[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Ave_f <- cbind(F0,F1,F2,F5,F15,F45)
setup.new <- setup[setup[,1]=="C",]
chi<- int[,rownames(setup.new)]
C0 <- apply(chi[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
C1 <- apply(chi[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
C2 <- apply(chi[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
C5 <- apply(chi[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
C15 <- apply(chi[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
C45 <- apply(chi[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Ave_chi <- cbind(C0,C1,C2,C5,C15,C45)
setup.new <- setup[setup[,1]=="H",]
water <- int[,rownames(setup.new)]
H0 <- apply(water[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
H1 <- apply(water[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
H2 <- apply(water[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
H5 <- apply(water[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
H15 <- apply(water[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
H45 <- apply(water[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Ave_water <- cbind(H0,H1,H2,H5,H15,H45)
Ave <- cbind(Ave_s,Ave_f,Ave_chi,Ave_water)%>%as.data.frame()
#saveRDS("../SFCH/Data/Ave_int_SFCH.rds",,object = Ave)

setup <- readRDS("../Rong_project2020/data/processed/setup_PCA_biorepmean.rds")
setup$treatment <- factor(setup$treatment,levels = c("S","F","C","H"))
setup <- setup[order(setup$treatment),]
setup$time.point <- factor(setup$time.point,levels = c("0","1","2","5","15","45"))
#saveRDS("../SFCH/Data/setup_ave.rds",object = setup)



