rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)
library(limma)
int <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
setup.new <- setup[setup[,1]=="S",] #S F C H
int <- int[,rownames(setup.new)]
levels(setup.new$time.point)# 0 as reference: Compare to time0
design <- model.matrix( ~ time.point, data = setup.new)
colSums(design)
fit <- lmFit(int, design)#fit a lm
fit <- eBayes(fit)#t-test and so on
#summary by topTable
res1 <- topTable(fit, coef = "time.point1", number = Inf,adjust.method = "BH")
res2 <- topTable(fit, coef = "time.point2", number = Inf,adjust.method = "BH")
res5 <- topTable(fit, coef = "time.point5", number = Inf,adjust.method = "BH")
res15 <- topTable(fit, coef ="time.point15", number = Inf,adjust.method = "BH")
res45 <- topTable(fit, coef ="time.point45", number = Inf,adjust.method = "BH")
test_df <- rbind(data.frame(res1,P_site_ID=rownames(res1),Sig_time_pointX="1"),
                 data.frame(res2,P_site_ID=rownames(res2),Sig_time_pointX="2"),
                 data.frame(res5,P_site_ID=rownames(res5),Sig_time_pointX="5"),
                 data.frame(res15,P_site_ID=rownames(res15),Sig_time_pointX="15"),
                 data.frame(res45,P_site_ID=rownames(res45),Sig_time_pointX="45"))
test_df$Sig_time_pointX <- factor(test_df$Sig_time_pointX, levels =c("1","2","5","15","45"))
test_df$treatment <- "systemin"
head(test_df)
rownames(test_df) <- NULL
#saveRDS("../SFCH/Data/limma/sys_limmatest_df.rds",object =test_df )
#saveRDS("../SFCH/Data/limma/flg22_limmatest_df.rds",object =test_df )
#saveRDS("../SFCH/Data/limma/chitin_limmatest_df.rds",object =test_df )
#saveRDS("../SFCH/Data/limma/H2O_limmatest_df.rds",object =test_df )

time1 <- na.omit(res1[res1[,"P.Value"]<0.05,])
time1$P_site_ID <- rownames(time1)
time1$Sig_time_pointX <- "1"
time1$regulation <- ifelse(time1$logFC>0,"Up","Down")# P and deP Trend determine by LogFC

time2 <- na.omit(res2[res2[,"P.Value"]<0.05,])
time2$P_site_ID <- rownames(time2)
time2$Sig_time_pointX <- "2"
time2$regulation <- ifelse(time2$logFC>0,"Up","Down")

time5 <- na.omit(res5[res5[,"P.Value"]<0.05,])
time5$P_site_ID <- rownames(time5)
time5$Sig_time_pointX <- "5"
time5$regulation <- ifelse(time5$logFC>0,"Up","Down")

time15 <- na.omit(res15[res15[,"P.Value"]<0.05,])
time15$P_site_ID <- rownames(time15)
time15$Sig_time_pointX <- "15"
time15$regulation <- ifelse(time15$logFC>0,"Up","Down")

time45 <- na.omit(res45[res45[,"P.Value"]<0.05,])
time45$P_site_ID <- rownames(time45)
time45$Sig_time_pointX <- "45"
time45$regulation <- ifelse(time45$logFC>0,"Up","Down")

head(time1)
df <-rbind(time1[,c(4:5,7:9)],time2[,c(4:5,7:9)],time5[,c(4:5,7:9)],time15[,c(4:5,7:9)],time45[,c(4:5,7:9)])  
df$Sig_time_pointX <- factor(df$Sig_time_pointX, levels =c("1","2","5","15","45"))
table(df$Sig_time_pointX)
df$treatment <- "systemin"
head(df)
rownames(df) <- NULL
df <- df[,c(3:6,1:2)]
table(df$regulation,df$Sig_time_pointX)
#saveRDS("../SFCH/Data/limma/sys_pvalue.rds",object = df)
#saveRDS("../SFCH/Data/limma/flg22_pvalue.rds",object = df)
#saveRDS("../SFCH/Data/limma/chitin_pvalue.rds",object = df)
#saveRDS("../SFCH/Data/limma/H2O_pvalue.rds",object = df)

sys <- readRDS("../SFCH/Data/limma/sys_pvalue.rds")
flg22<- readRDS("../SFCH/Data/limma/flg22_pvalue.rds")
chitin<- readRDS("../SFCH/Data/limma/chitin_pvalue.rds")
H2O <- readRDS("../SFCH/Data/limma/H2O_pvalue.rds")
sys$combine <- paste(sys$P_site_ID,sys$regulation,sys$Sig_time_pointX)
flg22$combine <- paste(flg22$P_site_ID,flg22$regulation,flg22$Sig_time_pointX)
chitin$combine <- paste(chitin$P_site_ID,chitin$regulation,chitin$Sig_time_pointX)
H2O$combine <- paste(H2O$P_site_ID,H2O$regulation,H2O$Sig_time_pointX)

nrow(sys%>%filter(combine%in%intersect(sys$combine, H2O$combine)))#-61
nrow(flg22%>%filter(combine%in%intersect(flg22$combine, H2O$combine)))#-47
nrow(chitin%>%filter(combine%in%intersect(chitin$combine, H2O$combine)))#-51
overlap <- c(sys%>%filter(combine%in%intersect(sys$combine, H2O$combine))%>%pull(P_site_ID),
             flg22%>%filter(combine%in%intersect(flg22$combine, H2O$combine))%>%pull(P_site_ID),
             chitin%>%filter(combine%in%intersect(chitin$combine, H2O$combine))%>%pull(P_site_ID))%>%unique()
H2O_only <- H2O%>%filter(!P_site_ID%in%overlap)
#saveRDS("../SFCH/Data/limma/H2O_only_sig_time0.rds",object = H2O_only)

p <- NULL
pp <- NULL
ppp <- NULL
a <- NULL
### now it's sys loop ####
for (j in 1:length(levels(sys$Sig_time_pointX))) {
  p <- sys%>%filter(Sig_time_pointX==levels(sys$Sig_time_pointX)[j])
  pp <-H2O$combine
  a <- p[p$combine %in% pp,]
  ppp <- rbind(ppp,p[!p$combine%in%a$combine,])
}
nrow(ppp)
#saveRDS("../SFCH/Data/limma/Rem_water_systemin_sig_time0.rds",object = ppp) #1106
#saveRDS("../SFCH/Data/limma/Rem_water_flg22_sig_time0.rds",object = ppp) #1025
#saveRDS("../SFCH/Data/limma/Rem_water_chitin_sig_time0.rds",object = ppp) #898

####### plot :significantly-phosphorylated sites over time (upon treatment after time 0) ####### 
# Remove Deg_P_sites with the same trend in water                   
sys <- readRDS("../SFCH/Data/limma/Rem_water_systemin_sig_time0.rds")
flg22<- readRDS("../SFCH/Data/limma/Rem_water_flg22_sig_time0.rds")
chitin<- readRDS("../SFCH/Data/limma/Rem_water_chitin_sig_time0.rds")    
#H2O <- readRDS("../SFCH/Data/limma/H2O_only_sig_time0.rds")
#Deg_time0 <- rbind(sys[,1:4],flg22[,1:4],chitin[,1:4],H2O[,1:4])
Deg_time0 <- rbind(sys[,1:4],flg22[,1:4],chitin[,1:4])
table(Deg_time0$treatment,Deg_time0$regulation)
treatment_share <- plyr::count(Deg_time0, vars =c("Sig_time_pointX","P_site_ID","regulation"))
a <- merge(treatment_share,Deg_time0)
colnames(a)[4] <- "treatment_share"
table(a$treatment)
table(a$regulation)
table(a$treatment_share)
table(a$Sig_time_pointX)
a$treatment <- factor(a$treatment, levels =c("systemin","flg22","chitin")) #,"H2O"
a$regulation <- factor(a$regulation, levels =c("Up","Down"))
a$treatment_share <- factor(a$treatment_share, levels =c("1","2","3"),labels = c("1 treatment","2 treatments","3 treatments"))
#saveRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds",object = a)
#saveRDS("../SFCH/Data/limma/Rem_water_SFC and H_only_a.rds",object = a)


##### SFCH: if water is included in the plot ######
sys <- readRDS("../SFCH/Result/sys_pvalue.rds")
flg22<- readRDS("../SFCH/Result/flg22_pvalue.rds")
chitin<- readRDS("../SFCH/Result/chitin_pvalue.rds")
H2O <- readRDS("../SFCH/Result/H2O_pvalue.rds")
Deg_time0 <- rbind(sys[,1:4],flg22[,1:4],chitin[,1:4],H2O[,1:4])
table(Deg_time0$treatment,Deg_time0$regulation)
treatment_share <- plyr::count(Deg_time0, vars =c("Sig_time_pointX","P_site_ID","regulation"))
a <- merge(treatment_share,Deg_time0)
colnames(a)[4] <- "treatment_share"
a$treatment <- factor(a$treatment, levels =c("systemin","flg22","chitin","H2O"))
a$regulation <- factor(a$regulation, levels =c("Up","Down"))
a$treatment_share <- factor(a$treatment_share, levels =c("1","2","3","4"),labels = c("1 treatment","2 treatments","3 treatments","4 treatments"))
#saveRDS("../SFCH/Result/SFCH_a.rds",object=a)

#### now make the plot:DEG_P_site number summarized according to their sig time (x-axis) ####
a <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")
#a <- readRDS("../SFCH/Data/limma/SFCH_a.rds")
count <- plyr::count(a, vars =c("Sig_time_pointX","treatment","regulation","treatment_share"))
count$freq <- ifelse(count$regulation=="Down",-1*count$freq,count$freq)
ggplot(data=count, aes(x=Sig_time_pointX, y=freq,fill=treatment,alpha =treatment_share))+
  geom_bar(stat="identity",position="stack",linewidth=0.3,color="lightgrey")+
  scale_alpha_discrete(range = c(0.2,1))+#scale_fill_brewer(palette = "Set1", name = "treatment",guide = "none")+
  scale_fill_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"),labels=c("Systemin","Flg22","Chitin","H2O"),guide = "none")+
  facet_grid(~as.factor(treatment))+theme_bw(base_size = 10)+
  labs(y= "Count differentially phosphrylated", x="Time (min)", alpha="# Patterns")+
  theme(axis.text=element_text(size=12,colour = "black"), axis.title=element_text(size=14),
        axis.line = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.background = element_blank(), plot.background = element_blank(), 
        legend.background = element_blank(), legend.box.background = element_blank(),
        panel.grid.major =element_blank() ,#element_line(colour = "white",size = 0.01)
        panel.grid.minor = element_blank(),panel.border = element_blank(), 
        legend.box=NULL, legend.key.size=unit(0.6, "cm"), legend.position ="top", 
        legend.text=element_text(size=14), legend.title = element_text(size=14), 
        strip.background = element_blank(),strip.text.x = element_text(size=14)) 
# no water
#ggsave("../SFCH/Figure/Diss_chap1/DEG_Psites_Up_and_Down.pdf", height=9, width=15, units="cm", bg="transparent")

dev.off()





