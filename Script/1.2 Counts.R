rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
count <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%mutate(ppep.site=rownames(.),.before = 1)%>%
          mutate(Protein.ID=substr(rownames(.),1,18),.before = 1)%>%filter(Count_Sum>0) #4567

########## prep Counts_SFCH4701.rds ################
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
intensities <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
setup.new <- setup[setup[,1]=="S",]
intensities.new <- intensities[,rownames(setup.new)]
S0 <- apply(intensities.new[,1:6], 1, function(x){(length(x)-sum(is.na(x)))})
S1 <- apply(intensities.new[,7:12], 1, function(x){(length(x)-sum(is.na(x)))})
S2 <- apply(intensities.new[,13:18], 1, function(x){(length(x)-sum(is.na(x)))})
S5 <- apply(intensities.new[,19:24], 1, function(x){(length(x)-sum(is.na(x)))})
S15 <- apply(intensities.new[,25:30], 1, function(x){(length(x)-sum(is.na(x)))})
S45 <- apply(intensities.new[,31:36], 1, function(x){(length(x)-sum(is.na(x)))})
Scount <- rbind(S0,S1,S2,S5,S15,S45)
Scount <- t(Scount)

setup.new <- setup[setup[,1]=="C",]
intensities.new <- intensities[,rownames(setup.new)]
C0 <- apply(intensities.new[,1:6], 1, function(x){(length(x)-sum(is.na(x)))})
C1 <- apply(intensities.new[,7:12], 1, function(x){(length(x)-sum(is.na(x)))})
C2 <- apply(intensities.new[,13:18], 1, function(x){(length(x)-sum(is.na(x)))})
C5 <- apply(intensities.new[,19:24], 1, function(x){(length(x)-sum(is.na(x)))})
C15 <- apply(intensities.new[,25:30], 1, function(x){(length(x)-sum(is.na(x)))})
C45 <- apply(intensities.new[,31:36], 1, function(x){(length(x)-sum(is.na(x)))})
Ccount <- rbind(C0,C1,C2,C5,C15,C45)
Ccount <- t(Ccount)

setup.new <- setup[setup[,1]=="F",]
intensities.new <- intensities[,rownames(setup.new)]
F0 <- apply(intensities.new[,1:5], 1, function(x){(length(x)-sum(is.na(x)))})
F1 <- apply(intensities.new[,6:10], 1, function(x){(length(x)-sum(is.na(x)))})
F2 <- apply(intensities.new[,11:15], 1, function(x){(length(x)-sum(is.na(x)))})
F5 <- apply(intensities.new[,16:20], 1, function(x){(length(x)-sum(is.na(x)))})
F15 <- apply(intensities.new[,21:25], 1, function(x){(length(x)-sum(is.na(x)))})
F45 <- apply(intensities.new[,26:30], 1, function(x){(length(x)-sum(is.na(x)))})
fcount <- rbind(F0,F1,F2,F5,F15,F45)
fcount <- t(fcount)

setup.new <- setup[setup[,1]=="H",]
intensities.new <- intensities[,rownames(setup.new)]
H0 <- apply(intensities.new[,1:5], 1, function(x){(length(x)-sum(is.na(x)))})
H1 <- apply(intensities.new[,6:10], 1, function(x){(length(x)-sum(is.na(x)))})
H2 <- apply(intensities.new[,11:15], 1, function(x){(length(x)-sum(is.na(x)))})
H5 <- apply(intensities.new[,16:20], 1, function(x){(length(x)-sum(is.na(x)))})
H15 <- apply(intensities.new[,21:25], 1, function(x){(length(x)-sum(is.na(x)))})
H45 <- apply(intensities.new[,26:30], 1, function(x){(length(x)-sum(is.na(x)))})
Hcount <- rbind(H0,H1,H2,H5,H15,H45)
Hcount <- t(Hcount)
allcount <- cbind(Scount,fcount,Ccount,Hcount)%>%as.data.frame()
alltreatmentcount <- data.frame(systemin=apply(allcount[,1:6], 1, function(x){sum(x)}),
                                flg22=apply(allcount[,7:12], 1, function(x){sum(x)}),
                                chitin=apply(allcount[,13:18], 1, function(x){sum(x)}),
                                H2O=apply(allcount[,19:24], 1, function(x){sum(x)}))
alltreatmentcount <- cbind(alltreatmentcount,allcount)%>%
  mutate(Sum=apply(alltreatmentcount[,1:4], 1,function(x){return(sum(x, na.rm = TRUE))}))%>%relocate(Sum)
colnames(alltreatmentcount) <- paste0("Count_",colnames(alltreatmentcount))
#saveRDS("../SFCH/Data/Counts_SFCH4701.rds",object = alltreatmentcount)

################# hypergeometeric test on counts ##############################
#### which treatment has the significantly higher phosphorylation level?/significantly de/phosphorylated in at least one condition
measure <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")
systemin=36
hyperg <- NULL
for(i in 1:nrow(measure)){
  m <- systemin
  n <- 132-m
  k <- measure[i,"Count_Sum"]
  x <- measure[i,"Count_systemin"]
  hyperg <- rbind(hyperg,phyper(x-1, m, n, k,lower.tail = FALSE))##accumulated p value=probability + probability of extreme case
}
colnames(hyperg) <- "Count_sys_p.value"
rownames(hyperg) <- rownames(measure)
hyperg <- cbind(hyperg,Count_sys_adj.p.val = p.adjust(hyperg[,"Count_sys_p.value"],method = "BH"))
sum(hyperg[,"Count_sys_adj.p.val"]<0.05)#290
#systemin specific due to meansurements on sys
sp_sys <- data.frame(hyperg,measure)%>%mutate(ppep.site=rownames(.),.before = 1)#%>%arrange(Count_sys_adj.p.val)

flg22=30
hyperg <- NULL
for(i in 1:nrow(measure)){
  m <- flg22
  n <- 132-m
  k <- measure[i,"Count_Sum"]
  x <- measure[i,"Count_flg22"]
  hyperg <- rbind(hyperg,phyper(x-1, m, n, k,lower.tail = FALSE))##accumulated p value=probability + probability of extreme case
}
rownames(hyperg) <- rownames(measure)
colnames(hyperg) <- "Count_flg22_p.value"
hyperg <- cbind(hyperg,Count_flg22_adj.p.val = p.adjust(hyperg[,"Count_flg22_p.value"],method = "BH"))
sum(hyperg[,"Count_flg22_adj.p.val"]<0.05)#70
sp_flg22 <-data.frame(hyperg,measure)%>%mutate(ppep.site=rownames(.),.before = 1)#%>%arrange(Count_flg22_adj.p.val)
#identical(rownames(sp_sys),rownames(sp_flg22))
chitin=36
hyperg <- NULL
for(i in 1:nrow(measure)){
  m <- chitin
  n <- 132-m
  k <- measure[i,"Count_Sum"]
  x <- measure[i,"Count_chitin"]
  hyperg <- rbind(hyperg,phyper(x-1, m, n, k,lower.tail = FALSE))##accumulated p value=probability + probability of extreme case
}
rownames(hyperg) <- rownames(measure)
colnames(hyperg) <- "Count_chi_p.value"
hyperg <- cbind(hyperg,Count_chi_adj.p.val = p.adjust(hyperg[,"Count_chi_p.value"],method = "BH"))
sum(hyperg[,"Count_chi_adj.p.val"]<0.05) #420
sp_chi <-data.frame(hyperg,measure)%>%mutate(ppep.site=rownames(.),.before = 1)#%>%arrange(Count_chi_adj.p.val)
#identical(rownames(sp_sys),rownames(sp_chi))
H2O=30
hyperg <- NULL
for(i in 1:nrow(measure)){
  m <- H2O
  n <- 132-m
  k <- measure[i,"Count_Sum"]
  x <- measure[i,"Count_H2O"]
  hyperg <- rbind(hyperg,phyper(x-1, m, n, k,lower.tail = FALSE))##accumulated p value=probability + probability of extreme case
}
rownames(hyperg) <- rownames(measure)
colnames(hyperg) <- "Count_H2O_p.value"
hyperg <- cbind(hyperg,Count_H2O_adj.p.val = p.adjust(hyperg[,"Count_H2O_p.value"],method = "BH"))
sum(hyperg[,"Count_H2O_adj.p.val"]<0.05) #0 !
sp_water <- data.frame(hyperg,measure)%>%mutate(ppep.site=rownames(.),.before = 1)#%>%arrange(Count_H2O_adj.p.val)
library(SuperExactTest)
#identical(rownames(sp_sys),rownames(sp_water))
sp <- cbind(sp_sys[,2:3],sp_flg22[,2:3],sp_chi[,2:3],sp_water[,2:3],measure)
sp <- sp%>%mutate(Count_specific_at=NA,.before = 1)
sp$Count_specific_at[which(sp$Count_sys_adj.p.val<0.05)] <- "systemin" 
table(sp$Count_specific_at) #290
sp$Count_specific_at[which(sp$Count_flg22_adj.p.val<0.05)] <- "flg22" 
table(sp$Count_specific_at) #70
sp$Count_specific_at[which(sp$Count_chi_adj.p.val<0.05)] <- "chitin" 
table(sp$Count_specific_at) #420
#sp$Count_specific_at[which(sp$Count_H2O_adj.p.val<0.05)]
sp <- sp%>%mutate(ppep.site=rownames(.),.before = 1) # no overlap among SFC
#saveRDS("../SFCH/Data/Count_specific_at_SFCH4701.rds",object=sp)

sp <- readRDS("../SFCH/Data/Count_specific_at_SFCH4701.rds")

################### Specificity measure (SPM): for the heatmap visualization ######################
##### SPM_Sum: firstly sum 6 reps then sum 6 timepoints #######
int <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
setup.new <- setup[setup[,1]=="S",]
s <- int[,rownames(setup.new)]
S0 <- apply(s[,1:6], 1, function(x){return(sum(x, na.rm = TRUE))})
S1 <- apply(s[,7:12], 1, function(x){return(sum(x, na.rm = TRUE))})
S2 <- apply(s[,13:18], 1, function(x){return(sum(x, na.rm = TRUE))})
S5 <- apply(s[,19:24], 1, function(x){return(sum(x, na.rm = TRUE))})
S15 <-apply(s[,25:30], 1, function(x){return(sum(x, na.rm = TRUE))})
S45 <-apply(s[,31:36], 1, function(x){return(sum(x, na.rm = TRUE))})
sum<- cbind(S0,S1,S2,S5,S15,S45)
Sum_s<- apply(sum, 1, function(x){return(sum(x, na.rm = TRUE))})
setup.new <- setup[setup[,1]=="F",]
f<- int[,rownames(setup.new)]
F0 <- apply(f[,1:5], 1, function(x){return(sum(x, na.rm = TRUE))})
F1 <- apply(f[,6:10], 1, function(x){return(sum(x, na.rm = TRUE))})
F2 <- apply(f[,11:15], 1, function(x){return(sum(x, na.rm = TRUE))})
F5 <- apply(f[,16:20], 1, function(x){return(sum(x, na.rm = TRUE))})
F15 <- apply(f[,21:25], 1, function(x){return(sum(x, na.rm = TRUE))})
F45 <- apply(f[,26:30], 1, function(x){return(sum(x, na.rm = TRUE))})
sum <- cbind(F0,F1,F2,F5,F15,F45)
Sum_f<- apply(sum, 1, function(x){return(sum(x, na.rm = TRUE))})
setup.new <- setup[setup[,1]=="C",]
chi<- int[,rownames(setup.new)]
C0 <- apply(chi[,1:6], 1, function(x){return(sum(x, na.rm = TRUE))})
C1 <- apply(chi[,7:12], 1, function(x){return(sum(x, na.rm = TRUE))})
C2 <- apply(chi[,13:18], 1, function(x){return(sum(x, na.rm = TRUE))})
C5 <- apply(chi[,19:24], 1, function(x){return(sum(x, na.rm = TRUE))})
C15 <- apply(chi[,25:30], 1, function(x){return(sum(x, na.rm = TRUE))})
C45 <- apply(chi[,31:36], 1, function(x){return(sum(x, na.rm = TRUE))})
sum <- cbind(C0,C1,C2,C5,C15,C45)
Sum_chi<- apply(sum, 1, function(x){return(sum(x, na.rm = TRUE))})
setup.new <- setup[setup[,1]=="H",]
water <- int[,rownames(setup.new)]
H0 <- apply(water[,1:5], 1, function(x){return(sum(x, na.rm = TRUE))})
H1 <- apply(water[,6:10], 1, function(x){return(sum(x, na.rm = TRUE))})
H2 <- apply(water[,11:15], 1, function(x){return(sum(x, na.rm = TRUE))})
H5 <- apply(water[,16:20], 1, function(x){return(sum(x, na.rm = TRUE))})
H15 <- apply(water[,21:25], 1, function(x){return(sum(x, na.rm = TRUE))})
H45 <- apply(water[,26:30], 1, function(x){return(sum(x, na.rm = TRUE))})
sum <- cbind(H0,H1,H2,H5,H15,H45)
Sum_water <- apply(sum, 1, function(x){return(sum(x, na.rm = TRUE))})
SPM_Sum<- as.data.frame(cbind(Sum_s,Sum_f,Sum_chi,Sum_water))
colnames(SPM_Sum) <- c("systemin","flg22","chitin","H2O")
SPM_Sum$Total <- apply(SPM_Sum, 1,function(x){return(sum(x, na.rm = TRUE))})#SPM_Sum <- SPM_Sum%>%filter(Total>0)#4567
#saveRDS("../SFCH/Data/Sum_int_SFCH.rds",object=SPM_Sum)

SPM_Sum <- readRDS("../SFCH/Data/Sum_int_SFCH.rds")
SPM_fraction <- NULL
for (i in 1:nrow(SPM_Sum)) {
  frac <- SPM_Sum[i,1:4]/SPM_Sum[i,5]
  SPM_fraction <- rbind(SPM_fraction,frac)
}
SPM_fraction <- na.omit(SPM_fraction)
cutoff=0 #how much difference among treatments-horizontal difference
# If P_sites with at least one SPM>0.33 (one treatment responsible for approximately 1/3 total phosphorylation intensity) were gathered
sys_spm <- na.omit(rownames(SPM_fraction)[SPM_fraction$systemin>cutoff])
length(sys_spm)
flg22_spm <- na.omit(rownames(SPM_fraction)[SPM_fraction$flg22>cutoff]) 
length(flg22_spm)
chitin_spm <-na.omit(rownames(SPM_fraction)[SPM_fraction$chitin>cutoff])
length(chitin_spm)
water_spm <- na.omit(rownames(SPM_fraction)[SPM_fraction$H2O>cutoff])
length(water_spm)
spm_df <- data.frame(ID=unique(c(sys_spm,flg22_spm,chitin_spm,water_spm)),
                     SPM_fraction[unique(c(sys_spm,flg22_spm,chitin_spm,water_spm)),],
                     SPM_treatment=colnames(SPM_fraction)[max.col(SPM_fraction[unique(c(sys_spm,flg22_spm,chitin_spm,water_spm)),])])
spm_df$SPM_treatment <- factor(spm_df$SPM_treatment, levels=c("systemin","flg22","chitin","H2O"))
table(spm_df$SPM_treatment)
### Order P_sites by max.SPM 
spm_df1 <- spm_df%>%filter(SPM_treatment=="systemin")
spm_df2 <- spm_df%>%filter(SPM_treatment=="flg22")
spm_df3 <- spm_df%>%filter(SPM_treatment=="chitin")
spm_df4 <- spm_df%>%filter(SPM_treatment=="H2O")
spm_df_plot <- rbind(spm_df1[order(spm_df1$systemin,decreasing = T),],
                     spm_df2[order(spm_df2$flg22,decreasing = T),],
                     spm_df3[order(spm_df3$chitin,decreasing = T),],
                     spm_df4[order(spm_df4$H2O,decreasing = T),])
table(spm_df_plot$SPM_treatment)
# count specificity by hypergeometeric test
sp <- readRDS("../SFCH/Data/Count_specific_at_SFCH4701.rds")%>%filter(!is.na(Count_specific_at))#780
spm_df_plot <- spm_df_plot%>%filter(ID%in%sp$ppep.site)
table(spm_df_plot$SPM_treatment)

############# wide to long for ggplot: the heatmap !! ##########################
library(reshape2)
long <- melt(spm_df_plot,id.vars = c("ID","SPM_treatment"),variable.name = c("treatment"), value.name = "SPM")
ggplot(long, aes(x=treatment, y=factor(ID, levels=rev(unique(ID))), fill=SPM))+ 
  theme_minimal(base_size=7) + geom_tile()+labs(x="", y="", fill="Specificity measure (SPM)")+
  scale_fill_viridis_c(guide=guide_colorbar(barwidth=unit(2.5, "cm"),barheight=unit(0.25, "cm")))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_text(color = "black",size = 8),rect=element_rect(fill="transparent",color=NA),
        panel.grid.major=element_blank(),legend.position="top",plot.margin=unit(c(0.1,0.4,0.1,0.4), "cm"),
        legend.box.spacing = unit(0.2, "cm"),legend.title = element_text(vjust = 1,color = "black",size =  8),
        legend.text =element_text(color = "black",size =8))
#ggsave("../SFCH/Figure/SPM_cutoff0.png", height=10, width=8, units='cm')
dev.off()

