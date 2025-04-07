rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(naniar)#as_shadow
### cutoff count>30 ###
df <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%filter(Count_Sum>30) # 2011
int <-readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%filter(rownames(.)%in%rownames(df))%>%select(c(1,7:138))
int_long <- int%>%pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(as_shadow(.))%>%select(c(1:3,6))

n_missing <- sum(is.na(int_long$Intensity)) # 105445
sd <- sd(int_long$Intensity,na.rm = T) #σ=1.831334 Standard Deviation;μ=23.97223
shrinked_sd <- 0.3 * sd # shrink sd width 0.549
downshifted_mean <- mean(int_long$Intensity,na.rm = T) - 1.8 * sd # shift mean of imputed values 20.67583
set.seed(123)
int_long[is.na(int_long$Intensity),"Intensity"] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
int_long$Intensity_NA <- factor(int_long$Intensity_NA,labels=c("Original data","imputed data"))
#saveRDS("../SFCH/Data/imputation/SFCH_imp30_int_long_Perseus.rds",object=int_long)

## Perseus-type ##
int_wide <- int_long%>%select(1:3)%>%pivot_wider(names_from = sample, values_from = Intensity)%>%
  mutate(Protein.ID=substr(ppep.site,1,18),.before = 1)%>%as.data.frame()
rownames(int_wide) <- int_wide$ppep.site
#saveRDS("../SFCH/Data/imputation/SFCH_imp30_Perseus.rds",object=int_wide)

Perseus <- readRDS("../SFCH/Data/imputation/SFCH_imp30_Perseus.rds")


## MissForest ##
library(missForest)#MAR/MCAR imp
set.seed(123)#in order to get reproducible results for people to repeat or check!!
int <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
int <- t(int[apply(int, 1, function(x){sum(!is.na(x))>30}),]) #2011
int.imp <- missForest(xmis = int,verbose = TRUE)
str(int.imp)
#saveRDS("../SFCH/Data/imputation/SFCH_int.imp.rds",object=int.imp)
imp30 <- t(int.imp$ximp)%>%as.data.frame()%>%mutate(ppep.site=rownames(.),.before = 1)%>%mutate(Protein.ID=substr(rownames(.),1,18),.before = 1)
#saveRDS("../SFCH/Result/SFCH_imp30.rds",object=imp30)


#### replace imputation of imp30(by MissForest) with Perseus ####
setup<- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")%>%mutate(name=paste0(treatment,time.point))
imp30 <- readRDS("../SFCH/Data/imputation/SFCH_imp30.rds")%>%select(3:134)#imputation among bio-reps
Perseus <- readRDS("../SFCH/Data/imputation/SFCH_imp30_Perseus.rds")%>%select(3:134)
i=NULL
j=NULL
# create an exmpty data.frame
new_imp <- data.frame(matrix(ncol =ncol(imp30), nrow = nrow(imp30)),row.names = rownames(imp30))
colnames(new_imp) <- colnames(imp30)
dim(new_imp)
################# Perseus2 (PS2) : if one time point has 0 or 1 measure(i.e.<2) ####################
for (i in 1:nrow(imp30)) {
  for (j in 1:ncol(df)) {
    if (df[which(rownames(imp30) == rownames(Perseus)[i]),j]<2) {
      aa <- rownames(setup)[setup$name==word(colnames(df)[j],2,sep = fixed("_"))]
      new_imp[which(rownames(imp30) == rownames(Perseus)[i]),aa] <- Perseus[rownames(Perseus)[i],aa]}
    else if (df[which(rownames(imp30) == rownames(Perseus)[i]),j]>=2){
      aa <- rownames(setup)[setup$name==word(colnames(df)[j],2,sep = fixed("_"))]
      new_imp[which(rownames(imp30) == rownames(Perseus)[i]),aa] <- imp30[rownames(Perseus)[i],aa]}
  }
}
#saveRDS("../SFCH/Data/imputation/new_imp_Perseus>2.rds",object=new_imp)

################# Missforest1 (MF1): if the only one measure has extremely high intensity >threshold ################# 
i=NULL
j=NULL
new_imp <-readRDS("../SFCH/Data/imputation/new_imp_Perseus>2.rds") 
NA_distribution <- readRDS("../SFCH/Data/imputation/SFCH_imp30_int_long_Perseus.rds")%>%filter(Intensity_NA=="imputed data")
threshold <-qnorm(p = 0.975,lower.tail = T,mean=mean(NA_distribution$Intensity),sd=sd(NA_distribution$Intensity,na.rm = T))#2.5%
for (i in 1:nrow(imp30)) {
  for (j in 1:ncol(df)) {
    if (df[which(rownames(imp30) == rownames(Perseus)[i]),j]==1) {
      aa <- rownames(setup)[setup$name==word(colnames(df)[j],2,sep = fixed("_"))]
      ifelse (int[rownames(imp30)[i],aa[!is.na(int[rownames(imp30)[i],aa])]]>threshold,## threshold of NA distribution
              new_imp[which(rownames(imp30) == rownames(Perseus)[i]),aa] <- imp30[rownames(Perseus)[i],aa],
              new_imp[which(rownames(imp30) == rownames(Perseus)[i]),aa] <- new_imp[rownames(Perseus)[i],aa])}
  }
}
new_imp <- new_imp%>%mutate(ppep.site=rownames(.),.before=1)
#saveRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds",object=new_imp)

new_imp <- readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")#%>%select(c(2:133))

##### imputation_evaluation_PS2MF1_threshold #######
new_imp_long <-readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")%>%pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")
new_imp_long$Intensity_NA <- int_long$Intensity_NA
new_imp_long$Intensity_NA <- factor(new_imp_long$Intensity_NA,labels=c("Original data","MissForest"))
levels(new_imp_long$Intensity_NA) <- c(levels(new_imp_long$Intensity_NA) ,"Perseus")

imp30_long <- readRDS("../SFCH/Data/imputation/SFCH_imp30.rds")%>%select(c(2:134))%>%pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")
## takes 12 min ..
job::job({
  for (i in 1:nrow(new_imp_long)) {
    ifelse(new_imp_long[i,"Intensity"]==imp30_long[i,"Intensity"],
           new_imp_long$Intensity_NA[i] <- new_imp_long$Intensity_NA[i],
           new_imp_long$Intensity_NA[i] <- "Perseus")}}) 
#saveRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold_evaluation.rds",object=new_imp_long)
table(new_imp_long$Intensity_NA)
new_imp_long <- readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold_evaluation.rds")

ggplot(new_imp_long,aes(x = Intensity,colour = Intensity_NA,fill =Intensity_NA))+geom_histogram(alpha=0.8,bins = 45)+
  scale_fill_grey(start = 0.1, end = 0.9)+scale_color_grey(start = 0.2, end = 0.8)+
  theme(strip.text = element_text(size = 14))+theme_classic()+facet_grid(.~Intensity_NA)
#+scale_fill_brewer(palette = "Dark2")+scale_color_brewer(palette = "Dark2")#

#ggsave("../SFCH/Figure/Diss_chap1/New_imp_Eva_histo.pdf",width =6,height =2) 

dev.off()

######### plot mean and SD before and after imputation ##########
SD_ori <- readRDS("../SFCH/Data/imputation/Ori_SD.rds")
SD_imp <- readRDS("../SFCH/Data/imputation/newimp_SD.rds")
SD_df <- rbind(SD_ori,SD_imp)
SD_df$type <- factor(SD_df$type,levels = c("Ori","Imp"))
SD_df$treatment <- factor(SD_df$treatment,levels = c("S","F","C","H"))
SD_df$attr <- "SD"

Mean_ori <- readRDS("../SFCH/Data/imputation/Ori_mean.rds")
Mean_imp <- readRDS("../SFCH/Data/imputation/newimp_mean.rds")
Mean_df <- rbind(Mean_ori,Mean_imp)
Mean_df$type <- factor(Mean_df$type,levels = c("Ori","Imp"))
Mean_df$treatment <- factor(Mean_df$treatment,levels = c("S","F","C","H"))
Mean_df$attr <- "Mean"

eva_df <- rbind(Mean_df,SD_df)
eva_df %>%as.data.frame() %>% ggplot(aes(x=type, y=Intensity, fill=type))+geom_boxplot(outlier.size =0.5,na.rm = TRUE,outlier.colour = "grey")+
  facet_grid(attr~treatment,scales = "free_y")+scale_fill_grey(start = 0.1, end = 0.9)+scale_color_grey(start = 0.2, end = 0.8)+theme_bw()
#ggsave("../SFCH/Figure/Diss_chap1/New_imp_Eva_SD&Mean_box.pdf",width =6,height =3) 


######## prep SD ########
df <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%filter(Count_Sum>30) # 2011
#Ori
#m <-readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%filter(rownames(.)%in%rownames(df))%>%select(c(1,7:138))
#new imp
m <- readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")#%>%select(c(2:133))
setup<- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")#%>%mutate(name=paste0(treatment,time.point))

### SD ###
a <- m[,rownames(setup[setup[,1]=="S",])]
S0 <- apply(a[,1:6], 1, function(x){sd(x,na.rm = TRUE)})
S1 <- apply(a[,7:12], 1, function(x){sd(x,na.rm = TRUE)})
S2 <- apply(a[,13:18], 1, function(x){sd(x,na.rm = TRUE)})
S5 <- apply(a[,19:24], 1, function(x){sd(x,na.rm = TRUE)})
S15 <- apply(a[,25:30], 1, function(x){sd(x,na.rm = TRUE)})
S45 <- apply(a[,31:36], 1, function(x){sd(x,na.rm = TRUE)})
Scount <- t(rbind(S0,S1,S2,S5,S15,S45))

a <- m[,rownames(setup[setup[,1]=="C",])]
C0 <- apply(a[,1:6], 1, function(x){return(sd(x, na.rm = TRUE))})
C1 <- apply(a[,7:12], 1, function(x){return(sd(x, na.rm = TRUE))})
C2 <- apply(a[,13:18], 1, function(x){return(sd(x, na.rm = TRUE))})
C5 <- apply(a[,19:24], 1, function(x){return(sd(x, na.rm = TRUE))})
C15 <- apply(a[,25:30], 1, function(x){return(sd(x, na.rm = TRUE))})
C45 <- apply(a[,31:36], 1, function(x){return(sd(x, na.rm = TRUE))})
Ccount <- cbind(C0,C1,C2,C5,C15,C45)

a <- m[,rownames(setup[setup[,1]=="F",])]
F0 <- apply(a[,1:5], 1, function(x){return(sd(x, na.rm = TRUE))})
F1 <- apply(a[,6:10], 1, function(x){return(sd(x, na.rm = TRUE))})
F2 <- apply(a[,11:15], 1, function(x){return(sd(x, na.rm = TRUE))})
F5 <- apply(a[,16:20], 1, function(x){return(sd(x, na.rm = TRUE))})
F15 <- apply(a[,21:25], 1, function(x){return(sd(x, na.rm = TRUE))})
F45 <- apply(a[,26:30], 1, function(x){return(sd(x, na.rm = TRUE))})
Fcount <- cbind(F0,F1,F2,F5,F15,F45)

a <- m[,rownames(setup[setup[,1]=="H",])]
H0 <- apply(a[,1:5], 1, function(x){return(sd(x, na.rm = TRUE))})
H1 <- apply(a[,6:10], 1, function(x){return(sd(x, na.rm = TRUE))})
H2 <- apply(a[,11:15], 1, function(x){return(sd(x, na.rm = TRUE))})
H5 <- apply(a[,16:20], 1, function(x){return(sd(x, na.rm = TRUE))})
H15 <- apply(a[,21:25], 1, function(x){return(sd(x, na.rm = TRUE))})
H45 <- apply(a[,26:30], 1, function(x){return(sd(x, na.rm = TRUE))})
Hcount <- cbind(H0,H1,H2,H5,H15,H45)
#allsd <- cbind(Scount,Fcount,Ccount,Hcount)%>%as.data.frame()%>%mutate(ppep.site=row.names(.),.before = 1)%>%
#pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(type="Ori",treatment=substr(sample, 1, 1))
#saveRDS("../SFCH/Data/imputation/Ori_SD.rds",object =allsd )
allsd <- cbind(Scount,Fcount,Ccount,Hcount)%>%as.data.frame()%>%mutate(ppep.site=row.names(.),.before = 1)%>%
  pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(type="Imp",treatment=substr(sample, 1, 1))
#saveRDS("../SFCH/Data/imputation/newimp_SD.rds",object=allsd )

###  Mean ###

a <- m[,rownames(setup[setup[,1]=="S",])]
S0 <- apply(a[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
S1 <- apply(a[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
S2 <- apply(a[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
S5 <- apply(a[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
S15 <-apply(a[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
S45 <-apply(a[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Scount <- cbind(S0,S1,S2,S5,S15,S45)

a <- m[,rownames(setup[setup[,1]=="C",])]
C0 <- apply(a[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
C1 <- apply(a[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
C2 <- apply(a[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
C5 <- apply(a[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
C15 <- apply(a[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
C45 <- apply(a[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Ccount <- cbind(C0,C1,C2,C5,C15,C45)

a <- m[,rownames(setup[setup[,1]=="F",])]
F0 <- apply(a[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
F1 <- apply(a[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
F2 <- apply(a[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
F5 <- apply(a[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
F15 <- apply(a[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
F45 <- apply(a[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Fcount <- cbind(F0,F1,F2,F5,F15,F45)

a <- m[,rownames(setup[setup[,1]=="H",])]
H0 <- apply(a[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
H1 <- apply(a[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
H2 <- apply(a[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
H5 <- apply(a[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
H15 <- apply(a[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
H45 <- apply(a[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Hcount <- cbind(H0,H1,H2,H5,H15,H45)
#allmean <- cbind(Scount,Fcount,Ccount,Hcount)%>%as.data.frame()%>%mutate(ppep.site=row.names(.),.before = 1)%>%
#pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(type="Ori",treatment=substr(sample, 1, 1))
#saveRDS("../SFCH/Data/imputation/Ori_mean.rds",object =allmean )

allmean <- cbind(Scount,Fcount,Ccount,Hcount)%>%as.data.frame()%>%mutate(ppep.site=row.names(.),.before = 1)%>%
  pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(type="Imp",treatment=substr(sample, 1, 1))
#saveRDS("../SFCH/Data/imputation/newimp_mean.rds",object=allmean )

