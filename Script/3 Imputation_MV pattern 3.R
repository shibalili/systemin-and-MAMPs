rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)


########## Histograms or density plot ###############
intensities <- readRDS(file = "../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(7:138)
count <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%select(6:29)
colnames(count) <- word(colnames(count),2,sep = fixed("_"))

setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
setup_ave <- readRDS("../SFCH/Data/setup_ave.rds")

#setup.new <- setup[setup[,1]=="S",]%>%filter(time.point=="1")
#setup_ave.new<- setup_ave[setup_ave[,1]=="S",]%>%filter(time.point=="1")


#count.new <- count[,rownames(setup_ave.new)]
#intensities.new <- intensities[,rownames(setup.new)]%>%mutate(ppep.site=row.names(.),count=count.new,.before = 1)#%>%filter(count>0)
#long <- intensities.new%>% pivot_longer(!c(ppep.site,count), names_to = "Sample", values_to = "detect")


long_df <- data.frame()
for (i in levels(setup$treatment)) {
  for (j in levels(setup$time.point)) {
    setup.new <- setup[setup[,1]==i,]%>%filter(time.point==j)
    setup_ave.new<- setup_ave[setup_ave[,1]==i,]%>%filter(time.point==j)
    count.new <- count[,rownames(setup_ave.new)]
    intensities.new <- intensities[,rownames(setup.new)]%>%mutate(ppep.site=row.names(.),count=count.new,.before = 1)
    long <- intensities.new%>% pivot_longer(!c(ppep.site,count), names_to = "Sample", values_to = "detect")
    long_df <- rbind(long,long_df)
  }
}

table(long_df$Sample)#4701
###### 6 biorep:systemin and chitin ##########
#sc <- long_df%>%filter(str_detect(Sample,"S|C"))
set <- setup[setup[,1]=="S",]%>%filter(time.point%in%c(0,1))
sc <- long_df%>%filter(Sample%in%rownames(set))
table(sc$Sample) #338472
table(sc$count) # 145914 NA

sc$MV <- paste0(6-sc$count,"/",6)
table(sc$MV)
### calculation #####
sc <- long_df%>%filter(str_detect(Sample,"S|C"),detect>0)
median <- median(sc$detect,na.rm = T)
sc_cal <- sc%>%filter(detect>=median,count>=3) #46163
49853/nrow(sc) # 46.3%
sc_cal <- sc%>%filter(detect>median,count<3) #4005
4005/nrow(sc) # 3.7%
sc_cal <- sc%>%filter(detect<=median,count>=3) #30875
39835/nrow(sc)# 37 %
sc_cal <- sc%>%filter(detect<median,count<3) #14023
14023/nrow(sc) #13 %



###### 5 biorep:flg22 and H2O ##########
#fh <- long_df%>%filter(str_detect(Sample,"F|H"))
fh <- long_df%>%filter(Sample%in%rownames(set))
set <- setup[setup[,1]=="F",]%>%filter(time.point%in%c(0,1))
table(fh$Sample) #282060
fh$MV <- paste0(5-fh$count,"/",5)


ggplot(sc, aes(x=detect, y=MV,colour = MV))+
  geom_point(size=2, shape=1,color="darkblue",na.rm = TRUE)+
  ggtitle("Missing Data Evaluation")+
  xlab("Normalized intensity (log2)")+ 
  ylab("Missingness%")+
  theme(axis.title= element_text(face="bold",colour="black", size=4))+
  scale_x_continuous(limits=c(18,34),breaks=seq(18, 34, by = 2))+
  geom_vline(aes(xintercept=median(detect,na.rm = TRUE)),color="red", linetype="dashed", size=1)+
  geom_segment(data = sc%>%filter(MV=="0/6"),aes(x = median(detect,na.rm = TRUE), y = "0/6", xend = median(detect,na.rm = TRUE), yend =0.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = sc%>%filter(MV=="1/6"),aes(x = median(detect,na.rm = TRUE), y = "1/6", xend = median(detect,na.rm = TRUE), yend =1.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = sc%>%filter(MV=="2/6"),aes(x = median(detect,na.rm = TRUE), y = "2/6", xend = median(detect,na.rm = TRUE), yend =2.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = sc%>%filter(MV=="3/6"),aes(x = median(detect,na.rm = TRUE), y = "3/6", xend = median(detect,na.rm = TRUE), yend =3.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = sc%>%filter(MV=="4/6"),aes(x = median(detect,na.rm = TRUE), y = "4/6", xend = median(detect,na.rm = TRUE), yend =4.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = sc%>%filter(MV=="5/6"),aes(x = median(detect,na.rm = TRUE), y = "5/6", xend = median(detect,na.rm = TRUE), yend =5.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  #geom_segment(data = sc%>%filter(MV=="6/6"),aes(x = mean(detect,na.rm = TRUE), y = "6/6", xend = mean(detect,na.rm = TRUE), yend =6.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))
ggsave("../SFCH/Figure/Diss_chap1/6biorep_NA1.pdf",width = 4,height =2.5)



ggplot(fh, aes(x=detect, y=MV,colour = MV))+
  geom_point(size=2, shape=1,color="darkblue",na.rm = TRUE)+
  ggtitle("Missing Data Evaluation")+
  xlab("Normalized intensity (log2)")+ 
  ylab("Missingness%")+
  theme(axis.title= element_text(face="bold",colour="black", size=4))+
  scale_x_continuous(limits=c(18,34),breaks=seq(18, 34, by = 2))+
  geom_vline(aes(xintercept=median(detect,na.rm = TRUE)),color="red", linetype="dashed", size=1)+
  geom_segment(data = fh%>%filter(MV=="0/5"),aes(x = median(detect,na.rm = TRUE), y = "0/5", xend = median(detect,na.rm = TRUE), yend =0.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = fh%>%filter(MV=="1/5"),aes(x = median(detect,na.rm = TRUE), y = "1/5", xend = median(detect,na.rm = TRUE), yend =1.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = fh%>%filter(MV=="2/5"),aes(x = median(detect,na.rm = TRUE), y = "2/5", xend = median(detect,na.rm = TRUE), yend =2.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = fh%>%filter(MV=="3/5"),aes(x = median(detect,na.rm = TRUE), y = "3/5", xend = median(detect,na.rm = TRUE), yend =3.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  geom_segment(data = fh%>%filter(MV=="4/5"),aes(x = median(detect,na.rm = TRUE), y = "4/5", xend = median(detect,na.rm = TRUE), yend =4.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  #geom_segment(data = fh%>%filter(MV=="5/5"),aes(x = median(detect,na.rm = TRUE), y = "5/5", xend = median(detect,na.rm = TRUE), yend =5.5),arrow = arrow(ends ="first",length = unit(0.2, "cm")))+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("../SFCH/Figure/Diss_chap1/5biorep_NA1.pdf",width =4,height = 2.5)
dev.off()






####### useful tool for NA detection #########
library(mice)
library(visdat)
library(naniar)
library(VIM)
library(gridExtra)

intensities <- readRDS(file = "../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(7:138)
#glimpse(intensities)
# 'aggr' plots the amount of missing/imputed values in each column(sample-wise)
#aggr(intensities[1:6,1:3])
# missing pattern visualization
# Visualises a data.frame to tell you what it contains.
ave_int <- readRDS(file = "../SFCH/Data/Ave_int_SFCH.rds")
#vis_dat(ave_int[1:500,]) # check the missing of 500 P_sites (average 6 biological reps)
#vis_miss(ave_int,sort_miss = F)#+coord_flip()
#gg_miss_case_cumsum(ave_int)#cumulative sum of missing values, reading the rows of the dataset from the top to bottom
#md.pattern(ave_int[1:36,] ,rotate.names = TRUE)
# explore missingness relationships
#ggplot(intensities, aes(x = Intensity.S01, y = Intensity.S02)) + geom_miss_point()
# samplewise-row
Tint_df <- as.data.frame(t(intensities))
# Two convenient counters of complete values and missings are n_miss() and n_complete(). 
n_miss(Tint_df)/(n_miss(Tint_df)+n_complete(Tint_df))
n_complete(Tint_df)/(n_miss(Tint_df)+n_complete(Tint_df))
#prop_miss_case(Tint_df)
#pct_miss_case(Tint_df)

# ppepwise/variablewise-columne:missing across samples
b <- miss_var_table(Tint_df)









