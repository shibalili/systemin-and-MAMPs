rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
table(summary_df$super_class)# 4701-3203=1498
summary_df$super_class <- as.character(summary_df$super_class) #remove factor
#### seperate MAPK from kinase ####
mapk <- summary_df%>%filter(str_detect(Protein_class,"MAPK"))%>%pull(Protein.ID)%>%unique()
summary_df[which(summary_df$Protein.ID %in% mapk),"super_class"] <- "MAPK"

pti <- c("Kinase","PPase","Transporter","Ca2+","ROS","MAPK","TF","hormone")

a <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")
colnames(a)[2] <- "ppep.site"
sys <- a%>%filter(treatment=="systemin",treatment_share=="1 treatment")
sys_m<- summary_df%>%filter(ppep.site%in%sys$ppep.site)%>%select(1,3:7,14,24:34)
table(sys_m$super_class)
s <- merge(sys,sys_m)%>%arrange(Sig_time_pointX)%>%relocate(super_class,super_class1,Protein_class,.before =Sig_time_pointX)%>%filter(super_class%in%pti)
table(s$super_class)
length(unique(s$Protein.ID)) #142
what_s <- plyr::count(s, vars =c("Sig_time_pointX","super_class","regulation","treatment"))%>%arrange(Sig_time_pointX,freq)
what_s$super_class <- factor(what_s$super_class,level=pti)

flg <- a%>%filter(treatment=="flg22",treatment_share=="1 treatment")
flg_m <- summary_df%>%filter(ppep.site%in%flg$ppep.site)%>%select(1,3:7,14,24:34)
f <- merge(flg,flg_m)%>%arrange(Sig_time_pointX)%>%relocate(super_class,super_class1,Protein_class,.before =Sig_time_pointX)%>%filter(super_class%in%pti)
table(f$super_class)
length(unique(f$Protein.ID)) #127
what_f <- plyr::count(f, vars =c("Sig_time_pointX","super_class","regulation","treatment"))%>%arrange(Sig_time_pointX,freq)
what_f$super_class <- factor(what_f$super_class,level=pti)

chi <- a%>%filter(treatment=="chitin",treatment_share=="1 treatment")
chi_m <- summary_df%>%filter(ppep.site%in%chi$ppep.site)%>%select(1,3:7,14,24:34)
ch <- merge(chi,chi_m)%>%arrange(Sig_time_pointX)%>%relocate(super_class,super_class1,Protein_class,.before =Sig_time_pointX)%>%filter(super_class%in%pti)
table(ch$super_class)
length(unique(ch$Protein.ID)) #99
what_c <- plyr::count(ch, vars =c("Sig_time_pointX","super_class","regulation","treatment"))%>%arrange(Sig_time_pointX,freq)
what_c$super_class <- factor(what_c$super_class,level=pti)

what <- rbind(what_s,what_f,what_c)
ggplot(what,aes(Sig_time_pointX,super_class))+geom_point(aes(color=regulation,size=freq),show.legend = T)+scale_y_discrete(limits=rev)+
  scale_x_discrete(limits =c("1","2","5","15","45") )+
  geom_text(position = "identity",aes(label=freq),colour="white",size =2.8)+
  facet_grid(regulation~treatment)+#facet_grid(~treatment)+
  scale_size_area(max_size = 8)+ #guides(size="none")+
  labs(y='',x=' ',size='P_sites number')+
  theme(axis.title=element_text(size=10,colour = 'black'), axis.text=element_text(size=10,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 0.5,fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA)) # get rid of legend panel bg
dev.off()
########  pie plot + dot plot #########
library(scatterpie)
what$freq <- as.numeric(what$freq)
pie <- what%>%pivot_wider(names_from = regulation, values_from = freq, values_fill = 0)%>% mutate(total = Up + Down )
pie$Sig_time_pointX <- as.numeric(pie$Sig_time_pointX )
pie$super_class <-  as.numeric(pie$super_class)
levels(what$super_class)# 1-8
label <- pie%>%pivot_longer(!c(1:3,6), names_to = "regulation", values_to = "freq")
#in the middle of the column (vjust = 0.5)
ggplot(data=pie) + geom_scatterpie(aes(x=Sig_time_pointX, y=super_class,r=0.025*(total)),#r=0.35
                                   cols=c("Up",'Down'), color=NA,legend_name = "Regulation") + 
  coord_equal()+facet_grid(~treatment)+
  geom_text(data = label,aes(x=Sig_time_pointX,y=super_class,group=regulation,label = freq),size = 1.5,position = position_dodge(width = 0.75), hjust = 0.5) +
  scale_x_continuous(breaks=0:5, labels = c("0","1","2","5","15","45"))+
  scale_y_continuous(trans = "reverse",breaks=1:8, labels = levels(what$super_class))+
  labs(y='',x=' ')+
  theme_minimal()+
  theme(axis.text=element_text(size=10,colour = 'black'),
        strip.text= element_text(size = 10, color = "black"),
        panel.grid.minor = element_blank())
ggsave("../SFCH/Figure/Diss_chap1/Upset_PTI.pdf",,width=8,height =5)

dev.off()


