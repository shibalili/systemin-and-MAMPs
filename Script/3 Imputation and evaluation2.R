rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(naniar)#as_shadow


# Ratio of MVs versus average signal intensity??

#df <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%select(-c(1:5))%>%mutate(ppep.site=row.names(.),.before=1)
#df_long <- df%>%pivot_longer(!ppep.site,names_to = "sample",values_to = "miss_ratio")
#df_long$miss_ratio <- as.numeric(df_long$miss_ratio/6)
#df_long$sample <- word(df_long$sample,2,sep = fixed("_"))
int <-readRDS("../SFCH/Data/Ave_int_SFCH.rds")%>%mutate(ppep.site=row.names(.),.before=1)
int_long <- int%>%pivot_longer(!ppep.site,names_to = "sample",values_to = "Intensity")%>%mutate(as_shadow(.))%>%select(c(1:3,6))
table(int_long$Intensity_NA)
#int_long$miss_ratio <- df_long$miss_ratio
#ddf <- plyr::count(int_long, vars =c("Intensity","Intensity_NA"))


### all detected P_sites intensity distribution ######  
ggplot(int_long,aes(x = Intensity))+geom_histogram(aes(y = ..density..),colour = 1, fill = "transparent")+
  geom_density(adjust=1,kernel = "epanechnikov",na.rm = TRUE,show.legend = FALSE)+scale_y_continuous(limits=c(0, 1))+
  theme(strip.text = element_text(size = 14))+theme_classic()
ggsave("../SFCH/Figure/Diss_chap1/imp_allsites_dist.pdf",width = 3,height = 2)
dev.off()

