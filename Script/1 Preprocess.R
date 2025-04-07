rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
dat <- read.table("../SFCH/Data/Phospho (STY)Sites.txt",header=TRUE,sep = "\t",quote="",fill = FALSE) #6533 1369
dat <- dat%>%mutate(ppep.site=paste(substr(dat[,4],1,18),dat[,'Positions.within.proteins']) ,.before = 1)# add col:ppep.site
dat <- dat[-grep("+",dat[,"Reverse"], fixed=TRUE),]
dat <- dat[-grep("+",dat[,"Potential.contaminant"], fixed=TRUE),]
dat <- dat[dat[,"Localization.prob"]>0.75,] #4701 P_site 
rownames(dat) <- paste(substr(dat[,5],1,18),dat[,'Positions.within.proteins']) 
#identical(rownames(dat),dat$ppep.site)
dim(dat) #4701 1370
grep("Intensity.C01",colnames(dat)) #689
grep("Intensity.S56",colnames(dat)) #820
int <- dat[,689:820]
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
int <- int[,rownames(setup)] # order treatment (S-F-C-H) & time point (0-1-2-5-15-45)  
int[int==0] <- NA
normalize.intensities <- function(x){
  x0 <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
  overall.mean <- mean(apply(x, 2, function(x){return(sum(x, na.rm = TRUE))}))#average the sum intensity of each column/sample
  print(overall.mean)
  for(i in 1:ncol(x0)){x0[,i] <- overall.mean*x[,i]/sum(x[,i], na.rm = TRUE)}
  colnames(x0) <- colnames(x)
  rownames(x0) <- rownames(x)
  return(x0)
}
intensities<- normalize.intensities(int)%>%as.data.frame()%>%log2() #73917901300
seq <- dat%>%select(c(542:544,547))%>%mutate(Flanking = gsub("\\;.*$", "", Sequence.window)%>%str_sub(9,23))
colnames(seq)[1] <- "No.STY"
colnames(seq)[4] <- "Probability"
all.equal(rownames(intensities),rownames(seq))
Summary4701 <-cbind(seq,intensities)%>%mutate(ppep.site=row.names(.),.before = 1)
#saveRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds",object =Summary4701)



