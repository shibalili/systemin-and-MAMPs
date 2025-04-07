rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
#### clustering dataframe format: col:P_site_s/f/c/H; row: 0/1/2/5/15/45 #####
###### Prep KMC: imp30_Ave 6 biorep ####
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
imp30 <-readRDS("../SFCH/Data/imputation/new_imp_PS2+MF1_threshold.rds")%>%select(c(2:133))
a <- imp30[,rownames(setup[setup[,1]=="S",])]
S0 <- apply(a[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
S1 <- apply(a[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
S2 <- apply(a[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
S5 <- apply(a[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
S15 <-apply(a[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
S45 <-apply(a[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Scount <- cbind(S0,S1,S2,S5,S15,S45)

a <- imp30[,rownames(setup[setup[,1]=="C",])]
C0 <- apply(a[,1:6], 1, function(x){return(mean(x, na.rm = TRUE))})
C1 <- apply(a[,7:12], 1, function(x){return(mean(x, na.rm = TRUE))})
C2 <- apply(a[,13:18], 1, function(x){return(mean(x, na.rm = TRUE))})
C5 <- apply(a[,19:24], 1, function(x){return(mean(x, na.rm = TRUE))})
C15 <- apply(a[,25:30], 1, function(x){return(mean(x, na.rm = TRUE))})
C45 <- apply(a[,31:36], 1, function(x){return(mean(x, na.rm = TRUE))})
Ccount <- cbind(C0,C1,C2,C5,C15,C45)

a <- imp30[,rownames(setup[setup[,1]=="F",])]
F0 <- apply(a[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
F1 <- apply(a[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
F2 <- apply(a[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
F5 <- apply(a[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
F15 <- apply(a[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
F45 <- apply(a[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Fcount <- cbind(F0,F1,F2,F5,F15,F45)

a <- imp30[,rownames(setup[setup[,1]=="H",])]
H0 <- apply(a[,1:5], 1, function(x){return(mean(x, na.rm = TRUE))})
H1 <- apply(a[,6:10], 1, function(x){return(mean(x, na.rm = TRUE))})
H2 <- apply(a[,11:15], 1, function(x){return(mean(x, na.rm = TRUE))})
H5 <- apply(a[,16:20], 1, function(x){return(mean(x, na.rm = TRUE))})
H15 <- apply(a[,21:25], 1, function(x){return(mean(x, na.rm = TRUE))})
H45 <- apply(a[,26:30], 1, function(x){return(mean(x, na.rm = TRUE))})
Hcount <- cbind(H0,H1,H2,H5,H15,H45)
allmean <- cbind(Scount,Fcount,Ccount,Hcount)%>%as.data.frame()%>%mutate(ppep.site=rownames(.),.before = 1)%>%mutate(Protein.ID=substr(rownames(.),1,18),.before = 1)
dim(allmean)#2011 26
#saveRDS("../SFCH/Data/clustering/Ave_Newimp_PS2+MF1_threshold.rds",object =allmean)
a <- readRDS("../SFCH/Data/clustering/Ave_Newimp_PS2+MF1_threshold.rds")%>%select(3:26)

## data prep for KMC ##
e <- a %>% select(contains("S"))%>%mutate(treatment="systemin")%>%mutate(ppep.site=rownames(.),.before=1)
colnames(e)[2:7] <- c("0","1","2","5","15","45")
rownames(e)<- paste0(rownames(e), "_sys")

g <- a %>% select(contains("F"))%>%mutate(treatment="flg22")%>%mutate(ppep.site=rownames(.),.before=1)
colnames(g)[2:7] <- c("0","1","2","5","15","45")
rownames(g)<- paste0(rownames(g), "_flg22")

b <- a %>% select(contains("C"))%>%mutate(treatment="chitin")%>%mutate(ppep.site=rownames(.),.before=1)
colnames(b)[2:7] <- c("0","1","2","5","15","45")
rownames(b)<- paste0(rownames(b), "_chitin")

d <- a %>% select(contains("H"))%>%mutate(treatment="H2O")%>%mutate(ppep.site=rownames(.),.before=1)
colnames(d)[2:7] <- c("0","1","2","5","15","45")
rownames(d)<- paste0(rownames(d), "_H2O")
kmc_df <- rbind(e,g,b,d)
table(kmc_df$treatment)
#saveRDS("../SFCH/Data/clustering/Newimp30_kmc_df_PS2+MF1_threshold.rds",object = kmc_df)

kmc_df <- readRDS("../SFCH/Data/clustering/Newimp30_kmc_df_PS2+MF1_threshold.rds")


