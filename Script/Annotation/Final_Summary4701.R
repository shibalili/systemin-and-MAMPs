rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)

summary_df <- readRDS("../SFCH/Data/annotation/summary4701/SFCH_summary4701_PSMF_threshold_Mfuzz_df.rds") #185 col
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/SFCH_summary4701_PSMF_threshold_Mfuzz_protein_class_df.rds") #190 col
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary4701_df_color.rds") #191 col
summary_df <- readRDS("../SFCH/Data/annotation/summary4701/summary4701_df_color_pti_group.rds")#191 col


int_norm <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")
########### number of measurement ########### 
#identical(rownames(int_norm),rownames(count))
count <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")
df <- bind_cols(count,int_norm)%>%mutate(Protein.ID=substr(ppep.site,1,18),.before = ppep.site)%>%relocate(Protein.ID,ppep.site)# 168 col

########### annotation ########### 
mapping <- readRDS("../SFCH/Data/annotation/phytozome_ITAG3.2_standard.rds")
ss <- mapping[unique(df$Protein.ID),c(1:10)]
df <- merge(ss,df,by="Protein.ID")%>%relocate(Protein.ID,ppep.site)
rownames(df) <- df$ppep.site
df <- df[rownames(int_norm),] #order rownames as the same order of int_norm

########### DEG should include significant P_sites in Mock(water) ########### 
DEG<- readRDS("../SFCH/Data/limma/Rem_water_SFC and H_only_a.rds")%>%mutate(combine = paste(treatment,regulation,Sig_time_pointX))
table(DEG$treatment,DEG$treatment_share)
Order_DEG <- DEG[order(DEG$Sig_time_pointX),] 
Order_DEG <- DEG[order(DEG$treatment,decreasing = F),] # systemin comes first
df <- df%>%mutate(DEG="no",.before = 1)
l <- list()
for (i in 1:nrow(df)) {
  l[[i]]<- Order_DEG[which(Order_DEG$P_site_ID==rownames(df)[i]),"combine"]
  names(l)[i] <- rownames(df)[i]}#identical(names(l),rownames(df))
df$DEG<- l%>%as.character()
df[df$DEG=="character(0)","DEG"] <- "no"
#the regex "^c\\(|\\)$" only affects the first/last chars in the string, remove "c(" and ")" from string
df$DEG <- df$DEG%>%gsub("\"", "", .)%>%gsub("^c\\(|\\)$", "",.) 

########### clustering ########### 
df <- df%>%mutate(H2O_cluster="filtered",.before = 1)
df <- df%>%mutate(chitin_cluster="filtered",.before = 1)
df <- df%>%mutate(flg22_cluster="filtered",.before = 1)
df <- df%>%mutate(sys_cluster="filtered",.before = 1)
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")
table(distribution$sys)
df[which(rownames(df) %in% rownames(distribution)),]$sys_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"sys"]
table(df$sys_cluster)

table(distribution$flg22)
df[which(rownames(df) %in% rownames(distribution)),]$flg22_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"flg22"]
table(df$flg22_cluster)

table(distribution$chitin)
df[which(rownames(df) %in% rownames(distribution)),]$chitin_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"chitin"]
table(df$chitin_cluster)

table(distribution$H2O)
df[which(rownames(df) %in% rownames(distribution)),]$H2O_cluster<- distribution[rownames(df[which(rownames(df) %in% rownames(distribution)),]),"H2O"]
table(df$H2O_cluster)

########### annotation: X4.2_symbol ########### 
X4.2_symbol <- readRDS("../SFCH/Data/annotation/phytozome_ITAG3.2_X4.2.rds")
df <- df%>%mutate(X4.2_symbol=X4.2_symbol[df$Protein.ID,c("X4.2_symbol")])%>%relocate(X4.2_symbol,.before = arabi.defline)
########### annotation: Localization ########### 
df$Ara_SUBA3 <- gsub("[\"]","",df$Ara_SUBA3)%>%noquote()
df$Localization <- "null"
df[which(str_detect(df$Ara_SUBA3, "plasma")),"Localization"] <- "PM"
df[which(str_detect(df$Ara_SUBA3, "cytosol")),"Localization"] <- "cytosol"
df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"]=="null","ER",df[which(str_detect(df$Ara_SUBA3, "endoplasmic")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"]=="null","nucleus",df[which(str_detect(df$Ara_SUBA3, "nucleus")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"]=="null","extracellular",df[which(str_detect(df$Ara_SUBA3, "extracellular")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"]=="null","plastid",df[which(str_detect(df$Ara_SUBA3, "plastid")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"]=="null","mitochondrion",df[which(str_detect(df$Ara_SUBA3, "mitochondrion")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"]=="null","golgi",df[which(str_detect(df$Ara_SUBA3, "golgi")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"]=="null","peroxisome",df[which(str_detect(df$Ara_SUBA3, "peroxisome")),"Localization"])
df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"] <- ifelse(df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"]=="null","vacuole",df[which(str_detect(df$Ara_SUBA3, "vacuole")),"Localization"])
sort(table(df$Localization))
df<- df%>%relocate(Localization,.before = Ara_SUBA3)

########### annotation:Protein.type (rough) ########### 
mapping <- readRDS("../SFCH/Data/annotation/phytozome_ITAG3.2_standard.rds")
a <- mapping[unique(df$Protein.ID),]
TF<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "transcription factor")))%>%
  unique()%>%pull(Protein.ID)
receptor <- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "receptor")))%>%
  unique()%>%pull(Protein.ID)
kinase<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "kinase")))%>%
  unique()%>%pull(Protein.ID)
fake_kinase <- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "kinase")))%>%
  filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "phosphatase")))
##### more true kinase than PPase so PPase first ####
phostase<- a%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern = "phosphatase")))%>%
  unique()%>%pull(Protein.ID)
df$Protein.type <- "unknown"
df[which(df$Protein.ID %in% TF),"Protein.type"] <- "TF"
df[which(df$Protein.ID %in% receptor),"Protein.type"] <- "receptor"
df[which(df$Protein.ID %in% phostase),"Protein.type"] <- "PPase"
df[which(df$Protein.ID %in% kinase),"Protein.type"] <- "kinase"
df[which(df$Protein.ID%in%c("Solyc03g118350.3.1","Solyc10g005640.3.1")),"Protein.type"] <- "PPase"
table(df$Protein.type)
df<- df%>%relocate(Protein.type,.after = Protein.ID)
colnames(df)

df1 <- df[,c(8,6:7,49:53,5,1:4,20:24,9:19,25:48,54:185)]
#saveRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds",object=df1) # 185 col

########### add post-defined elicitor specificity feature and shorten bin names ########### 
## shorten bin names
summary_df <-readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")%>%mutate(Mapman=paste("Bin",word(MapmanX4.2_Bincode,sep = fixed("."))),.before =MapmanX4.2_Bincode )
table(summary_df$Mapman)
summary_df[which(summary_df$Mapman=="Bin 0"),"Mapman"] <- 0
table(summary_df$Mapman)
## elicitor specificity feature
summary_df <-summary_df%>%mutate(specify=0,.before = ppep.site)
library(ComplexHeatmap)
a <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")
sys <- a%>%filter(treatment=="systemin")
flg22 <- a%>%filter(treatment=="flg22")
chitin <- a%>%filter(treatment=="chitin")
list <- list(systemin=unique(sys$P_site_ID),flg22=unique(flg22$P_site_ID),chitin=unique(chitin$P_site_ID))
allcomb <- make_comb_mat(list,mode = "distinct")#all possible intersections/combinations
#c("100","010","001","110","101","011","111")
#c("sys only","flg22 only","chitin only","sys&flg22","sys&chitin","flg22&chitin","sys&flg22&chitin"))
table(summary_df$specify)
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"100")),"specify"] <- "S"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"010")),"specify"] <- "F"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"001")),"specify"] <- "C"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"110")),"specify"] <- "SF"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"101")),"specify"] <- "SC"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"011")),"specify"] <- "FC"
summary_df[which(summary_df$ppep.site%in%extract_comb(allcomb,"111")),"specify"] <- "SFC"
summary_df$specify <- factor(summary_df$specify,levels = c("S","F","C","SF","SC","FC","SFC","0"))
table(summary_df$specify)
#saveRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds",object=summary_df) # 187 col

########### annotation:protein_class (specify pti components) ########### 
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")# 187 col
mapping <- readRDS("../SFCH/Data/annotation/mapping_protein_class.rds")%>%mutate(super_class=word(Protein_class,sep = fixed(":")),.before = Protein_class)%>%select(1:4)
ss <- merge(mapping,summary_df,by ="Protein.ID",no.dups =T )
rownames(ss) <- ss$ppep.site
df <- ss[rownames(summary_df),] # 190 col
table(df$super_class)
#saveRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds",object=df) # 190 col


summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")%>%mutate(super_class1=0,.after = super_class)# 191 col
sort(table(summary_df$super_class),decreasing = T)
length(unique(summary_df$super_class)) #28=27+1
summary_df$super_class[summary_df$super_class%in%c("transport","Mechano")] <- "Transporter"
summary_df$super_class[summary_df$super_class%in%c("PPase","PPase/BR")] <- "PPase"
summary_df$super_class[summary_df$super_class%in%c("RLK","RLCK","RLP","PK","MAPK")] <- "Kinase"
summary_df$super_class[summary_df$super_class%in%c("ABA","ET","TF/ET","Auxin","BR","GA","JA")] <- "hormone"
sort(table(summary_df$super_class),decreasing = T)
length(unique(summary_df$super_class)) #16=15+1
complete <- c("Kinase","PPase","Ca2+","ROS","TF","hormone","Transporter","cellwall","PR","ETI","NLR","cyclin","MLO","TPR","14-3-3","0")
summary_df$super_class <- factor(summary_df$super_class,levels = complete)

############ arrange pti components #################
## PTI = Ca2+, hormone, MAPK, NLRs, PRs, PPase,RLPs, RLKs, RLCKs, and TFs.
complete <- c("Kinase","PPase","Ca2+","ROS","TF","hormone","Transporter","cellwall","PR","ETI","NLR","cyclin","MLO","TPR","14-3-3","0")
pti <- c("Kinase","PPase","Transporter","Ca2+","ROS","TF","hormone")
pti_early <- c("Kinase","PPase","Ca2+","ROS")
#df$Protein.type <- factor(df$Protein.type,levels=c("Kinase","PPase","Ca2+","ROS","ETI","NLR","TF","hormone","cyclin","MLO","TPR","Transporter","cellwall"))
#filter(!super_class%in%c("0","cyclin","TPR","ETI","PK","NLR"))
#filter(!super_class%in%c("0","TPR","ETI","MLO","cyclin","Transporter","TF","Kinase","Ca2+","hormone","cellwall","ASP","NLR","TF:WRKY","PPase","TF/ET","14-3-3"))

###### add one more col:super_class1 
dd <- summary_df%>%filter(super_class=="Ca2+")
table(dd$Protein_class)
summary_df[which(str_detect(summary_df$Protein_class,"ACA")),"super_class1"] <- "Ca2+:ACA"
summary_df[which(str_detect(summary_df$Protein_class,"CIPK")),"super_class1"] <- "Ca2+:CIPK"
summary_df[which(str_detect(summary_df$Protein_class,"CML")),"super_class1"] <- "Ca2+:CML"
summary_df[which(str_detect(summary_df$Protein_class,"CNGC")),"super_class1"] <- "Ca2+:CNGC"
summary_df[which(str_detect(summary_df$Protein_class,"CPK")),"super_class1"] <- "Ca2+:CPK"
summary_df[which(str_detect(summary_df$Protein_class,"IQD")),"super_class1"] <- "Ca2+:IQD"
summary_df[which(str_detect(summary_df$Protein_class,"OSCA")),"super_class1"] <- "Ca2+:OSCA"
table(summary_df$super_class1)

dd <- summary_df%>%filter(super_class=="hormone")
summary_df[which(str_detect(summary_df$Protein_class,"SnRK")),"super_class1"] <- "ABA:SnRK"
summary_df[which(str_detect(summary_df$Protein_class,"PIN")),"super_class1"] <- "Auxin:PIN"
summary_df[which(str_detect(summary_df$Protein_class,"BR:SlSK|BSK|GSK")),"super_class1"] <- "BR:SK"
summary_df[which(str_detect(summary_df$Protein_class,"CTR")),"super_class1"] <- "ET:CTR"
summary_df[which(str_detect(summary_df$Protein_class,"ETR")),"super_class1"] <- "ET:ETR"

dd <- summary_df%>%filter(super_class=="Kinase")
summary_df[which(str_detect(summary_df$Protein_class,"MAPKKK")),"super_class1"] <- "MAP3K"
summary_df[which(str_detect(summary_df$Protein_class,"MAP4K")),"super_class1"] <- "MAP4K"
summary_df[which(summary_df$Protein_class%in%c("MAPK:SlMAPK2","MAPK:SlMAPK4","MAPK:SlMAPK5")),"super_class1"] <- "MAPK"
summary_df[which(summary_df$Protein_class%in%c("MAPK:SlMAPKK1")),"super_class1"] <- "MAP2K"
summary_df[which(str_detect(summary_df$Protein_class,"PBL|BIK|PBS")),"super_class1"] <- "RLCK:PBL"
summary_df[which(str_detect(summary_df$Protein_class,"LRK10")),"super_class1"] <- "RLK:WAK"
summary_df[which(str_detect(summary_df$Protein_class,"CrRLK")),"super_class1"] <- "RLK:CrRLK1L"
summary_df[which(str_detect(summary_df$Protein_class,"LecRK")),"super_class1"] <- "RLK:LecRK"
summary_df[which(str_detect(summary_df$Protein_class,"LYK")),"super_class1"] <- "RLK:LYK"

dd <- summary_df%>%filter(super_class=="PPase")
summary_df[which(str_detect(summary_df$Protein_class,"PP2C")),"super_class1"] <- "PP2C"
#summary_df[which(str_detect(summary_df$Protein_class,"PLL")),"super_class1"] <- "PP2C:PLL"
summary_df[which(str_detect(summary_df$Protein_class,"BSL")),"super_class1"] <- "BR:BSL"
summary_df[which(str_detect(summary_df$Protein_class,"PTP")),"super_class1"] <- "PTP"
summary_df[which(str_detect(summary_df$Protein_class,"PP2A")),"super_class1"] <- "PPP:PP2A"
summary_df[which(str_detect(summary_df$Protein_class,"PP4")),"super_class1"] <- "PPP:PP4"
summary_df[which(str_detect(summary_df$Protein_class,"PP6")),"super_class1"] <- "PPP:PP6"
summary_df[which(str_detect(summary_df$Protein_class,"ASP")),"super_class1"] <- "ASP"

dd <- summary_df%>%filter(super_class=="ETI")
summary_df[which(str_detect(summary_df$Protein_class,"RIN4")),"super_class1"] <- "ETI:RIN4"

dd <- summary_df%>%filter(super_class=="ROS")
summary_df[which(str_detect(summary_df$Protein_class,"Rboh")),"super_class1"] <- "ROS:Rboh"

dd <- summary_df%>%filter(str_detect(super_class,"TF"))
summary_df[which(str_detect(summary_df$Protein_class,"CAMTA")),"super_class1"] <- "TF:CAMTA"
summary_df[which(str_detect(summary_df$Protein_class,"WRKY")),"super_class1"] <- "TF:WRKY"
summary_df[which(str_detect(summary_df$Protein_class,"MYB")),"super_class1"] <- "TF:MYB"

dd <- summary_df%>%filter(str_detect(super_class,"Transporter"))
summary_df[which(str_detect(summary_df$Protein_class,"Mechano")),"super_class1"] <- "Transport:Mechano"
summary_df[which(str_detect(summary_df$Protein_class,"ABC")),"super_class1"] <- "Transport:ABC"
summary_df[which(str_detect(summary_df$Protein_class,"PIP")),"super_class1"] <- "Transport:PIP"
summary_df[which(str_detect(summary_df$Protein_class,"LHA")),"super_class1"] <- "H+-ATPase"
summary_df[which(str_detect(summary_df$Protein_class,"CLC")),"super_class1"] <- "Transport:CLC"
summary_df[which(str_detect(summary_df$Protein_class,"SOS")),"super_class1"] <- "Transport:SOS"
summary_df[which(str_detect(summary_df$Protein_class,"SLAC")),"super_class1"] <- "Transport:SLAC"
summary_df[which(str_detect(summary_df$Protein_class,"HAK")),"super_class1"] <- "Transport:HAK"

#saveRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds",object=summary_df) # 191 col
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")

####### for network version ###########
####### add WGCNA ###########
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")%>%mutate(colors="no",.before = super_class)
# color
module_df <- readRDS("../SFCH/Data/WGCNA/netwk13/netwk13_module_df.rds")
table(module_df$colors)
summary_df$colors <- module_df[summary_df$ppep.site,"colors"]
table(summary_df$colors)
# module hub
hub <- readRDS("../SFCH/Data/WGCNA/netwk13/hub.rds") #303
summary_df <- summary_df%>%mutate(WGCNA_hub="N",.after=colors)
table(summary_df$WGCNA_hub)
summary_df[hub$ppep.site,"WGCNA_hub"] <- "Y"
# GRN hub
grn_hubs <- readRDS("../SFCH/Data/WGCNA/genie3/grn_hubs.rds") #102
summary_df <- summary_df%>%mutate(GRN_hub="N",.after=WGCNA_hub)
table(summary_df$GRN_hub)
summary_df[grn_hubs$ppep.site,"GRN_hub"] <- "Y"

#saveRDS("../SFCH/Data/annotation/SFCH_summary4701_df_Net.rds",object=summary_df) # 194 col

sumNet <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df_Net.rds")




