rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
####### batch effect:samplewise ##########
library(limma)#batch effect
library(ggplot2)
library(reshape2)
library(stringr)
library(dplyr)

int <- readRDS("../SFCH/Data/Intensity_4701P_sites_Norm.rds")%>%select(c(7:138))
setup <- readRDS("../SFCH/Data/setup_allorder_SL3.0.rds")
###### df prep######
rep_sum <- NULL
for (i in 1:6) {
  setup.new <- setup[setup[,1]=="S"&setup[,3]==i,]# change treatment accordingly
  s <- int[,rownames(setup.new)]
  rep_sum<-cbind(rep_sum,apply(s, 1, function(x){return(sum(x, na.rm = TRUE))}))
}
colnames(rep_sum) <- c("systemin_rep1","systemin_rep2","systemin_rep3","systemin_rep4","systemin_rep5","systemin_rep6")
sys_sum_rep6 <- rep_sum
colnames(rep_sum) <- c("flg22_rep1","flg22_rep2","flg22_rep3","flg22_rep4","flg22_rep5","flg22_rep6")
flg22_sum_rep5 <- rep_sum[,-2]
colnames(rep_sum) <- c("chitin_rep1","chitin_rep2","chitin_rep3","chitin_rep4","chitin_rep5","chitin_rep6")
chitin_sum_rep6 <- rep_sum
colnames(rep_sum) <- c("H2O_rep1","H2O_rep2","H2O_rep3","H2O_rep4","H2O_rep5","H2O_rep6")
H2O_sum_rep5 <- rep_sum[,-6]
df <- as.data.frame(cbind(sys_sum_rep6,flg22_sum_rep5,chitin_sum_rep6,H2O_sum_rep5))
#saveRDS("../SFCH/Data/Aggreg_4701P_22col.rds",object = df)

df <- readRDS("../SFCH/Data/Aggreg_4701P_22col.rds")
df1 <- na.omit(t(scale(t(df))))%>%as.data.frame()
df <- df[rownames(df1),]
boxplot(df)
setup_batch <- readRDS("../SL3.0_run5/Data_SL3.0/processed/Preprocessed data/Aggregation_Setup_Batch.rds")
setup_batch$batch <- factor(setup_batch$batch)
rm_batch <- as.data.frame(removeBatchEffect(df[,1:22], batch = setup_batch$batch, design = model.matrix( ~ treatment, data = setup_batch)))
boxplot(rm_batch)
library(factoextra)
df1 <- t(na.omit(t(scale(t(rm_batch)))))
hc.res <- eclust(df1, "hclust", k =4, hc_metric = "pearson", hc_method = "complete",graph = FALSE,seed=123)
fviz_dend(hc.res,horiz = TRUE, k_colors = c("#FC4E07","#00AFBB","#2E9FDF","#E7B800"),color_labels_by_k = TRUE, # color labels by groups
          rect_fill = T,rect=T,lwd = 0.5,cex=0.8, main = "Dendrogram of Sample")
#ggsave("../SFCH/Figure/Diss_chap1/hc_rmbatch.pdf",width =5,height=5)
dev.off()
