rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ComplexHeatmap)
library(hypeR)
a <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")
sys <- a%>%filter(treatment=="systemin")
flg22 <- a%>%filter(treatment=="flg22")
chitin <- a%>%filter(treatment=="chitin")
list <- list(systemin=unique(sys$P_site_ID),flg22=unique(flg22$P_site_ID),chitin=unique(chitin$P_site_ID))
allcomb <- make_comb_mat(list,mode = "distinct")#all possible intersections/combinations
#extract_comb(allcomb,"100")
set_size(allcomb)
comb_size <- data.frame(code=names(comb_size(allcomb)),gene.size=comb_size(allcomb))
#### background list ####
back.ID <- readRDS("../SFCH/Data/Counts_SFCH4701.rds")%>%mutate(ppep.site=rownames(.),.before = 1)%>%
  mutate(Protein.ID=substr(rownames(.),1,18),.before = 1)%>%pull(Protein.ID)%>%unique() #2119 proteins

gl <- list()
for(x in names(comb_size(allcomb))){gl[[x]]<- unique(substr(extract_comb(allcomb,x),1,18))}
names(gl) <- names(comb_size(allcomb))
lengths(gl)

anno <- readRDS("../SFCH/Data/annotation/SFCH_bin1.rds")
#anno <- readRDS("../SFCH/Data/annotation/SFCH_bin1.1.rds")
IDs.sum <- apply(anno[,-c(1:4)], 2, sum)
str(IDs.sum)
relevant.bin <-names(IDs.sum)[IDs.sum>0]#how many bins(functional groups) in background
gl_bin <- list()
for(k in names(comb_size(allcomb))){
  gene.list <- gl[[k]]
  gene.set <- list() 
  anno.gene.list <- list() 
  for(i in 1:length(relevant.bin)){
    bin <- relevant.bin[i]
    gene.set[[i]] <- rownames(anno)[anno[,bin]==1]#as long as background is set, universe geneset is fixed
    anno.gene.list[[i]] <- intersect(gene.set[[i]], gene.list)}# anno genelist is varied due to diff clusters
  names(gene.set) <- relevant.bin #combine protein ID (from background list) with bincode
  names(anno.gene.list) <- relevant.bin
  gl_bin [[k]]<- names(anno.gene.list[lapply(anno.gene.list, length) > 0 ])}
names(gl_bin) <- names(gl)

b <- list() #map gene ID from background to binname
annotation <- readRDS("../SL3.0_run5/Data_SL3.0/Annotation source/ITAG3.2 mapman_Ensambl44.rds")
for(j in 1:length(gl_bin)){
  b[[j]] <- gene.set[gl_bin[[j]]]
  bin.name <- annotation[,2:3] %>% filter(BINCODE %in% gl_bin[[j]]) %>%unique()
  names(b[[j]]) <-paste0("bin ",bin.name[,1]," ",bin.name[,2])}
names(b) <- names(gl)

hyp_obj <- list()
for (m in names(comb_size(allcomb))) {hyp_obj[[m]] <- hypeR(signature=gl[[m]], genesets=b[[m]],test = "hypergeometric",background=length(back.ID))}
names(hyp_obj) <- names(gl) 
names(hyp_obj)

hyp_to_excel <- NULL
for (z in names(comb_size(allcomb))) {
  df <- data.frame(hyp_obj[[z]]$data,cluster=z)
  hyp_to_excel <- rbind(hyp_to_excel,df)}
rownames(hyp_to_excel) <- NULL

y <- NULL
for (n in 1:length(names(comb_size(allcomb)))) {
  s <- strsplit(hyp_obj[[n]]$data[["hits"]], split = ",")
  y <-rbind(y,data.frame(sig_cluster=rep(names(hyp_obj)[n], sum(lengths(s))),
                         bin.name= rep(hyp_obj[[n]]$data[["label"]], sapply(s, length)),
                         pvalue= rep(hyp_obj[[n]]$data[["pval"]], sapply(s, length)),
                         Protein.ID = unlist(s)%>%gsub("\\s+","", .))) }# remove whitespace in a string
## map protein ID to ppep ID
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_PSMF_threshold_Mfuzz_df.rds")
ppep.y1 <- NULL
for (r in names(comb_size(allcomb))) {
  site_kn <- summary_df%>%filter(ppep.site %in% extract_comb(allcomb,r))
  y_kn <- y%>%filter(sig_cluster==r)
  ppep.y <- subset(site_kn,site_kn$Protein.ID %in% y_kn$Protein.ID) 
  ppep.y1 <- rbind(ppep.y1,merge(y_kn,ppep.y,by ="Protein.ID"))%>%arrange(sig_cluster,pvalue)}

ppep.y1$sig_cluster <- factor(ppep.y1$sig_cluster,
                              levels =c("100","010","001","110","101","011","111"),
                              labels = c("sys only","flg22 only","chitin only","sys&flg22","sys&chitin","flg22&chitin","sys&flg22&chitin"))
#saveRDS("../SFCH/Data/limma/ORA_Upset_ppep.y1_bin1-1.rds",object=ppep.y1)
#saveRDS("../SFCH/Data/limma/ORA_Upset_ppep.y1_bin1.rds",object=ppep.y1)

ppep.y1 <- readRDS("../SFCH/Data/limma/ORA_Upset_ppep.y1_bin1.rds")%>%filter(pvalue<0.05)
#ppep.y1 <- readRDS("../SFCH/Result/ORA_Upset_ppep.y1.rds")%>%filter(pvalue<0.05,sig_cluster=="sys only"|sig_cluster=="chitin only",bin.name=="bin 18 Protein modification")
table(ppep.y1$bin.name)
bin_count <- plyr::count(ppep.y1, vars =c("bin.name","sig_cluster","pvalue"))%>%
             mutate(bin.code=word(bin.name,1,2,sep = fixed(" ")))%>%filter(pvalue<0.05)%>%arrange(freq)
#levels(bin_count$sig_cluster) <- c("S","F","C","SF","SC","FC","SFC")
## change order to align with the upset plot :flg22_only ahead
bin_count$sig_cluster <- factor(bin_count$sig_cluster,levels =c("flg22 only","sys only","chitin only","sys&chitin","sys&flg22","sys&flg22&chitin","flg22&chitin"))
levels(bin_count$sig_cluster) <- c("F","S","C","SC","SF","SFC","FC")
ggplot(bin_count,aes(sig_cluster,bin.name))+geom_point(aes(color=pvalue,size=freq),show.legend = T)+#scale_y_discrete(limits=rev)+
  scale_color_gradientn(name ="pvalue",colours=c("#371ea3","#46bac2", "#b3eebe"),guide=guide_colorbar(reverse=T))+
  scale_size_area(max_size = 8)+ #guides(size="none")+
  labs(y='',x=' ',size='P_sites number')+
  theme(axis.title=element_text(size=10,colour = 'black'), axis.text=element_text(size=10,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 0.5,fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA)) # get rid of legend panel bg
ggsave("../SFCH/Figure/Diss_chap1/ORA_Upset_bin1.pdf",width=6,height =3.5)

dev.off()




