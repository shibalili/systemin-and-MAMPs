rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(ggplot2)
library(hypeR)
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
back.ID <- unique(summary_df$Protein.ID)# 2119
anno <- readRDS("../SFCH/Data/annotation/SFCH_bin1.rds")

sp <- readRDS("../SFCH/Data/sp_sys1563.rds") # sys: cluster2&5 both enrich bin18,24 
gl <- list()
for(x in 1:7){gl[[x]]<- sp%>%filter(sys_cluster==x)%>%pull(Protein.ID)%>%unique()}

sp <- readRDS("../SFCH/Data/sp_flg1332.rds")
gl <- list()      
for(x in 1:7){gl[[x]]<- sp%>%filter(flg22_cluster==x)%>%pull(Protein.ID)%>%unique()}

sp <- readRDS("../SFCH/Data/sp_chi1499.rds")
gl <- list()      
for(x in 1:7){gl[[x]]<- sp%>%filter(chitin_cluster==x)%>%pull(Protein.ID)%>%unique()}

###### ORA loop ###########
names(gl) <- c("k1","k2","k3","k4","k5","k6","k7")
lengths(gl)
IDs.sum <- apply(anno[,-c(1:4)], 2, sum)
str(IDs.sum)
relevant.bin <-names(IDs.sum)[IDs.sum>0]#how many bins(functional groups) in background
gl_bin <- list()
for(k in 1:7){
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
for (m in 1:7) {hyp_obj[[m]] <- hypeR(signature=gl[[m]], genesets=b[[m]],test = "hypergeometric",background=length(back.ID))}
names(hyp_obj) <- names(gl) 
names(hyp_obj)
# briefly check ORA results
dd=1
hyp_dots(hyp_obj[[dd]],val = "pval",abrv = 35,title = paste("cluster",dd))
#hyp_show(hyp_obj[[1]])
dev.off()

hyp_to_excel <- NULL
for (z in 1:7) {
  #hyp_to_excel$cluster=z
  df <- data.frame(hyp_obj[[z]]$data,cluster=z)
  hyp_to_excel <- rbind(hyp_to_excel,df)}
rownames(hyp_to_excel) <- NULL
y <- NULL
for (n in 1:7) {
  s <- strsplit(hyp_obj[[n]]$data[["hits"]], split = ",")
  y <-rbind(y,data.frame(sig_cluster=rep(names(hyp_obj)[n], sum(lengths(s))),
                         bin.name= rep(hyp_obj[[n]]$data[["label"]], sapply(s, length)),
                         pvalue= rep(hyp_obj[[n]]$data[["pval"]], sapply(s, length)),
                         Protein.ID = unlist(s)%>%gsub("\\s+","", .))) }# remove whitespace in a string
summary_df <-summary_df%>%filter(ppep.site%in%sp$ppep.site)
ppep.y1 <- NULL
for (r in 1:7) {
  site_kn <- summary_df%>%filter(chitin_cluster==r) #!!!!
  y_kn <- y%>%filter(sig_cluster==names(gl)[r])
  ppep.y <- subset(site_kn,site_kn$Protein.ID %in% y_kn$Protein.ID) ##map protein ID to ppep ID
  ppep.y1 <- rbind(ppep.y1,merge(y_kn,ppep.y,by ="Protein.ID"))%>%arrange(sig_cluster,pvalue)}

#saveRDS("../SFCH/Data/sp_sys1563_ora.rds",object=ppep.y1)
#saveRDS("../SFCH/Data/sp_flg1332_ora.rds",object=ppep.y1)
#saveRDS("../SFCH/Data/sp_chi1499_ora.rds",object=ppep.y1)

###### prepare Top5 rank bin #########
sp_s <- readRDS("../SFCH/Data/sp_sys1563_ora.rds")%>%select(1:4,6:8,10) #1338
sp_f <- readRDS("../SFCH/Data/sp_flg1332_ora.rds")%>%select(1:4,6:8,10) #1166
sp_chi <- readRDS("../SFCH/Data/sp_chi1499_ora.rds")%>%select(1:4,6:8,10) #1166

top <- NULL
for (uu in c("k1","k2","k3","k4","k5","k6","k7")) {
  sg <- sp_f %>%filter(sig_cluster==uu) %>%select(c(2:4))%>%unique()%>%
    mutate(Rank=paste0(1:length(bin.name)," (",length(bin.name),")"),
           Rank_short=1:length(bin.name),sig_treatment="flg22",.before = 1)
  sgg <-merge(sp_f %>%filter(sig_cluster==uu),sg%>%select(c(1:3,5)),by = "bin.name" )
  top <- rbind(top,sgg)%>%arrange(sig_cluster,pvalue)}
#saveRDS("../SFCH/Data/sp_sys1563_ora_rank.rds",object=top)
#saveRDS("../SFCH/Data/sp_flg1332_ora_rank.rds",object=top)
#saveRDS("../SFCH/Data/sp_chi1499_ora_rank.rds",object=top)

###################### plot: ORA - plot PTI class ##############################
Top_s <- readRDS("../SFCH/Data/sp_sys1563_ora_rank.rds")
Top_f <- readRDS("../SFCH/Data/sp_flg1332_ora_rank.rds")
Top_chi <- readRDS("../SFCH/Data/sp_chi1499_ora_rank.rds")
sp <- rbind(Top_s,Top_f,Top_chi)%>%filter(Rank_short<=5)
sp$super_class <- as.character(sp$super_class) # remove factor
table(sp$bin.name)
#### plot y-axis: bin groups and rank #####
sp_df <- plyr::count(sp, vars =c("bin.name","sig_cluster","Rank_short","pvalue","sig_treatment"))%>%
  mutate(bin.code=word(bin.name,1,2,sep = fixed(" ")))%>%arrange(freq)%>%
  filter(!bin.name%in%c("bin 35 not assigned","bin 15 RNA biosynthesis","bin 12 Chromatin organisation","bin 7 Coenzyme metabolism",
                        "bin 27 Multi-process regulation","bin 22 Vesicle trafficking","bin 20 Cytoskeleton organisation",
                        "bin 1 Photosynthesis","bin 16 RNA processing","bin 25 Nutrient uptake","bin 14 DNA damage response","bin 5 Lipid metabolism",
                        "bin 2 Cellular respiration","bin 6 Nucleotide metabolism","bin 4 Amino acid metabolism","bin 3 Carbohydrate metabolism"))
sp_df$sig_treatment <- factor(sp_df$sig_treatment,levels = c("systemin","flg22","chitin"))
ggplot(sp_df,aes(sig_cluster,bin.name))+geom_point(aes(color=Rank_short,size=freq),show.legend = T)+scale_y_discrete(limits=rev)+scale_x_discrete(labels=c("1","2","3","4","5","6","7"))+
  facet_grid(~as.factor(sig_treatment))+
  scale_color_gradientn(name ="Rank",colours=c("#371ea3","#46bac2", "#b3eebe"),guide=guide_colorbar(reverse=T))+
  scale_size_area(max_size = 8)+labs(y='',x='Phosphorylation clusters',size='P_sites number')+
  theme(axis.title=element_text(size=10,colour = 'black'), axis.text=element_text(size=10,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 0.5,fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA)) # get rid of legend panel bg
dev.off()


#### plot y-axis: PTI protein groups #####
#### seperate MAPK from kinase ####
mapk <- sp%>%filter(str_detect(Protein_class,"MAPK"))%>%pull(Protein.ID)%>%unique()
sp[which(sp$Protein.ID %in% mapk),"super_class"] <- "MAPK"
sp_df <- plyr::count(sp, vars =c("super_class","sig_cluster","sig_treatment"))%>%filter(!super_class%in%c("0","TPR","ETI","MLO","cyclin","NLR","14-3-3","cellwall"))
sort(table(sp_df$super_class),decreasing = T) # 8 components
sp_df$super_class <- factor(sp_df$super_class,levels = c("Kinase","PPase","Transporter","ROS","Ca2+","MAPK","TF","hormone"))
sp_df$sig_treatment <- factor(sp_df$sig_treatment,levels = c("systemin","flg22","chitin"))

ggplot(sp_df,aes(sig_cluster,super_class))+geom_point(aes(color=freq,size=freq),show.legend = T)+scale_y_discrete(limits=rev)+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7"))+facet_grid(~as.factor(sig_treatment))+
  guides(color= guide_legend(), size=guide_legend())+scale_size_area(max_size = 8)+
  labs(y='',x='Phosphorylation clusters',size='P_sites number',color='P_sites number')+
  theme(axis.title=element_text(size=10,colour = 'black'), axis.text=element_text(size=10,colour = 'black'), 
        panel.background = element_rect(color='black',linewidth = 0.5,fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_line(colour = "gray",linetype = "dotted"),
        panel.grid.major.x = element_blank(),
        legend.key = element_blank(),legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent",color = NA)) # get rid of legend panel bg
#ggsave("../SFCH/Figure/Diss_chap1/elicitor_sp_ORA_pti.pdf",width=6.6,height =2.5)
dev.off()

