rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
summary_df <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_df.rds")
count_sp <- readRDS("../SFCH/Data/Count_specific_at_SFCH4701.rds")%>%filter(!is.na(Count_specific_at))# count specificity
############### systemin specific : black, brown, yellow/ cluster 2 and 5 ##################
## cluster specificity : this is rather loose criteria, further filter is needed!
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")%>%
                filter(sys!=flg22,sys!=chitin,sys!=H2O) # 1132
## DEG ##
sys <- summary_df%>%filter(specify=="S") # 400
## count ##
count <- count_sp%>%filter(Count_specific_at=="systemin")
# count+cluster+DEG = sys_sp 1563 P_sites
sys_sp <-summary_df[unique(c(distribution$ppep.site,sys$ppep.site,count$ppep.site)),] 
table(sys_sp$specify)
#saveRDS("../SFCH/Data/sp_sys1563.rds",object=sys_sp) # 1308 with cluster profile

############### flg22 specific :mainly red/cluster1 ##################
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")%>%
                filter(sys!=flg22,flg22!=chitin,flg22!=H2O) # 1005
flg22 <- summary_df%>%filter(specify=="F") # 415
count <- count_sp%>%filter(Count_specific_at=="flg22")#70
# count+cluster+DEG = flg22_sp 1332 P_sites
flg22_sp <- summary_df[unique(c(distribution$ppep.site,flg22$ppep.site,count$ppep.site)),]
table(flg22_sp$specify)
#saveRDS("../SFCH/Data/sp_flg1332.rds",object=flg22_sp) # 1155 with cluster profile

############### chitin specific  ##################
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")%>%
                filter(sys!=chitin,flg22!=chitin,chitin!=H2O) # 1010
chitin <- summary_df%>%filter(specify=="C") # 282
count <- count_sp%>%filter(Count_specific_at=="chitin")#420
# count+cluster+DEG = chitin_sp 1499 P_sites
chitin_sp <- summary_df[unique(c(distribution$ppep.site,chitin$ppep.site,count$ppep.site)),] 
#saveRDS("../SFCH/Data/sp_chi1499.rds",object=chitin_sp) # 1146 with cluster profile

############### SFC share ##################
share <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")%>%filter(treatment_share=="3 treatments")%>%pull(P_site_ID)%>%unique() # 56
distribution <- readRDS("../SFCH/Data/clustering/SFCH_a_k7_NewPS2+MF1_threshold_Mfuzz_time_arg_order_dist.rds")%>%
  filter(sys!=H2O,flg22!=H2O,chitin!=H2O)%>%filter(sys==flg22,sys==chitin) #109
ns <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds")%>%filter(treatment_share!="3 treatments")%>%pull(P_site_ID)%>%unique()
# cluster-(nsDEG+count)+shareDEG = share_sp 91 P_sites
share_sp <- summary_df%>%filter(ppep.site%in%unique(c(share,setdiff(distribution$ppep.site,c(count_sp$ppep.site,ns)))))
table(share_sp$specify)
#saveRDS("../SFCH/Data/share91_sfc.rds",object=share_sp)


###### plot vertical cluster profile : col:treatment;row:1-7 clusters ###########
sp <- readRDS("../SFCH/Data/sp_sys1563.rds")
table(sp$sys_cluster)
sp <- readRDS("../SFCH/Data/sp_flg1332.rds")
table(sp$flg22_cluster)
sp <- readRDS("../SFCH/Data/sp_chi1499.rds")
table(sp$chitin_cluster)

## share ##
sp <- readRDS("../SFCH/Data/share91_sfc.rds")


Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%filter(ppep.site%in%sp$ppep.site)%>%dplyr::select(2:25)

# systemin specific
Input_df<- data.frame(Input)%>%mutate(gene_id = row.names(.))%>%pivot_longer(-gene_id)%>%mutate(module = summary_df[gene_id,]$sys_cluster)
# flg22 specific
Input_df<- data.frame(Input)%>%mutate(gene_id = row.names(.))%>%pivot_longer(-gene_id)%>%mutate(module = summary_df[gene_id,]$flg22_cluster)
# chitin specific
Input_df<- data.frame(Input)%>%mutate(gene_id = row.names(.))%>%pivot_longer(-gene_id)%>%mutate(module = summary_df[gene_id,]$chitin_cluster)

table(Input_df$module)/24
Input_df$treat <- substr(Input_df$name,1,1)
Input_df$treat <- factor(Input_df$treat,levels = c("S","F","C","H"))
Input_df$time <- substr(Input_df$name,2,3)
Input_df$time <- factor(Input_df$time,levels = c("0","1","2","5","15","45"))
Input_df$name <- factor(Input_df$name,levels = unique(Input_df$name))
Input_df$module <- factor(Input_df$module)

se <- function(x) sqrt(var(x)/length(x))#SE calculation
Percluster <- Input_df %>% group_by(module, treat, time) %>% summarize(mean=mean(value), SE=se(value))
ggplot(data = Input_df, aes(x = time, y = value,group=gene_id))+labs(y= "Normalized intensity", x = "Time (min)")+
  theme(legend.position = "none")+scale_color_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A","#E7B800"),labels=c("Systemin","Flg22","Chitin","H2O"))+
  geom_line(aes(color=treat),linewidth =0.6,alpha=0.02)+
  facet_grid(as.factor(module)~as.factor(treat))+
  geom_line(Percluster,mapping=aes(y =mean,group=interaction(module,treat)),alpha=0.5)+
  geom_errorbar(Percluster,mapping=aes(x = time,y =mean,ymin=mean-SE, ymax=mean+SE,group=interaction(module,treat)),width=0.2)+
  theme(axis.text = element_text(colour = "black",size=13),
        axis.title=element_text(size=16,colour = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.text= element_text(size =14,colour = "black"),
        panel.grid = element_blank(),strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(fill = "transparent",color = NA),
        legend.title=element_blank(),
        legend.text=element_blank())+
  geom_text(data=plyr::count(Input_df, vars =c("treat","module")),
            aes(x=1.8, y=1.8, label=paste0("N=",freq/6)),colour="black", inherit.aes=FALSE, parse=FALSE,size = 3.5)
#ggsave("../SFCH/Figure/Diss_chap1/Verti_sys_sp_plot.pdf",width =5.6,height =7.5)
#ggsave("../SFCH/Figure/Diss_chap1/Verti_flg22_sp_plot.pdf",width =5.6,height =7.5)
#ggsave("../SFCH/Figure/Diss_chap1/Verti_chitin_sp_plot.pdf",width =5.6,height =7.5)
dev.off()

############# plot: the heatmap of vertical cluster profile ####################   
library(pheatmap)
library(viridis)#viridis color
save_pheatmap_pdf <- function(x, filename, width=4, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()}
# systemin specific
Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%filter(ppep.site%in%sp$ppep.site)%>%mutate(module=summary_df[ppep.site,]$sys_cluster)%>%arrange(module)%>%dplyr::select(2:26)
# flg22 specific
Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%filter(ppep.site%in%sp$ppep.site)%>%mutate(module=summary_df[ppep.site,]$flg22_cluster)%>%arrange(module)%>%dplyr::select(2:26)
# chitin specific
Input <- readRDS("../SFCH/Data/clustering/newdf_scale.rds")%>%filter(ppep.site%in%sp$ppep.site)%>%mutate(module=summary_df[ppep.site,]$chitin_cluster)%>%arrange(module)%>%dplyr::select(2:26)

k <- Input%>% select(1:24)
p1 <- pheatmap(k,color=viridis(3),main = "",cluster_rows = F,cluster_cols=F,show_rownames = FALSE,gaps_col=c(6,12,18),
               labels_col=colnames(k),angle_col="45",fontsize = 5,legend = T)#gaps_row=gaps_row

#save_pheatmap_pdf(p1, "../SFCH/Figure/Diss_chap1/pheatmap_sys_sp.pdf")
#save_pheatmap_pdf(p1, "../SFCH/Figure/Diss_chap1/pheatmap_flg22_sp.pdf")
#save_pheatmap_pdf(p1, "../SFCH/Figure/Diss_chap1/pheatmap_chitin_sp.pdf")
#save_pheatmap_pdf(p1, "../SFCH/Figure/Diss_chap1/pheatmap_share_sp.pdf")












