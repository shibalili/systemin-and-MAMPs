rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)#list_to_matrix
library(rlang)#parse_quo
library(scales)#labels comma
library(ggpubr) #get_legend
library(cowplot)#draw plot
library(ggforce) # for theme_no_axes
a <- readRDS("../SFCH/Data/limma/Rem_water_SFC_a.rds") # exclude water
sys <- a%>%filter(treatment=="systemin")
flg22 <- a%>%filter(treatment=="flg22")
chitin <- a%>%filter(treatment=="chitin")
list <- list(systemin=unique(sys$P_site_ID),flg22=unique(flg22$P_site_ID),chitin=unique(chitin$P_site_ID))
binary <- as.data.frame(list_to_matrix(list))%>%rownames_to_column("P_sites")#turn categorical value to binary
binary_long <- binary%>%pivot_longer(!P_sites,names_to=c("group"), values_to="binary")
#generates all possible combinations of A&B&C, A&B&!C, A&!B&C, etc, as strings, then reduces them to non-redundant set
combins <- unique(binary_long$group) %>%c(paste0("!", .)) %>%combn(length(.)/2) %>%t() %>%
  as.data.frame() %>%
  filter(apply(., 1,function(x) length(unique(gsub("!", "", x)))== ncol(.)&!(length(grep("!",x))%in% c(0,ncol(.))))) %>%
  unite("expressions",names(.), sep = " & ")
combins <- rbind(combins,paste(unique(binary_long$group),sep = ' & ', collapse=' & '))
### attach setsize
filter_sets = function(filter_expr){
  binary_long %>%spread(group, binary) %>%filter(!!parse_quo(env = global_env(),filter_expr))%>%nrow()
}
combins$value = as.numeric(sapply(combins$expressions, filter_sets))
separated <- combins %>% separate(expressions, into = c('a','b','c'), sep=' & ')#prep for deviation calculation
###different set:number induced by each individually regardless of overlaps
setsizes <- data.frame(systemin=as.numeric(length(unique(sys$P_site_ID))),flg22=as.numeric(length(unique(flg22$P_site_ID))),chitin=as.numeric(length(unique(chitin$P_site_ID))))
n <- as.numeric(nrow(binary)) #total number DEG
# deviation from set size predicted by random mixing
calculate_deviation = function(elicitors){
  use <- colnames(setsizes)[colnames(setsizes) %in% elicitors]
  yes <- setsizes[,use]
  yes <- yes/n
  yes <- prod(yes)
  dont <- colnames(setsizes)[!colnames(setsizes) %in% elicitors]
  if(length(dont) >0){
    no <- setsizes[,dont]
    no <- no/n
    no <- 1-no
    no <- prod(no)#multiplication results
  }else{no <- 1}
  product <- prod(yes, no)#expected value
  cardinality <- as.numeric(elicitors['value'])/n #exclusive intersection setsize value
  deviation <- cardinality-product ##deviation:the relative sizes of these sets departed from that of a random assortment of genes among treatments
  deviation
}
separated$deviation <- apply(separated, 1, calculate_deviation)
recombined <- separated %>% unite(expression, a,b,c)
#for a given elicitor, make a binary presence/absence column/ make a column for each elicitor: 1 if included in set, 0 if not
elicitor_columns=function(deviationdata, elicitor){
  deviationdata$new <- 1
  deviationdata$new[grep(paste0("!",elicitor),deviationdata$expression)] <- 0
  colnames(deviationdata)[colnames(deviationdata)=="new"] <- elicitor
  return(deviationdata)}
doit <- elicitor_columns(recombined,"systemin")
doit <- elicitor_columns(doit,"flg22")
doit <- elicitor_columns(doit,"chitin")
doit$degree <-  rowSums(doit[4:6])#calculate the degree
doit <- arrange(doit, desc(value))
doit$start_angle <- seq(0,2*pi-(2*pi/7), length.out=7)#ggforce needs start and stop agle of each "chunk" of pie, then calculate and add
doit$end_angle <- seq(2*pi/7, 2*pi, length.out=7)
#this is mostly just for making the fills different colors
doit$systemin[doit$systemin==1] <- "systemin"
doit$flg22[doit$flg22==1] <- "flg22"
doit$chitin[doit$chitin==1] <- "chitin"
sort(doit$deviation,decreasing = T)##set legend bar range
#make plots for vertical set sizes bars
vertical <- ggplot(doit[1:7,], aes(x=as.factor(start_angle), y=value, fill=deviation))+
  geom_bar(stat="identity", color="grey", linewidth=0.2)+geom_text(aes(label=value), vjust=-0.3, size=3)+
  theme_classic(base_size=10)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),legend.position="none",
        rect=element_rect(fill="transparent",color=NA), plot.background = element_blank(), panel.background = element_blank())+
  scale_y_continuous(label=comma, expand=c(0,0), limits=c(0, 500), name="P_sites in set")+
  scale_fill_distiller(palette="PRGn", direction=-1, limits=c(-0.2, 0.2))
vert <- ggplot(doit[1:7,], aes(x=as.factor(start_angle), y=value, fill=deviation))+labs(fill=expression("Deviation\n"))+geom_bar(stat="identity")+
  scale_fill_distiller(palette="PRGn",direction=-1,limits=c(-0.2, 0.2),breaks = c(-0.2,-0.1,0,0.1,0.2),
                       guide=guide_colorbar(barwidth=unit(2.5, "cm"),barheight=unit(0.25, "cm")))+theme_minimal(base_size=6)+
  theme(legend.position = "top",legend.title = element_text(vjust =0.2,color="black",size =8),legend.text=element_text(size=6),legend.key.size=unit(0.5,"cm"))
deviat <- ggpubr::get_legend(vert)#%>%as_ggplot() # only the first legend is returned.
#make plot for horizontal elicitor sizes bars
setsizes_long <-setsizes%>% pivot_longer(everything(),names_to="elicitor", values_to="ngenes")
setsizes_long$elicitor <- factor(setsizes_long$elicitor,levels=c("systemin","flg22","chitin"))
horizontal <- ggplot(setsizes_long, aes(x=ngenes/100, y=reorder(elicitor, -ngenes), fill=elicitor))+
  geom_bar(stat="identity")+theme_classic(base_size=8)+
  theme(axis.title.y=element_blank(), axis.text.y = element_blank(), 
        legend.position="none", axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        plot.margin = margin(0,0,0,0, "pt"),
        panel.border = element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_reverse(expand=c(0,0), limit=c(8,0), name="Responsive P_sites (x100)")+
  scale_fill_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A"),labels=c("Systemin","Flg22","Chitin"))
#scale_fill_brewer(palette = "Set1")
#make a simple version, just to grab the legend
horiz <- ggplot(setsizes_long, aes(x=ngenes/100, y=reorder(elicitor, -ngenes), fill=elicitor))+
  geom_bar(stat="identity")+labs(fill="Treatment")+
  scale_fill_manual(values= c("#E41A1C", "#377EB8", "#4DAF4A"),labels=c("Systemin","Flg22","Chitin"))+
  theme_minimal(base_size=10)+
  theme(legend.key.size=unit(0.35,"cm"),
        panel.border = element_blank(),panel.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank())
elic <- ggpubr::get_legend(horiz)%>%as_ggplot()
#make plots for which elicitor sets dots
doit_long <- doit %>% pivot_longer(c(systemin, flg22, chitin), names_to="elicitor", values_to="present")
doit_long$elicitor <- factor(doit_long$elicitor, levels=levels(reorder(setsizes_long$elicitor, -setsizes_long$ngenes)))
doit_long$present <- factor(doit_long$present, levels = c("systemin","flg22","chitin","0"))
definedvibrant <- c("0"="white", "systemin"="#E41A1C", "flg22"="#377EB8", "chitin"="#4DAF4A")
dots <- ggplot(doit_long[1:nrow(doit_long),], aes(x=start_angle, y=elicitor, fill=present))+
  geom_point(shape=21, size=3.5,stroke = 1)+scale_fill_manual(values=definedvibrant)+
  theme_no_axes()+theme(panel.background = element_blank(), plot.background = element_blank(), panel.border=element_blank(),
                        legend.position="none", plot.margin = margin(0,0,0,0, "pt"))
#Combine all the bits
upset <- ggdraw()+
  draw_plot(vertical,x=0.25, y=0.38, width=0.5, height=0.62)+draw_plot(horizontal,x=0.02, y=0.08, width=0.3, height=0.3)+
  draw_plot(dots, x=0.345, y=0.15, height=0.24, width=0.38)+
  draw_plot(deviat, x=0.65, y=0.85, height=0.1, width=0.2)+
  draw_plot(elic, x=0.03, y=0.55, height=0.05, width=0.2)
#deviation color need to adjust in ilustrator
#ggsave("../SFCH/Figure/Diss_chap1/upset.pdf",width =5,height =2.8) 
dev.off()
