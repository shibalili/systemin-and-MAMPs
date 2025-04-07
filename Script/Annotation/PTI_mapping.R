rm(list = ls())
setwd("~/Documents/HOH_260/EXP/Phospho/SFCH/")
getwd()
library(tidyverse)
library(magrittr)#%$%
##### Prepare datafram:add Protein_class to annotation used in this dataset #####
mapping <- readRDS("../SFCH/Data/annotation/phytozome_ITAG3.2_standard.rds")%>%
  mutate(simp.id=str_split_i(Protein.ID, '\\.', 1),Protein_class=0,.after = Protein.ID)
# subset annotation for my dataset,4701 P_sites mapping 2119 protein 
df_PSMF <- readRDS("../SFCH/Data/annotation/SFCH_summary4701_PSMF_threshold_Mfuzz_df.rds")%>%
  mutate(simp.id=str_split_i(Protein.ID, '\\.', 1),Protein_class=0,.before = Protein.ID)
mapping <-mapping%>%filter(simp.id%in%df_PSMF$simp.id) #2119

### Protein_class reference list ####
RLK2012 <- readRDS("../SFCH/Data/pti/RLKs.rds") # Sakamoto2012:contain all RLKs

# HMMER results(from chirs):for transmembrane domains ##
hmm_all <- read_table("../SFCH/Data/pti/chris/degs_sigonly_hmmscan.txt", col_names = FALSE, skip = 3) %$% .[1:(nrow(.)-10),]
hmm_all <- hmm_all[hmm_all$X5 < sort(hmm_all$X5)[15916],]  # calculated from chris::brkpnt_fun

# PPases ###############################################################################################################
# Solyc10g005640.3.1 AtCRPK1 homologue, COLD-RESPONSIVE PROTEIN KINASE 1 : very interesting, harbour both Kinase and PP2C PPase domains 
pp2c <- read.delim("../SFCH/Data/pti/SlPP2C.txt") # 13 hits
pp2c$simp.id <- gsub("\\s+", "", pp2c$simp.id)
pp2c$name <- gsub("\\s+", "", pp2c$name)
rownames(pp2c) <- pp2c$simp.id
mapping[which(mapping$simp.id %in% pp2c$simp.id),"Protein_class"] <- paste0("PPase:",pp2c[mapping[which(mapping$simp.id %in% pp2c$simp.id),"simp.id"],"name"])

ppp_ref <- read.delim("../SFCH/Data/pti/SlPPP.txt") #16 hits+1 hit against Protein phosphatase inhibitor 2
rownames(ppp_ref)<- ppp_ref$simp.id
mapping[which(mapping$simp.id %in% ppp_ref$simp.id),"Protein_class"] <- paste0("PPase:",ppp_ref[mapping[which(mapping$simp.id %in% ppp_ref$simp.id),"simp.id"],2])
## PPKL (protein phosphatase with Kelch-like repeat domains) e.g. BSU (see below)

# RLKs #################################################################################################################
## CRINKLY4-like: no hits
## CrRLKL1 - malectins:8 hits=3+5
malec_ref <- read.delim("../SFCH/Data/pti/SlCrRLK1.txt")
rownames(malec_ref) <- malec_ref$simp.id
mapping[which(mapping$simp.id %in% malec_ref$simp.id),"Protein_class"] <- malec_ref[mapping[which(mapping$simp.id %in% malec_ref$simp.id),"simp.id"],2]
malectin<-c(read.delim("../SFCH/Data/pti/SlCrRLK1.txt")[[1]],RLK2012$locus[RLK2012$RLK_subfamily=="CrRLK"],hmm_all$X3[hmm_all$X1 == "Malectin"],
            mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                                any_vars(str_detect(., pattern="CrRLK")))%>%unique()%>%pull(simp.id))%>%unique()%>%sort()
malectin<-mapping[which(mapping$simp.id %in% malectin),"simp.id"] #8
mapping[which(mapping$simp.id %in% malectin),"Protein_class"] <- ifelse(mapping[which(mapping$simp.id %in% malectin),"Protein_class"]!=0,
                                                                        mapping[which(mapping$simp.id %in% malectin),"Protein_class"] <- mapping[which(mapping$simp.id %in% malectin),"Protein_class"],
                                                                        mapping[which(mapping$simp.id %in% malectin),"Protein_class"] <- "SlCrRLK1L")%>%paste0("RLK:",.)

## DUF26 domain:5 hits
# DUF26-containing cell surface receptors, e.g.Cysteine-rich Receptor-like Kinases (CRKs),are predicted to act as carbohydrate receptors
duf26<-c(mapping%>%filter(simp.id%in% unique(RLK2012$locus[RLK2012$RLK_subfamily == "DUF26"]))%>%pull(simp.id),
         mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                             any_vars(str_detect(., fixed("DUF26",ignore_case=TRUE))))%>%pull(simp.id))%>%unique()%>% sort()
mapping[which(mapping$simp.id %in% duf26),"Protein_class"] <-RLK2012[duf26,"RLK_subfamily"]%>%paste0("RLK:",.)
## Lectin receptor kinases (LecRKs):C-type LecRKs (C-LecRKs), L-type LecRKs (L-LecRKs) and G-type LecRKs (G-LecRKs):5 hits
lec_ref <- read.delim("../SFCH/Data/pti/SlLecRKs.txt")
rownames(lec_ref) <- lec_ref$simp.id
lec_ref$name <- paste0("SlLecRK-",word(lec_ref$Type,sep = fixed("-")))
lec <- mapping%>%filter(simp.id%in% unique(lec_ref$simp.id))%>%pull(simp.id) #5
mapping[which(mapping$simp.id %in% lec),"Protein_class"] <- lec_ref[mapping[which(mapping$simp.id %in% lec),"simp.id"],2]%>%paste0("RLK:",.)
## LEAF RUST 10 DISEASE-RESISTANCE LOCUS RECEPTOR-LIKE PROTEIN KINASE-LIKE2
## no clear HMMER hits, only hit is a WAK
## wall-associated kinases:all included : 4 hits
lrk10l2<- RLK2012$locus[str_detect(RLK2012$RLK_subfamily,pattern = "LRK1" )] %>% unique() %>% sort()
mapping[which(mapping$simp.id %in% lrk10l2),"Protein_class"] <- RLK2012[mapping[which(mapping$simp.id %in% lrk10l2),"simp.id"] ,"RLK_subfamily"]%>%paste0("RLK:",.)
#Solyc12g036325.1.1 in ITAG3.2= Solyc12g036330.1.1 in ITAG4
mapping["Solyc12g036325.1.1","Protein_class"] <- RLK2012["Solyc12g036330","RLK_subfamily"]%>%paste0("RLK:",.)

## LysM:5 hits
lyk_ref <- read.delim("../SFCH/Data/pti/SlLYKs.txt")
rownames(lyk_ref)<- lyk_ref$simp.id
lyk_ref$name<-gsub("\\s+", "", lyk_ref$name)
mapping[which(mapping$simp.id %in% lyk_ref$simp.id),"Protein_class"] <- lyk_ref[mapping[which(mapping$simp.id %in% lyk_ref$simp.id),"simp.id"] ,2]%>%paste0("RLK:",.)

## PERK(Extensin):9 hits
perk <- c(mapping%>%filter(simp.id%in% unique(RLK2012$locus[str_detect(RLK2012$RLK_subfamily,pattern = "Extensin" )]))%>%pull(simp.id),
          mapping%>%filter(simp.id%in% unique(RLK2012$locus[str_detect(RLK2012$RLK_subfamily,pattern = "PERK" )]))%>%pull(simp.id))%>%
  unique()%>%sort()
mapping[which(mapping$simp.id %in% perk),"Protein_class"] <- RLK2012[mapping[which(mapping$simp.id %in% perk),"simp.id"] ,"RLK_subfamily"]%>%paste0("RLK:",.)
ext_ref <- read.delim("../SFCH/Data/pti/SlEXTs.txt")
rownames(ext_ref)<- ext_ref$simp.id
ext <- setdiff(mapping%>%filter(simp.id%in% ext_ref$simp.id)%>%pull(simp.id)%>%sort(),perk) # extra 2 hits:Formin-homolog (FH) chimeric EXTs
mapping[which(mapping$simp.id %in% ext),"Protein_class"] <- ext_ref[mapping[which(mapping$simp.id %in% ext),"simp.id"],2]%>%paste0("RLK:",.)

## SD:4 hits but 3 are Leck, one is Solyc02g079590;one hit here
sd <- mapping%>%filter(simp.id%in% unique(RLK2012$locus[str_detect(RLK2012$RLK_subfamily,pattern = "SD" )]))%>%pull(simp.id)
mapping[which(mapping$simp.id %in% setdiff(sd,lec_ref$simp.id)),"Protein_class"] <- "RLK:SD-1"
## FOR stem diameter:NET/KIP1: kinase interacting (KIP1-like) family protein:
# Solyc09g082510 SlSD1 (same abbrev different finction)
kip <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(., fixed("NET1",ignore_case=TRUE))))%>%pull(simp.id)%>%sort()
mapping[which(mapping$simp.id %in% kip),"Protein_class"] <- "RLK:NET1"

### LRR-RLK ##########
## Receptor
# STRUBBELIG (SUB) also known as SCRAMBLED (SCM);LRR-V/STRUBBELIG-RECEPTOR FAMILY (SRF)
srf_ref <- read.delim("../SFCH/Data/pti/SlSRFs.txt")#4 hits+1 hit=5 hits
rownames(srf_ref) <- srf_ref$simp.id
srf <- mapping%>%filter(simp.id%in% srf_ref$simp.id)%>%pull(simp.id)
mapping[which(mapping$simp.id %in% srf),"Protein_class"] <- srf_ref[mapping[which(mapping$simp.id %in% srf),"simp.id"],2]%>%paste0("RLK:",.)
SRF <- setdiff(mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                                   any_vars(str_detect(., fixed("STRUBBELIG",ignore_case=TRUE))))%>%pull(simp.id)%>%sort(),srf)
mapping[which(mapping$simp.id %in% SRF),"Protein_class"] <- "RLK:SlSRF"
# RGF-peptide receptor(RGFR):Solyc09g091400
mapping[which(mapping$simp.id %in% "Solyc09g091400"),"Protein_class"] <- "RLK:SlRGFR"
# PSKR2 :Solyc07g063000
mapping[which(mapping$simp.id %in% "Solyc07g063000"),"Protein_class"] <- "RLK:SlPSKR2"
# BAM1: Solyc02g091840
mapping[which(mapping$simp.id %in% "Solyc02g091840"),"Protein_class"] <- "RLK:SlBAM1"
# SlMIK2 (best hit against AtMIK2:At4g08850,but 44% indentity)
mapping[which(mapping$simp.id %in% "Solyc09g072810"),"Protein_class"] <- "RLK:SlMIK2"
# SlBRI1 Solyc04g051510 (literature) but 2 hits in the datasets have sequence similarity (50% and 47%) with SlBRI1
bri <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(., fixed("brassinosteroid receptor protein",ignore_case=TRUE))))%>%pull(simp.id)%>%sort()
mapping[which(mapping$simp.id %in% bri),"Protein_class"] <- "RLK:SlBRI1"
# IMK2 (inflorescence meristem receptor-like kinase 2):1 nit Solyc09g015170
mapping[which(mapping$simp.id %in% "Solyc09g015170"),"Protein_class"] <- "RLK:SlIMK2"
# SlSIRK1:SUCROSE-INDUCED RECEPTOR KINASE 1 (best hit against AtSIRK1 AT5G10020.1):Solyc04g054200
mapping[which(mapping$simp.id %in% "Solyc04g054200"),"Protein_class"] <- "RLK:SlSIRK1"
# SlGHR1:GUARD CELL HYDROGEN PEROXIDE-RESISTANT1 (best hit against AtGHR1:AT4G20940.1): two hits but sequence identity were ~40%
ghr <- mapping%>%filter(Best.hit.arabi.name=="AT4G20940.1")%>%pull(simp.id)%>%sort()
mapping[which(mapping$simp.id %in% ghr),"Protein_class"] <- "RLK:SlGHR1"
# chitin receptor might be SlLYK4 Solyc02g089900: no hits; co-receptor SlCERK1 Solyc07g049180 is mentioned in lym
# HERK1&2:hercules receptor kinase: 3 hits mentioned in CrRLKL1
# FER:1 hit Solyc09g015830 also mentioned in CrRLKL1
mapping[which(mapping$simp.id %in% "Solyc09g015830"),"Protein_class"] <- "RLK:SlFER1"
# no hits:SlFLS2/SlFLS3/SYR1/SYR2/PORK1
## Co-receptor
# SlSERK3B Solyc01g104970
mapping[which(mapping$simp.id %in% "Solyc01g104970"),"Protein_class"] <- "RLK:SlSERK3B"
# SlBIR2 Solyc02g087460 (best hit against AtBIR2 At3g28450),BIR1/3/4 no hits
mapping[which(mapping$simp.id %in% "Solyc02g087460"),"Protein_class"] <- "RLK:SlBIR2"
# NIK (NSP-INTERACTING KINASE):interacts with nuclear shuttle protein (NSP) of geminivirus during infection:one hit
nik_ref <- read.delim("../SFCH/Data/pti/SlNIKs.txt")
rownames(nik_ref) <- nik_ref$simp.id
nik <- mapping%>%filter(simp.id%in% nik_ref$simp.id)%>%pull(simp.id)# one hit
mapping[which(mapping$simp.id %in% nik),"Protein_class"] <- nik_ref[nik,2]%>%paste0("RLK:",.)
# RKL1:RECEPTOR-LIKE KINASE1
# A LRR-RLK homologous to TARK1 (Tomato Atypical Receptor Kinase1):Solyc11g011020
mapping[which(mapping$simp.id %in% "Solyc11g011020"),"Protein_class"] <- "RLK:SlTARK1L"
# SlPXC1:a regulator of secondary wall formation (best hit against AtPXC1:At2g36570)
mapping[which(mapping$simp.id %in% "Solyc09g008860"),"Protein_class"] <- "RLK:SlPXC1"

# Transmembrane kinase involved in auxin signaling:TMK
# TMK1 (At1g66150): 2 hits; TMK2 (At1g24650); TMK3 (At2g01820) and TMK4 (At3g23750):1 hit
mapping[mapping%>%filter(Best.hit.arabi.name%in%c("AT1G66150.1"))%>%pull(Protein.ID),"Protein_class"] <- "RLK:SlTMK1"
mapping[mapping%>%filter(Best.hit.arabi.name%in%c("AT3G23750.1"))%>%pull(Protein.ID),"Protein_class"] <- "RLK:SlTMK4"
## other LRR-RLKs
lrr2015<-readRDS("../SFCH/Data/pti/LRR-RLK2015.rds")%>%arrange(Gene.ID)%>%pull(1) # Wei2015
dplrr<-read_csv("../SFCH/Data/pti/chris/degs_sigonly_rlk.txt", col_names = FALSE, skip = 448)[[1]] %>% sort()# DeepLRR
lrr_hmm<-intersect(unique(hmm_all$X3[grep("^LRR", hmm_all$X1)]),unique(hmm_all$X3[hmm_all$X1 == "PK_Tyr_Ser-Thr"])) %>% unique() %>% sort()
lrr_unver<-read.table("../SFCH/Data/pti/chris/rlk_lrr_unverified.txt", quote="\"", comment.char="")[[1]] #33
rlk_lrr<-c(lrr2015,dplrr,lrr_hmm,lrr_unver,RLK2012%>%filter(str_detect(RLK_subfamily,pattern="LRR"))%>%pull(1),
           mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                               any_vars(str_detect(., fixed("LRR",ignore_case=TRUE))))%>%unique()%>%pull(simp.id))
# remove above annotated LRR-RLKs
lrr<-mapping[which(mapping$simp.id %in% rlk_lrr),]%>%filter(Protein_class=="0")%>%pull(simp.id)%>%sort()
mapping[which(mapping$simp.id %in% lrr),"Protein_class"] <- "RLK:SlLRR-RLK"
#table(mapping$Protein_class) #2119-1990=129 annotated
# RLPs #################################################################################################################
rlp_ref <- read.delim("../SFCH/Data/pti/SlRLP.txt")[[1]] # 5 hits but 3 hits are AtRLP homologues
rlp <- mapping[which(mapping$simp.id %in% rlp_ref),]%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% rlp$simp.id),"Protein_class"] <- "RLP:SlRLP"
# RLCK (receptor like cytopasmic kinases) ##############################################################################
# SlFir1:Flagellin-sensing 2 (Fls2)/Fls3 interacting RLCK 1:Solyc12g099830 (BSK-type)
mapping[which(mapping$simp.id %in% "Solyc12g099830"),"Protein_class"] <- "RLCK:SlFir1"
# SlACIK1(Avr9/Cf-9 induced kinase 1) or (auxin-regulated dual specificity cytosolic kinase):Solyc06g062920
# PBL12 homologue c("Solyc06g062920")
mapping[which(mapping$simp.id %in% "Solyc06g062920"),"Protein_class"] <- "RLCK:SlACIK1"
# PATTERN-TRIGGERED IMMUNITY (PTI) COMPROMISED RECEPTOR-LIKE CYTOPLASMIC KINASE 1 (AtPCRK1/2 homologue):Solyc09g008010
mapping[which(mapping$simp.id %in% "Solyc09g008010"),"Protein_class"] <- "RLCK:SlPCRK"
## AtPBL homologue
# PBL19/20 homologue c("Solyc08g077560","Solyc01g088690")
mapping[which(mapping$simp.id %in% c("Solyc08g077560","Solyc01g088690")),"Protein_class"] <- "RLCK:SlPBL19"
# ATSIK homologue c("Solyc10g075040"):Encodes an osmotic stress-inducible kinase that functions as a negative regulator of osmotic stress signaling in plants.
mapping[which(mapping$simp.id %in% "Solyc10g075040"),"Protein_class"] <- "RLCK:SlSIK"
# PBL31 homologue c("Solyc04g082500")
mapping[which(mapping$simp.id %in% "Solyc04g082500"),"Protein_class"] <- "RLCK:SlPBL31"
# PBL30 homologue c("Solyc03g032150")
mapping[which(mapping$simp.id %in% "Solyc03g032150"),"Protein_class"] <- "RLCK:SlPBL30"
# PBL2/3/4/18/29 homologue c("Solyc11g062400","Solyc05g007140")
mapping[which(mapping$simp.id %in% c("Solyc11g062400","Solyc05g007140")),"Protein_class"] <- "RLCK:SlPBL3"
# PBL11 homologue c("Solyc09g010850","Solyc10g084770")
mapping[which(mapping$simp.id %in% c("Solyc09g010850","Solyc10g084770")),"Protein_class"] <- "RLCK:SlPBL11"
# PBL1/9/10/BIK1 homologue c("Solyc04g011520","Solyc06g005500")
mapping[which(mapping$simp.id %in% c("Solyc04g011520","Solyc06g005500")),"Protein_class"] <- "RLCK:SlBIK1"
# PBL8 homologue c("Solyc10g074710")
mapping[which(mapping$simp.id %in% "Solyc10g074710"),"Protein_class"] <- "RLCK:SlPBL8"
# PBL34/35/36 homologue c("Solyc01g010780")
mapping[which(mapping$simp.id %in% "Solyc01g010780"),"Protein_class"] <- "RLCK:SlPBL34"
# PBL7 homologue c("Solyc10g085990")
mapping[which(mapping$simp.id %in% "Solyc10g085990"),"Protein_class"] <- "RLCK:SlPBL7"
# PBS1 homologue c("Solyc05g024290")
mapping[which(mapping$simp.id %in% "Solyc05g024290"),"Protein_class"] <- "RLCK:SlPBS1"
# PBS2 homologue c("Solyc11g072290")
mapping[which(mapping$simp.id %in% "Solyc11g072290"),"Protein_class"] <- "RLCK:SlPBS2"
# other RLKs
rlk <- mapping[which(mapping$simp.id %in% RLK2012$locus),]%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% rlk$simp.id),"Protein_class"] <- RLK2012[mapping[which(mapping$simp.id %in% rlk$simp.id),"simp.id"],2]%>%paste0("RLCK:",.)
## RKF3 is LRR-RLK 
mapping$Protein_class[str_detect(mapping$Protein_class,"RKF3")] <-"RLK:RKF3"
# NLRs #################################################################################################################
nlr_ref <- read.csv("../SFCH/Data/pti/chris/NLRs.csv")# 1 hit
rownames(nlr_ref) <- nlr_ref$SequenceID
# 1 hit in nlr_ref Solyc06g008480;3 hits in hmm_all(NB-ARC domain) but not accurate
nlr<- mapping%>%filter(simp.id%in% nlr_ref$SequenceID)
mapping[which(mapping$simp.id %in% nlr$simp.id),"Protein_class"] <- nlr_ref[nlr$simp.id,2]%>%paste0("NLR:",.)
## ETI related 
# SlRIN4(RPM1-INTERACTING PROTEIN4):a guardee protein required for the recognition of bacterial effectors
rin4_ref <-read.delim("../SFCH/Data/pti/SlRIN4s.txt")%>%suppressWarnings()
rownames(rin4_ref) <- rin4_ref$simp.id
rin4 <-  mapping%>%filter(simp.id%in% rin4_ref$simp.id)%>%pull(simp.id)%>%sort() #3 hits
mapping[which(mapping$simp.id %in% rin4_ref$simp.id),"Protein_class"] <- rin4_ref[rin4,2]%>%paste0("ETI:",.)
## SlRIN4-like (Best hit against AtRIN4-like):One hit Solyc07g045340
mapping[which(mapping$simp.id %in% "Solyc07g045340"),"Protein_class"] <- "ETI:SlRIN4-like"

# RBOHs ################################################################################################################
rboh_ref <- read.delim("../SFCH/Data/pti/SlRBOHs.txt")# 3 hits
rownames(rboh_ref)<- rboh_ref$simp.id
mapping[which(mapping$simp.id %in% rboh_ref$simp.id),"Protein_class"] <- rboh_ref[mapping[which(mapping$simp.id %in% rboh_ref$simp.id),"simp.id"],2]%>%paste0("ROS:",.)

# Ca2+ ##########################################################################################################
## CPK/CDPKs(Calcium-dependent protein kinases) & CRK(CDPK-related kinase)
# Solyc03g031670 LeCPK1
cpk_ref <- read.delim("../SFCH/Data/pti/SlCDPK.txt") #12 hits=9 CPK+3 CDPK-related kinases
rownames(cpk_ref)<- cpk_ref$simp.id
mapping[which(mapping$simp.id %in% cpk_ref$simp.id),"Protein_class"] <- cpk_ref[mapping[which(mapping$simp.id %in% cpk_ref$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)
cdpk <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                            any_vars(str_detect(., fixed("cpk",ignore_case=TRUE))))%>%filter(Protein_class=="0") # 2 hits
mapping[which(mapping$simp.id %in% cdpk$simp.id),"Protein_class"] <- "Ca2+:SlCPK"
## CIPKs(calcineurin B-like interacting protein kinase)/SNF1-related protein kinase (SnRK3):4 hits
## CBLs (calcineurin B-like protein) :Sl-Cbl3 one hit
cipk_ref <- read.delim("../SFCH/Data/pti/SlCIPK.txt") # 5 hits
rownames(cipk_ref) <- cipk_ref$simp.id
mapping[which(mapping$simp.id %in% cipk_ref$simp.id),"Protein_class"] <- cipk_ref[mapping[which(mapping$simp.id %in% cipk_ref$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)

## SlCMLs(calmodulin-like protein) & SlCaMs (Calmodulin):5 hits
cmlcam_ref <- read.delim("../SFCH/Data/pti/SlCMLCAM.txt")
rownames(cmlcam_ref) <- cmlcam_ref$simp.id
mapping[which(mapping$simp.id %in% cmlcam_ref$simp.id),"Protein_class"] <- cmlcam_ref[mapping[which(mapping$simp.id %in% cmlcam_ref$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)

## SlCRK (Cysteine-Rich Receptor-Like Protein Kinase) :5 hits
crk_ref <- read.delim("../SFCH/Data/pti/SlCRKs.txt")
rownames(crk_ref) <- crk_ref$simp.id
crk <- mapping%>%filter(simp.id%in%crk_ref$simp.id)%>%pull(simp.id)
mapping[which(mapping$simp.id %in%crk) ,"Protein_class"] <- crk_ref[crk,2]%>%paste0("RLK:",.)

## cyclin-dependent kinases (CDKs):9 hits
cdk <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(., fixed("cdk",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in%cdk$simp.id) ,"Protein_class"] <- "cyclin:SlCDK"
## cyclin genes
cyc_ref <- read.delim("../SFCH/Data/pti/Slcyclin.txt") #5 hits
rownames(cyc_ref) <- cyc_ref$simp.id
mapping[which(mapping$simp.id %in%cyc_ref$simp.id) ,"Protein_class"] <- cyc_ref[mapping[which(mapping$simp.id %in%cyc_ref$simp.id) ,"simp.id"],2]%>%paste0("cyclin:",.)

# MAP kinases ##########################################################################################################
# The MAP kinases are components of a linear cascade of three consecutively acting protein kinases: MAPKK kinases (MAPKKKs), MAPK kinases (MAPKKs), and MAPKs
# 16 SlMPK;5 MPKK;33 MEKK(MAPKKK) + 11 ZIK + 40 RAF
## MAP4K:2 hits
map4k <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., fixed("MAP4K",ignore_case=TRUE))))%>%filter(Protein_class=="0")%>%pull(simp.id)
mapping[which(mapping$simp.id %in% map4k),"Protein_class"]<-"MAPK:SlMAP4K"
## MAP3K:18 hits
mekk_ref <- read.delim("../SFCH/Data/pti/SlMAPKKK.txt")
rownames(mekk_ref) <- mekk_ref$simp.id
mekk <- mapping[which(mapping$simp.id %in% mekk_ref$simp.id),]%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% mekk$simp.id),"Protein_class"] <- mekk_ref[mapping[which(mapping$simp.id %in% mekk$simp.id),"simp.id"],2]%>%paste0("MAPK:",.)
# 9 hits
map3k <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                             any_vars(str_detect(., fixed("MAP3K",ignore_case=TRUE))))%>%filter(Protein_class=="0")%>%pull(simp.id)
mapping[which(mapping$simp.id %in% map3k),"Protein_class"]<-"MAPK:SlMAPKKK"
## MAPKK:one hit
mpkk <- read.delim("../SFCH/Data/pti/SlMPKK.txt")
rownames(mpkk) <- mpkk$simp.id
mapping[which(mapping$simp.id %in% mpkk$simp.id),"Protein_class"] <- mpkk[mapping[which(mapping$simp.id %in% mpkk$simp.id),"simp.id"],2]%>%paste0("MAPK:",.)
## MAPK:3 hits
mapk_ref <-read.delim("../SFCH/Data/pti/SlMAPK.txt") 
rownames(mapk_ref) <- mapk_ref$simp.id
mapping[which(mapping$simp.id %in% mapk_ref$simp.id),"Protein_class"] <- mapk_ref[mapping[which(mapping$simp.id %in% mapk_ref$simp.id),"simp.id"],2]%>%paste0("MAPK:",.)

##### stress.biotic.signalling #################################
## SlTFT(tomato 14-3-3):1 hit
tft_ref <- read.delim("../SFCH/Data/pti/Sl1433.txt")
rownames(tft_ref) <- tft_ref$simp.id
mapping[which(mapping$simp.id %in% tft_ref$simp.id),"Protein_class"] <- tft_ref[mapping[which(mapping$simp.id %in% tft_ref$simp.id),"simp.id"],2]%>%paste0("14-3-3:",.)

## SlTPR (Tetratricopeptide repeat):16 hits=2+14
tpr_ref <- read.delim("../SFCH/Data/pti/SlTPR.txt")# 2 hits
rownames(tpr_ref) <- tpr_ref$simp.id
mapping[which(mapping$simp.id %in% tpr_ref$simp.id),"Protein_class"] <- tpr_ref[mapping[which(mapping$simp.id %in% tpr_ref$simp.id),"simp.id"],2]%>%paste0("TPR:",.)
# 14 hits against AtTPRs
tpr <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(., pattern ="repeat \\(TPR\\)")))%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% tpr$simp.id),"Protein_class"] <- "TPR:SlTPR"

## Mildew resistance Locus O (MLO) : 3 hits+1 hit=4 hits
#The MLO genes encode plant transmembrane proteins which typically span across the plasma membrane seven times and end in the cytoplasm with a C-terminal domain.
mlo_ref <- read.delim("../SFCH/Data/pti/SlMLOs.txt") # 3hits+1 hit=4 hits
rownames(mlo_ref) <- mlo_ref$simp.id
mlo <- mapping%>%filter(simp.id%in% mlo_ref$simp.id)%>%pull(simp.id) # 3 hits
mapping[which(mapping$simp.id %in% mlo),"Protein_class"] <- mlo_ref[mapping[which(mapping$simp.id %in% mlo),"simp.id"],2]%>%paste0("MLO:",.)
MLO <- setdiff(mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                                   any_vars(str_detect(., fixed("MLO",ignore_case=TRUE))))%>%pull(simp.id)%>%sort(),mlo) # 1 hit
mapping[which(mapping$simp.id %in% MLO),"Protein_class"] <- "MLO:SlMLO"

## IQD/SUN proteins(The tomato fruit shape genes SUN,IQ67 domain):12 hits=11+1
sun <- read.delim("../SFCH/Data/pti/SlSUN.txt")
rownames(sun) <- sun$simp.id
mapping[which(mapping$simp.id %in% sun$simp.id),"Protein_class"] <- sun[mapping[which(mapping$simp.id %in% sun$simp.id),"simp.id"],3]%>%paste0("Ca2+:",.)
# Solyc06g053455 also has IQ-domain
mapping[which(mapping$simp.id %in% "Solyc06g053455"),"Protein_class"] <- "Ca2+:SlIQD"
## Pathogenesis-related :3 tomato PR-1 + PR1-17 from chris: 1 hit
pr_ref <- read.delim("../SFCH/Data/pti/SlPR-1.txt")
rownames(pr_ref) <- pr_ref$simp.id
mapping[which(mapping$simp.id %in% pr_ref$simp.id),"Protein_class"] <- pr_ref[mapping[which(mapping$simp.id %in% pr_ref$simp.id),"simp.id"],2]%>%paste0("PR:",.)

## ABC transporter: 18hits+2hits=20 hits
abc_ref <- read.delim("../SFCH/Data/pti/SlABCtransporter.txt")
rownames(abc_ref) <- abc_ref$simp.id
abc_ref$name <- gsub("\\s+", "", abc_ref$name)
abc <- mapping%>%filter(simp.id%in%abc_ref$simp.id)%>%pull(simp.id)# 18 hits
mapping[which(mapping$simp.id %in% abc),"Protein_class"] <- abc_ref[mapping[which(mapping$simp.id %in% abc),"simp.id"],2]%>%paste0("transport:",.)
ABC <- setdiff(mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., fixed("ABC",ignore_case=TRUE))))%>%pull(simp.id),abc)%>%sort()
# only c("Solyc03g114960","Solyc06g074960") could be ABC-transporters but neither contain relevant domains
mapping[which(mapping$simp.id %in% c("Solyc03g114960","Solyc06g074960")),"Protein_class"] <- c("SlABCB","SlABC-2")%>%paste0("transport:",.)
## S (slow) type anionic channels (SLACs):2 hits
## QUICK- ACTIVATING ANION CHANNEL 1 (QUAC1) :no hits
slac_ref <- read.delim("../SFCH/Data/pti/SlSLAC.txt")
rownames(slac_ref) <- slac_ref$simp.id
mapping[which(mapping$simp.id %in% slac_ref$simp.id),"Protein_class"] <- slac_ref[mapping[which(mapping$simp.id %in% slac_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)
## PM H+-ATPases:Solyc03g113405 LHA1;Solyc07g017780 LHA4
mapping[which(mapping$simp.id %in% c("Solyc03g113405","Solyc07g017780")),"Protein_class"] <- c("SlLHA1","SlLHA4")%>%paste0("transport:",.)
## P5-type cation-transporting ATPase (MIA):Solyc01g096830
mapping[which(mapping$simp.id %in% "Solyc01g096830"),"Protein_class"] <- "transport:SlMIA"
## Vacuolar ATPase (V-type ATPase):VHA 4 hits against AtVHAs + "Solyc12g055800"
vha <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., fixed("V-type",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% c(vha$simp.id,"Solyc12g055800")),"Protein_class"] <- "transport:SlVHA"
## Active calcium transporter:3 hits
aca <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., fixed("calcium cation-",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% aca$simp.id),"Protein_class"] <- c("SlACA2","SlACA8","SlECA4")%>%paste0("Ca2+:",.)
## The Ca2+/cation antiporter (SlCaCA) superfamily 
# 16 CaCAs were grouped into four families: six CAXs, four NCLs, five CCXs, and one MHX
# Ca2+/H+ exchanger:CAX : no hits
# Na+/Ca2+ exchanger:NCL:2 hits
# H+/Ca2+ exchangers (CCXs): no hits
# Mg2+/H+ exchangers (MHXs): no hits
cax_ref <- read.delim("../SFCH/Data/pti/SlCAXs.txt")
rownames(cax_ref) <- cax_ref$simp.id
mapping[which(mapping$simp.id %in% cax_ref$simp.id),"Protein_class"] <- cax_ref[mapping[which(mapping$simp.id %in% cax_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)
## Cation/Proton Antiporter (CPA) 
# SlSOS(Na+/H+ antiporter) 2 hits ; NO hits against CHX; RLKs are all covered
sos_ref <- read.delim("../SFCH/Data/pti/Slchx+sos+rlk.txt") %>%filter(str_detect(name, pattern ="SOS"))
rownames(sos_ref) <- sos_ref$simp.id
mapping[which(mapping$simp.id %in% sos_ref$simp.id),"Protein_class"] <- sos_ref[mapping[which(mapping$simp.id %in% sos_ref$simp.id),"simp.id"] ,2]%>%paste0("transport:",.)
# Sodium(Na+) cation(H+) antiporter (NHX) :3 hits but two are SOS1
# potassium cation antiporter (KEA) #1 hit
cpa_ref <- read.delim("../SFCH/Data/pti/SlCPA.txt")
rownames(cpa_ref) <- cpa_ref$simp.id
mapping[which(mapping$simp.id %in% setdiff(cpa_ref$simp.id,sos_ref$simp.id)),"Protein_class"] <- cpa_ref[mapping[which(mapping$simp.id %in% setdiff(cpa_ref$simp.id,sos_ref$simp.id)),"simp.id"],2]%>%paste0("transport:",.)
## proton antiporter (CLC): 5 hits
clc_ref <- read.delim("../SFCH/Data/pti/SlCLC.txt") # 1/3/6/8 Antiporters;2/4/5/7 channels
rownames(clc_ref) <- clc_ref$simp.id
mapping[which(mapping$simp.id %in% clc_ref$simp.id),"Protein_class"] <- clc_ref[mapping[which(mapping$simp.id %in% clc_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)
## metal cation transporter (MRS/MGT): 3 hits
mrs_ref <- read.delim("../SFCH/Data/pti/SlMRS.txt")
rownames(mrs_ref) <- mrs_ref$simp.id
mapping[which(mapping$simp.id %in% mrs_ref$simp.id),"Protein_class"] <- mrs_ref[mapping[which(mapping$simp.id %in% mrs_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)
## potassium(K+) cation(H+) transporter (HAK/KUP/KT):8 hits
kt_ref <- read.delim("../SFCH/Data/pti/SlKT.txt")
rownames(kt_ref) <- kt_ref$simp.id
mapping[which(mapping$simp.id %in% kt_ref$simp.id),"Protein_class"] <- kt_ref[mapping[which(mapping$simp.id %in% kt_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)
# 2 hits against Arabidopsis potassium cation channel(TPK/KCO)/(AKT/SKOR/GORK)
mapping[which(mapping$simp.id %in% c("Solyc10g006010","Solyc12g006850")),"Protein_class"] <- c("SlTPK1","SlAKT1")%>%paste0("transport:",.)

## plasma membrane intrinsic protein (PIP):aquaporins 8 hits
aqua_ref <- read.delim("../SFCH/Data/pti/SlAquaporins.txt")
rownames(aqua_ref)<- aqua_ref$simp.id
mapping[which(mapping$simp.id %in% aqua_ref$simp.id),"Protein_class"] <- aqua_ref[mapping[which(mapping$simp.id %in% aqua_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)

## cyclic nucleotide-binding transporter(CNGCs):4 hits
cngc_ref<- read.delim("../SFCH/Data/pti/SlCNGC.txt")
rownames(cngc_ref)<- cngc_ref$simp.id
mapping[which(mapping$simp.id %in% cngc_ref$simp.id),"Protein_class"] <- cngc_ref[mapping[which(mapping$simp.id %in% cngc_ref$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)
## Glutamate receptor like receptors(GLRs):1 hit
glr_ref<- read.delim("../SFCH/Data/pti/SlGLR.txt")
rownames(glr_ref) <- glr_ref$simp.id
mapping[which(mapping$simp.id %in% glr_ref$simp.id),"Protein_class"] <- glr_ref[mapping[which(mapping$simp.id %in% glr_ref$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)
## Mechanosensitive channel
# Mechanosensitive channel of Small conductance Like(MSL)
msl_ref <-  read.delim("../SFCH/Data/pti/SlMSLs.txt") # 2 hits against msl
rownames(msl_ref) <- msl_ref$simp.id
mapping[which(mapping$simp.id %in% msl_ref$simp.id),"Protein_class"] <- "Mechano:SlMSL"
# Mid1-Complementing Activity(MCA): best hit against AtMCA1 AT4G35920
mapping[which(mapping$simp.id %in% "Solyc02g083540"),"Protein_class"] <- "Mechano:SlMCA1"
# Piezo, a channel first identified from a mechanical screen of mouse neuro- blastoma cell line and TPK for Two-Pore K þ channel.
# ITAG3.2 Solyc09g074113 BUT ITAG4 Solyc09g074110 : 1 hit
mapping[which(mapping$simp.id %in% "Solyc09g074113"),"Protein_class"] <- "Mechano:SlPiezo1"
# calcium-permeable channel (OSCA):4 hits
#HYPERPOLARIZATION-ACTIVATED Ca2+ CHANNELs (HACCs) were first found in Solanum lycopersicum (tomato) responding to fungal infection
osca <- read.delim("../SFCH/Data/pti/SlOSCA.txt")
rownames(osca) <- osca$simp.id
mapping[which(mapping$simp.id %in% osca$simp.id),"Protein_class"] <- osca[mapping[which(mapping$simp.id %in% osca$simp.id),"simp.id"],2]%>%paste0("Ca2+:",.)
## Two-pore channel(TPCs) (best hit against AtTPC1:AT4G03560)
mapping[which(mapping$simp.id %in% "Solyc07g053970"),"Protein_class"] <- "Ca2+:SlTPC1"
## left calmodulin-binding family protein:15 hits
cam <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(., fixed("calmodulin",ignore_case=TRUE))))%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% cam$simp.id),"Protein_class"] <- "Ca2+:calmodulin"
table(mapping$Protein_class)# 2119-1732=387

# ABA ##################################################################################################################
## SnRK2 kinases - ABA-activated protein kinases
# Solyc08g077780.3.1 AtOST1 homologue
snrk_ref <- read.delim("../SFCH/Data/pti/SlSnRKs.txt") # 8 hits including 4 hits are CIPKs
rownames(snrk_ref) <- snrk_ref$simp.id
snrk <- mapping[which(mapping$simp.id %in% snrk_ref$simp.id),]%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% snrk$simp.id),"Protein_class"] <- snrk_ref[mapping[which(mapping$simp.id %in% snrk$simp.id),"simp.id"],2]%>%paste0("ABA:",.)
# 3 hits
snrk <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                            any_vars(str_detect(., fixed("SNF1",ignore_case=TRUE))))%>%filter(Protein_class=="0")%>%pull(simp.id) 
mapping[which(mapping$simp.id %in% snrk),"Protein_class"] <- "ABA:SlSnRK"

# BR (Brassinosteroid) ######################################################################################################
# BSU: nuclear-localized serine–threonine protein phosphatase with an N-terminal Kelch-repeat domain
bsu <-  mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),any_vars(str_detect(., pattern="BSU")))
mapping[which(mapping$simp.id %in% bsu$simp.id),"Protein_class"] <- paste0("PPase/BR:Sl",mapping[which(mapping$simp.id %in% bsu$simp.id),"arabi.symbol"])
# GSK3-type protein kinase/shaggy-like kinase (SK):BIN2(BR-insensitive 2):4 hits=3+1
sk_ref <-read.delim("../SFCH/Data/pti/SlGSK.txt")
rownames(sk_ref) <- sk_ref$simp.id
mapping[which(mapping$simp.id %in% sk_ref$simp.id),"Protein_class"] <- sk_ref[mapping[which(mapping$simp.id %in% sk_ref$simp.id),"simp.id"] ,2]%>%paste0("BR:",.)
# Solyc06g082640 has 87% seq similarity with SlBIN2
mapping[which(mapping$simp.id %in% "Solyc06g082640"),"Protein_class"] <- "BR:SlSK"
# BSK (AtRLCKXII): 5 hits
bsk_ref <-read.delim("../SFCH/Data/pti/Slbsk.txt")# here it didn't include SlFir1 but replace RLCKXII into SlBSK
rownames(bsk_ref) <- bsk_ref$simp.id
mapping[which(mapping$simp.id %in% bsk_ref$simp.id),"Protein_class"] <- bsk_ref[mapping[which(mapping$simp.id %in% bsk_ref$simp.id),"simp.id"],2]%>%paste0("BR:",.)
## BR biosynthesis:c("Solyc01g111830","Solyc05g008680","Solyc06g035870")
mapping[which(mapping$simp.id %in%c("Solyc01g111830","Solyc05g008680","Solyc06g035870")),"Protein_class"] <- c("BR:C-24","BR:C-24","BR:SlMSBP")

# ethylene #############################################################################################################
# ACC synthase (ACS): 1 hit
# ethylene receptor ETR: 3 hits
# CTR:2 hits
# EIN2:1 hit
# TF:EIL: 1 hit
# ERF: no hits
eth_ref <- read.delim("../SFCH/Data/pti/Sleth.txt") # 8 hits 
rownames(eth_ref) <- eth_ref$simp.id
mapping[which(mapping$simp.id %in% eth_ref$simp.id),"Protein_class"] <- eth_ref[mapping[which(mapping$simp.id %in% eth_ref$simp.id),"simp.id"],2]%>%paste0("ET:",.)

# auxin ################################################################################################################
## TIR1/AFB proteins (Transport Inhibitor Response 1/Auxin Signaling F-Box proteins): 3 hits
tirafb <- c(hmm_all$X3[grep("f-box", hmm_all$X1, ignore.case = TRUE)],
            mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                                any_vars(str_detect(., pattern = "TIR1/AFB")))%>%pull(simp.id))%>%unique()%>% sort()
mapping[which(mapping$simp.id %in% tirafb),"Protein_class"] <- "Auxin:TIR1/AFB"
# Aux/IAA proteins:no hits
# ARF: one hit (in TF)
# small auxin up RNA(saur):no hits
## auxin transporter
# auxin transporter (PIN): 6 hits & auxin transporter (AUX/LAX):1 hit
pin_ref <- read.delim("../SFCH/Data/pti/Slpin.txt")
rownames(pin_ref) <- pin_ref$simp.id
mapping[which(mapping$simp.id %in% pin_ref$simp.id),"Protein_class"] <- pin_ref[mapping[which(mapping$simp.id %in% pin_ref$simp.id),"simp.id"],2]%>%paste0("Auxin:",.)
## TOPLESS (TPL) and TOPLESS-Related (TPR) proteins:
tpl <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(.,fixed("TOPLESS",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% tpl),"Protein_class"] <- "Auxin:SlTPL"

# gibberellins #########################################################################################################
ga <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                          any_vars(str_detect(.,fixed("gibberellin",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% ga$simp.id),"Protein_class"] <- "GA:SlGA20OX"
# jasmonate ############################################################################################################
## JA.biosynthesis
ja <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                          any_vars(str_detect(.,fixed("action.jasmonic",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% ja$simp.id),"Protein_class"] <- "JA:13-LOX"
## JAZ (JASMONATE ZIM-DOMAIN) proteins: no hits
## NINJA (Novel Interactor of JAZ): no hits
# SA #######################################################################################################
## NPR1 (NON-EXPRESSER OF PR GENES1) and its paralogs NPR3 and NPR4: no hits
# TF(transcription factors) ################################################################################################
# Calmodulin Binding Transcription Activator (CAMTA): 5 hits
camta_ref <- read.delim("../SFCH/Data/pti/SlCAMTA.txt")
rownames(camta_ref) <- camta_ref$simp.id
mapping[which(mapping$simp.id %in% camta_ref$simp.id),"Protein_class"] <- camta_ref[mapping[which(mapping$simp.id %in% camta_ref$simp.id),"simp.id"],2]%>%paste0("TF:",.)

# TF:BES1/BZR1 (BRI1-EMS-SUPPRESSOR1/BRASSINAZOLE-RESISTANT1): 1 hit
bes_ref <- read.delim("../SFCH/Data/pti/SlBES1.txt")
rownames(bes_ref)<- bes_ref$simp.id
mapping[which(mapping$simp.id %in% bes_ref$simp.id),"Protein_class"] <- bes_ref[mapping[which(mapping$simp.id %in% bes_ref$simp.id),"simp.id"] ,2]%>%paste0("TF:",.)

# SlWRKY:2 hits
wrky_ref <- read.delim("../SFCH/Data/pti/SlWRKY.txt")
rownames(wrky_ref) <- wrky_ref$simp.id
mapping[which(mapping$simp.id %in% wrky_ref$simp.id),"Protein_class"] <- wrky_ref[mapping[which(mapping$simp.id %in% wrky_ref$simp.id),"simp.id"],2]%>%paste0("TF:",.)

# SlbHLH:3 hits
bhlh_ref <-read.delim("../SFCH/Data/pti/SlbHLH.txt")
rownames(bhlh_ref) <- bhlh_ref$simp.id
mapping[which(mapping$simp.id %in% bhlh_ref$simp.id),"Protein_class"] <- bhlh_ref[mapping[which(mapping$simp.id %in% bhlh_ref$simp.id),"simp.id"],2]%>%paste0("TF:",.)

# other TFs
TF_ref <- read.delim("../SFCH/Data/pti/chris/Sly_TF_list.txt")#%>%filter(Family=="bHLH")
rownames(TF_ref) <- TF_ref$TF_ID
TF_ref$Family <- paste0("TF:",TF_ref$Family)
#table(TF_ref$Family)
tf <- mapping[which(mapping$simp.id %in% TF_ref$TF_ID),]%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% tf$simp.id),"Protein_class"] <- TF_ref[mapping[which(mapping$simp.id %in% tf$simp.id),"simp.id"],3]

# Solyc01g014480 TF/ET:EIL
mapping[which(mapping$simp.id %in% "Solyc01g014480"),"Protein_class"] <- mapping[which(mapping$simp.id %in% "Solyc01g014480"),"Protein_class"]%>%paste0("TF/",.)

## bin 18.4 Protein modification.phosphorylation: 18 protein kinases (PK)
gg <- mapping%>%filter(str_detect(`MapmanX4.2_Bin name`,"Protein modification.phosphorylation"),Protein_class=="0")
# Extract info inside parenthesis
mapping[which(mapping$simp.id %in% gg$simp.id),"Protein_class"]<-gsub("[\\(\\)]","",regmatches(gg$`MapmanX4.2_Bin name`,gregexpr("\\(.*?\\)",gg$`MapmanX4.2_Bin name`)))%>%paste0("PK:",.)

## bin 24 Solute transport: 24.1 primary active transport/24.2 carrier/24.3 channels/24.4 porins: 64
## Major facilitator superfamily (MFS) 
# monosaccharide transporter (AZT)
# glucose transporter
jj <- mapping%>%filter(str_detect(`MapmanX4.2_Bin name`,"Solute transport"),Protein_class=="0")
mapping[which(mapping$simp.id %in% jj$simp.id),"Protein_class"]<-gsub("[\\(\\)]","",regmatches(jj$`MapmanX4.2_Bin name`,gregexpr("\\(.*?\\)",jj$`MapmanX4.2_Bin name`)))%>%paste0("transport:",.)

## SWEET:Sugars will eventually be exported transporter :no hits
## Sucrose transporter(SUT/SUC);Tonoplast monosaccharide transporter(TMT);Vacuolar glucose transporter(VGT);INOSITOL TRANSPORTER(INT); SUGAR FACILITATOR PROTEIN family (SFPs)==EARLY RESPONSE TO DEHYDRATION 6-like (ERD6-like)
swt_ref <- read.delim("../SFCH/Data/pti/SlSWEET.txt")
rownames(swt_ref) <- swt_ref$simp.id
mapping[which(mapping$simp.id %in% swt_ref$simp.id),"Protein_class"] <- swt_ref[mapping[which(mapping$simp.id %in% swt_ref$simp.id),"simp.id"],2]%>%paste0("transport:",.)


###### cell wall integrity (CWI) ################################
## bin 21 Cell wall organisation
# cellulose synthase: cellulose synthase (CesA) and cellulose synthase-like (Csl) 5 hits
CesA_ref <- read.delim("../SFCH/Data/pti/SlCesA.txt") 
rownames(CesA_ref) <- CesA_ref$simp.id
mapping[which(mapping$simp.id %in% CesA_ref$simp.id),"Protein_class"] <- "cellwall:SlCesA/Csl"
# Callose Synthases (CalS): 3 hits+ 2 hits including Solyc07g053980 SlPMR4 (callose synthase)
CalS_ref <- read.delim("../SFCH/Data/pti/SlCalS.txt") 
rownames(CalS_ref) <- CalS_ref$simp.id
mapping[which(mapping$simp.id %in% CalS_ref$simp.id),"Protein_class"] <- CalS_ref[mapping[which(mapping$simp.id %in% CalS_ref$simp.id),"simp.id"],2]%>%paste0("cellwall:",.)
cals <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                            any_vars(str_detect(.,fixed("callose",ignore_case=TRUE))))%>%filter(Protein_class=="0")# 2 hits
mapping[which(mapping$simp.id %in% cals$simp.id),"Protein_class"] <- "cellwall:SlCalS"

# Arabinosyltransferases:4 hits+2 hits
ara_ref <- read.delim("../SFCH/Data/pti/SlArabinosyltransferases.txt")%>%suppressWarnings()
rownames(ara_ref) <- ara_ref$simp.id
mapping[which(mapping$simp.id %in% ara_ref$simp.id),"Protein_class"] <- ara_ref[mapping[which(mapping$simp.id %in% ara_ref$simp.id),"simp.id"],2]%>%paste0("cellwall:",.)
mur <- mapping%>%filter_at(vars(`MapmanX4.2_Bin name`,MapmanX4.2_DESCRIPTION,arabi.symbol:Ara_Mapman_DESCRIPTION),
                           any_vars(str_detect(.,fixed("mur",ignore_case=TRUE))))%>%filter(Protein_class=="0")
mapping[which(mapping$simp.id %in% mur$simp.id),"Protein_class"] <- "cellwall:SlMUR"

# Best hit against AT1G16860.1 SHOU4L : 4 hits
shou <- mapping%>%filter_at(vars(MapmanX4.2_DESCRIPTION:Ara_Mapman_DESCRIPTION),any_vars(str_detect(.,fixed("AT1G16860",ignore_case=TRUE))))
mapping[which(mapping$simp.id %in% shou$simp.id),"Protein_class"] <- "cellwall:SlSHOU4L"
# SlKOR1(KORRIGAN) Solyc01g102580
mapping[which(mapping$simp.id %in% "Solyc01g102580"),"Protein_class"] <- "cellwall:SlKOR1"
# STELLO-type cellulose synthase CSC-interactive protein
mapping[which(mapping$simp.id %in% c("Solyc01g087600","Solyc04g015740")),"Protein_class"] <- "cellwall:STELLO"
# microtubule-interacting component CC of cellulose synthase complex:Solyc04g078750
mapping[which(mapping$simp.id %in% "Solyc04g078750"),"Protein_class"] <- "cellwall:CC"
# left 19 hits
yy <- mapping%>%filter(str_detect(`MapmanX4.2_Bin name`,"Cell wall organisation"),Protein_class=="0")
yy$Protein_class<- gsub("[\\(\\)]","",regmatches(yy$`MapmanX4.2_Bin name`,gregexpr("\\(.*?\\)",yy$`MapmanX4.2_Bin name`)))
yy1 <- yy[setdiff(yy$simp.id,yy$simp.id[yy$Protein_class=="character0"]),]
mapping[which(mapping$simp.id %in% yy1$simp.id),"Protein_class"] <- gsub("[\\(\\)]","",regmatches(yy1$`MapmanX4.2_Bin name`,gregexpr("\\(.*?\\)",yy1$`MapmanX4.2_Bin name`)))%>%paste0("cellwall:",.)
## manually anotate
yy2 <- yy[yy$simp.id[yy$Protein_class=="character0"],]
mapping[which(mapping$simp.id %in% yy2$simp.id),"Protein_class"] <- c("pectin","cutin","cutin","xyloglucan","cutin","expansin","xylan")%>%paste0("cellwall:",.)
## SlPME:Pectin Methylesterases
pme_ref <- read.delim("../SFCH/Data/pti/SlPME.txt")[[1]]
mapping[which(mapping$simp.id %in% pme_ref),"Protein_class"] <- "cellwall:SlPME"

sort(table(mapping$Protein_class))#2119-1493=626
#saveRDS("../SFCH/Data/annotation/mapping_protein_class.rds",object=mapping)

mapping <- readRDS("../SFCH/Data/annotation/mapping_protein_class.rds")

