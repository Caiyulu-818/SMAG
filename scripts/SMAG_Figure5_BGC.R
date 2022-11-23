#BGC
#bgc length
setwd("/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/r1")
gbk<-list.files("/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/r1",pattern = "*gbk$")
length.tab<-NULL

for (i in gbk){
  region<-read_table(i)
  start<-as.numeric(region[13,4])
  end<-as.numeric(region[14,4])
  length.tab<-rbind(length.tab,c(i,start,end))
}
write.csv(length.tab,"length.csv")
min(contig$length)
SRR12347719<-subset(sample_info,sample_info$sample_name=="SRR12347719")
SRR12347719bin.1<-filter(maginfo,mag=="SRR12347719bin.1")
P__<-subset(maginfo,maginfo$phylum=="p__")
#筛选>5kb的bgc
ocean<-read.csv("/Users/caiyulu/Desktop/MAGcode/BGC/ocean_Microbiomics Website.csv")
library(readr)
mag<-read_tsv("/Users/caiyulu/Desktop/MAGcode/BGC/Network_Annotations_Full_x1.tsv")
mag<-data.frame(mag[c(-1:-1918),])
length<-read.csv("/Users/caiyulu/Desktop/MAGcode/BGC/length.csv")
colnames(length)[2]<-"region"
colnames(length)[3]<-"start"
colnames(length)[4]<-"end"
length$length<-length$end-length$start
colnames(mag)[1]<-"region"
length$region<-gsub(".gbk","",length$region)
mag_length<-merge(length,mag,by="region")
write.csv(mag_length,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/r1/mag_bgc_length.csv")
bg_5kb<-filter(mag_length,length >=5000)
bg_50kb<-filter(mag_length,length >=50000)#742
bg_100kb<-filter(mag_length,length >=100000)#34
nrps_100kb<-filter(bg_100kb,BiG.SCAPE.class=="NRPS")
nrps_100kb<-separate(nrps_100kb,Description,into = c("bin",NA),sep="k141")
colnames(maginfo)[4]<-"bin"
nrps_100kb<-merge(nrps_100kb,maginfo,by="bin")
#S379 contig edge

write.csv(nrps_100kb,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/nrps——longest/nrps_100kb.csv")

bg_30kb<-filter(mag_length,length >=30000)#4772
bg_100kb<-filter(mag_length,length >=100000)






mag1<-separate(mag,Description,into = c("bin",NA),sep="k141")
mag1_genome_phylum<-aggregate(mag1$bin,by=list(mag1$BiG.SCAPE.class,mag1$bin),length)
genome_bgc<-aggregate(mag1_genome_phylum$x,by=list(mag1_genome_phylum$Group.2),sum)
genome_bgc_phylum<-merge(genome_bgc,maginfo,by="user_genome",all=TRUE)
genome_bgc_phylum<-genome_bgc_phylum[,c(1,2,20)]
genome_bgc_phylum<-subset(genome_bgc_phylum,genome_bgc_phylum$region_num!="NA")
phylum_num<-data.frame(table(genome_bgc_phylum$phylum))
bgcnum_phylum<-aggregate(genome_bgc_phylum$region_num,by=list(genome_bgc_phylum$phylum),sum)
bgc_num_phylum<-cbind(phylum_num,bgcnum_phylum)
colnames(bgc_num_phylum)[1]<-"phylum"
colnames(bgc_num_phylum)[4]<-"region_num"
bgc_num_phylum<-bgc_num_phylum[,-3]
bgc_num_phylum$bgc_per_genome<-bgc_num_phylum$region_num/bgc_num_phylum$genome_num.1

write.csv(bg_5kb,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/r1/mag_bgc_5kb.csv")

NRPS<-filter(bg_5kb,BiG.SCAPE.class=="NRPS")#10277 49
RiPPs<-filter(bg_5kb,BiG.SCAPE.class=="RiPPs")#9632 69
Terpene<-filter(bg_5kb,BiG.SCAPE.class=="Terpene")#7671  45
PKSI<-filter(bg_5kb,BiG.SCAPE.class=="PKSI")#1790 28phyla
Others<-filter(bg_5kb,BiG.SCAPE.class=="Others")#8599 55
PKS_NRP_Hybrids<-filter(bg_5kb,BiG.SCAPE.class=="PKS-NRP_Hybrids")#1664 23
PKSother<-filter(bg_5kb,BiG.SCAPE.class=="PKSother")#3496 47
saccharides<-filter(bg_5kb,BiG.SCAPE.class=="Saccharides")#40 8

#正确
bg_5kb<-separate(bg_5kb,Description,into = c("bin",NA),sep="k141")
bin_bgc<-aggregate(bg_5kb$BiG.SCAPE.class,by=list(bg_5kb$bin,bg_5kb$BiG.SCAPE.class),length)
colnames(bin_bgc)[1]<-"mag"
colnames(bin_bgc)[2]<-"BiG.SCAPE.class"
colnames(bin_bgc)[3]<-"count"
colnames(maginfo)[4]<-"mag"
mag_bgc<-merge(maginfo,bin_bgc,by="mag",all=TRUE)
mag_bgcinfo<-mag_bgc[,c(1,5:6,18:24,43:44)]
mag_bgcinfo<-mag_bgcinfo[,c(1,5,11,12)]
mag_bgcinfo<-filter(mag_bgcinfo,BiG.SCAPE.class!="NA")
mag_bgcinfo$phylum<-gsub("p__","",mag_bgcinfo$phylum)
genome_num<-data.frame(table(mag_bgcinfo$phylum))
genome_num_200<-subset(genome_num,genome_num$Freq>200)#筛选前15
bgcnum_phylum<-aggregate(mag_bgcinfo$count,by=list(mag_bgcinfo$phylum,mag_bgcinfo$BiG.SCAPE.class),sum)
bgcnum_phylum_filter<-bgcnum_phylum%>%filter(Group.1%in%genome_num_200$Var1)
colnames(genome_num)[2]<-"genome_num"
colnames(bgcnum_phylum)[1]<-"phylum"
colnames(bgc_genome_num_nofilter)[3]<-"genome_num"
bgc_num_nofilter<-aggregate(bgcnum_phylum$BGC_count,by=list(bgcnum_phylum$phylum),sum)
bgc_genome_num_nofilter<-merge(bgc_num_nofilter,genome_num,by="phylum")
bgc_genome_num_nofilter$bgc_per_genome<-bgc_genome_num_nofilter$BGC_count/bgc_genome_num_nofilter$genome_num
bgc_num<-aggregate(bgcnum_phylum_filter$BGC_count,by=list(bgcnum_phylum_filter$phylum),sum)
colnames(bgc_genome_num)[3]<-"genome_count"
bgc_genome_num<-merge(bgc_num,genome_num_200,by="phylum")
bgc_genome_num$bgc_per_genome<-bgc_genome_num$bgc_count/bgc_genome_num$genome_count
bgc_fromphylum<-aggregate(bgcnum_phylum$phylum,by=list(bgcnum_phylum$BGC_class),length)
nofilter_bgc<-aggregate(bgcnum_phylum$BGC_count,by=list(bgcnum_phylum$BGC_class),sum)
bgc_from_phyla<-merge(bgc_fromphylum,nofilter_bgc,by="Group.1")
colnames(nofilter_bgc)[1]<-"Group.1"
colnames(bgc_from_phyla)[2]<-"cross_phylum"
colnames(bgc_from_phyla)[3]<-"BGC_count"
write.csv(bgc_from_phyla,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/r1/bgc_from_phyla_r2.csv")
#查看那个基因组编码最多bgc
genome1<-aggregate(mag_bgcinfo$count,by=list(mag_bgcinfo$mag),sum)
SRR6185369rbin.37_bgcnum<-filter(bin_bgc,mag=="SRR6185369rbin.37") #最长
num1<-filter(bin_bgc,mag=="SRR8558721rbin.5") #最多

log<-filter(genome1,Group.1=="SRR11836025bin.3")
#基因组编码最多bgc:SRR8558721rbin.5
#最长的region 对应的基因组 SRR11836025bin.3k141_651121  3453015bp
longest<-filter(bg_5kb,bin=="SRR6185369rbin.37")
num1<-filter(bin_bgc,mag=="SRR8558721rbin.5")
SRR6185369rbin.37<-filter(maginfo,bin=="SRR6185369rbin.37") #最长
#p__Acidobacteriota c__Thermoanaerobaculia o__UBA5704 f__UBA5704 g__UBA5704 s
SRR8558721rbin.5<-filter(maginfo,bin=="SRR8558721rbin.5")# 最多bgc
#p__Acidobacteriota #c__Thermoanaerobaculia #o__UBA5704 #f__UBA5704 #g__#s__
maginfo<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/mag95done.csv")
maginfo$genome<-gsub(".fa","",maginfo$genome)
maginfo$genome<-gsub(".fa","",maginfo$genome)
colnames(maginfo)[4]<-"mag"
ucgb<-filter(maginfo,order=="o__")
genome_num<-data.frame(table(bg_5kb$bin))
genome_num_info<-maginfo%>%filter(user_genome%in%genome_num$Var1)
genome_num_phylum<-data.frame(table(genome_num_info$phylum))
genome_num_phylum$phylum<-gsub("p__","",genome_num_phylum$phylum)
colnames(genome_num_phylum)[1]<-"phylum"
colnames(genome_num_phylum)[2]<-"genome_num"
genome_num_filter<-subset(genome_num_phylum,genome_num_phylum$genome_num>90)

write.csv(genome_num_filter,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_network/genome_num_filter.csv")
mag_bgc<-aggregate(bg_5kb$bin,by=list(bg_5kb$bin,bg_5kb$`Product Prediction`,bg_5kb$`BiG-SCAPE class`),length)
colnames(mag_bgc)[1]<-"mag"
colnames(mag_bgc)[2]<-"Product.Prediction"
colnames(mag_bgc)[3]<-"BiG.SCAPE.class"
colnames(mag_bgc)[4]<-"count"
maginfo<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/mag95done.csv")

mag16<-subset(mag_g_,mag_g_$X16S<1|mag_g_$X23S<1|mag_g_$X5s<1)
magg<-subset(mag16,mag16$X23S==0&mag16$X16S==0&mag16$X5s==0)
#head(maginfo)
maginfo$genome<-gsub(".fa","",maginfo$genome)
colnames(maginfo)[4]<-"mag"
mag_bgc<-merge(maginfo,mag_bgc,by="mag",all=TRUE)
mag_bgc[is.na (mag_bgc)] <- 0
mag_bgcinfo<-mag_bgc[,c(1,5:6,18:24,43:45)]
mag_bgcinfo$phylum<-gsub("p__","",mag_bgcinfo$phylum)
mag_bgcinfo<-aggregate(mag_bgcinfo$mag,by=list(mag_bgcinfo$mag,mag_bgcinfo$phylum,mag_bgcinfo$Product.Prediction,mag_bgcinfo$BiG.SCAPE.class),length)
mag_bgcinfo<-subset(mag_bgcinfo,mag_bgcinfo$Group.3!=0)
write.csv(mag_bgcinfo_phylum,"/Users/caiyulu/Desktop/MAGcode/BGC/mag_bgcinfo_5kb_phylum.csv")
sample<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/sample_ecoinfo_r8.csv")

#5kb bgc——info
mag_bgcinfo<-read.csv("/Users/caiyulu/Desktop/MAGcode/BGC/mag_bgcinfo_5kb.csv")

genome_phylum<-aggregate(mag_bgcinfo$BGC.class,by=list(mag_bgcinfo$mag,mag_bgcinfo$phylum),length)
num_phylum<-data.frame(table(genome_phylum$Group.2))
num_bgc<-aggregate(genome_phylum$x,by=list(genome_phylum$Group.2),sum)
num_bgc_phylum<-cbind(num_bgc,num_phylum)
colnames(num_bgc_phylum)[4]<-"num_genomes"
colnames(num_bgc_phylum)[2]<-"num_bgc"
num_bgc_phylum<-num_bgc_phylum[,-3]
bgc_per_genome<-num_bgc_phylum$num_bgc/num_bgc_phylum$num_genomes
num_bgc_phylum$bgc_per_genome<-bgc_per_genome
write.csv(num_bgc_phylum,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/bgc_per_genome_filter.csv")
mag_bgcinfo_phylum<-aggregate(mag_bgcinfo$phylum,by=list(mag_bgcinfo$phylum,mag_bgcinfo$BGC.class),length)

#mag_bgcinfo_phylum_numkuo<-aggregate(mag_bgcinfo$mag,by=list(mag_bgcinfo$phylum),length)
bgc_100_info<-mag_bgcinfo_phylum%>%filter(Group.1%in%genome_num_filter$phylum)
bgc_num_phylum<-aggregate(bgc_100_info$BGC.count,by=list(bgc_100_info$phylum),sum)
bgc_phy_geno<-merge(bgc_num_phylum,genome_num_filter,by="phylum")
bgc_phy_geno$bgc_per_genome<-bgc_phy_geno$BGC.count/bgc_phy_geno$genome_num

write.csv(bgc_phy_geno,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/bgc_phy_geno.csv")
bgc_group_by<-aggregate(bgc_100_info$BGC.count,by=list(bgc_100_info$phylum),sum)
m<-sum(mag_bgcinfo_phylum1$x)
colnames(bgc_num_phylum)[1]<-"phylum"
colnames(bgc_phy_geno)[2]<-"BGC.count"
colnames(bgc_100_info)[1]<-"phylum"
colnames(bgc_100_info)[2]<-"BGC.class"
colnames(bgc_100_info)[3]<-"BGC.count"

phylum <- phylum[order(phylum$BGC_count),]
m<-(mag_bgcinfo_phylum$BGC_count)
relative_abundance<-mag_bgcinfo_phylum$BGC_count/m
mag_bgcinfo_phylum$relative_abundance<-relative_abundance
mag_bgcinfo_phylum <- mag_bgcinfo_phylum[order(mag_bgcinfo_phylum$relative_abundance),]
mag_bgcinfo_phylum$relative_abundance<-factor(mag_bgcinfo_phylum$phylum,levels = mag_bgcinfo_phylum$relative_abundance,ordered = T)
library(ggplot2)
bgcinfo<-data.frame(t(mag_bgcinfo_phylum))
#genome plot
library(forcasts)
p<-genome_num_200%>%mutate(phylum = fct_reorder(phylum, Freq))%>%
  ggplot(aes(phylum,Freq))+
  geom_bar(stat="identity",position = "stack", width=0.9,size=0.25,color='grey')+
  geom_text(aes(label=Freq),size=2.9,position = position_stack(vjust=1),color="black")+
  coord_flip()+
  labs(x="",y="",title="Numbers of MAGs")+
  theme(plot.title = element_text(hjust = 0.5),size=2) +
  scale_y_continuous(breaks = NULL)+ #隐藏yaxis
  theme_bw()+
  theme(aspect.ratio = 1.8,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank (),plot.margin = unit(rep(1,8),'lines'))
p
position = position_stack(vjust=0.98)
x=3, y=30
panel.border = element_line(colour = "black")
p<-p+scale_x_discrete(breaks=NULL)
#p<-p+theme(plot.margin=unit(rep(3,5),'cm'))

ggsave("/Users/caiyulu/Desktop/MAGcode/BGC/pshot.pdf",height = 4, width =8)

bgcnum_phylum_filter$BGC_class<- factor(bgcnum_phylum_filter$BGC_class,levels=c("Saccharides","PKS-NRP_Hybrids","PKSI","Terpene","PKSother","NRPS","Others","RiPPs"),ordered=T)
bgcnum_phylum_filter$phylum<-factor(bgcnum_phylum_filter$phylum,levels=c( "Thermoproteota","Nitrospirota","Desulfobacterota_B","Firmicutes","Cyanobacteria","Methylomirabilota","Planctomycetota","Verrucomicrobiota","Gemmatimonadota","Chloroflexota","Bacteroidota","Myxococcota","Acidobacteriota","Actinobacteriota","Proteobacteria"),ordered = T)    
                                                          
c( "Firmicutes(146)","Cyanobacteria(239)","Thermoproteota(108)","Desulfobacterota_B(179)","Nitrospirota(184)","Methylomirabilota(240)","Planctomycetota(380)","Verrucomicrobiota(443)","Gemmatimonadota(578)","Myxococcota(982)","Chloroflexota(573)","Bacteroidota(701)","Acidobacteriota(3,023)","Actinobacteriota(3,206)","Proteobacteria(3.808)")    
#bgc_100_info%>%mutate(phylum = fct_reorder(phylum, BGC.count))%>%
library(tidyverse)
p1<-ggplot(data=bgcnum_phylum_filter,aes(phylum,BGC_count,fill=BGC_class))+
  geom_bar(stat="identity",position="fill", width=0.9,size=0.25)+
  coord_flip()+
  xlab("phylum")+
  ylab("BGC_count(%)")+
  scale_fill_manual(values=c("#9ec4b6","#FADE9A","#059f86","#7ca72f","#455D1C","#2a7465","#7b838b","#38415b"),name="BGC")+
  theme_bw()+
  theme(aspect.ratio=1.5,legend.position = "right",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_text(size=9),legend.text=element_text(size=8,colour = "black"))

p1
p1<- p1 +scale_x_discrete(breaks=c( "Thermoproteota","Nitrospirota","Desulfobacterota_B","Firmicutes","Cyanobacteria","Methylomirabilota","Planctomycetota","Verrucomicrobiota","Gemmatimonadota","Chloroflexota","Bacteroidota","Myxococcota","Acidobacteriota","Actinobacteriota","Proteobacteria"),
                          labels=c( "Thermoproteota(244)","Nitrospirota(427)","Desulfobacterota_B(418)","Firmicutes(455)","Cyanobacteria(667)","Methylomirabilota(632)","Planctomycetota(864)","Verrucomicrobiota(1,060)","Gemmatimonadota(2,046)","Chloroflexota(1,457)","Bacteroidota(1,862)","Myxococcota(2,633)","Acidobacteriota(9,116)","Actinobacteriota(9,575)","Proteobacteria(9,823)"))

#p1<-p1+theme(plot.margin=unit(rep(3,5),'cm'))

ggsave("/Users/caiyulu/Desktop/MAGcode/BGC/figure5b1.pdf",height = 4, width = 8)
plot_grid(p,p1,align="v")
library(cowplot)
class<-table(mag_bgcinfo_phylum[2])
class
bin<-data.frame(mag[,5])
colnames(bin)[1]<-"mag"
library(tidyr)
bin<-separate(bin,mag,into = c("bin",NA),sep = '.fa')
write.table(bin,"/Users/caiyulu/Desktop/MAGcode/MAG/drep/undomag.txt")
write.csv(mag_bgcinfo_phylum,"/Users/caiyulu/Desktop/MAGcode/BGC/mag_bgcinfo_phylum.csv")


#bgc_usgb
mag_bgc<-read.csv("/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/mag_bgcinfo_5kb.csv")
usgb<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/mag95info.csv")
usgb$genome<-gsub(".fa","",usgb$genome)
colnames(usgb)[2]<-"mag"
mag_bgc1<-usgb%>%filter(mag%in%mag_bgc$mag)
mag_bgc1<-select(mag_bgc1,c(mag,info))
mag_bgc<-merge(mag_bgc,mag_bgc1,by = "mag",all = TRUE)
write.csv(mag_bgc,"/Users/caiyulu/Desktop/MAGcode/BGC/bgc_annotation/mag_bgc_usgbinfo.csv")
usgbinfo<-select(mag_bgc,c(BGC.class,info))
usgbinfo<-aggregate(usgbinfo$info,by=list(usgbinfo$BGC.class,usgbinfo$info),length)
usgbinfo$info<-gsub("unkSGBs","uSGBs",usgbinfo$info)

colnames(usgbinfo)[1]<-"BGC"
colnames(usgbinfo)[2]<-"info"
colnames(usgbinfo)[3]<-"value"
library(tidyverse)
# 查看颜色
library(RColorBrewer)
library(ggsci)
        # Load

p2<-ggplot(usgbinfo,aes(value,BGC,fill=info))+
  geom_bar(stat="identity",position="fill", width=0.7,size=0.25)+
  xlab("Percentage of BGC")+
  #scale_color_aaas(color = pal_aaas(alpha = .5)(3)[2])+
  scale_fill_manual(values=c("#009075","#B76262"))+
  theme_bw()+
  theme(aspect.ratio = 2,legend.direction = "horizontal",legend.position = "top",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_text(size=9),legend.text=element_text(size=8,colour = "black"))
        
p2
ggsave("/Users/caiyulu/Desktop/MAGcode/BGC/p2_.pdf",height = 5, width = 6)

install.packages("ggeasy") # 安装包
library(ggeasy)
??ggeasy

sample_info<-read.csv("/Users/caiyulu/Desktop/MAGcode/MAG/drep/sample_ecoinfo_r8.csv")

sample_bin<-select(sample_info,mag_name,sample_name)
write_csv(sample_bin,"/Users/caiyulu/Desktop/MAGcode/MAG/drep/sample_bin.csv")

