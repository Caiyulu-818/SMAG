library(readr)
setwd("/home/ec2-user/data1/LUCY/sediment/sedi_gbk")
gbk<-list.files("/home/ec2-user/data1/LUCY/sediment/sedi_gbk",pattern = "*gbk$")
l.tab<-NULL
for (i in gbk){
  region<-read_table(i)
  name<-region[1,2]
  start<-as.numeric(region[13,4])
  end<-as.numeric(region[14,4])
  l.tab<-rbind(l.tab,c(i,name,start,end))
}
l.tab<-data.frame(l.tab)
colnames(l.tab)[1]<-"region"
colnames(l.tab)[2]<-"mag"
colnames(l.tab)[3]<-"start"
colnames(l.tab)[4]<-"end"
l.tab$region<-as.character(l.tab$region)
l.tab$mag<-as.character(l.tab$mag)
l.tab$start<-as.numeric(l.tab$start)
l.tab$end<-as.numeric(l.tab$end)
l.tab$len<-l.tab$end-l.tab$start
write.csv(l.tab,"bgc_length.csv")