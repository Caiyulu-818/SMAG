#kegg&go enrichment

#eggNOG result
1. first download the obo file of GO, using go-basic.obo resolve the obo file
python parse_go_obofile.py -i go-basic.obo \
                           -o go.tb
                         
2. egg.annotation result
python3 parse_eggNOG.py -i /function/mag.emapper.annotations -g go.tb -o /egg_result


#Rscript
#csv from the egg_result
library(clusterProfiler)
library(KEGGREST)
library(plyr)
source("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RbioRXN-1.5.1/R/get.kegg.all.R")
source("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/RbioRXN-1.5.1/R/get.kegg.byId.R")
#keggall<-get.kegg.all()
#save(keggall,file="/Users/caiyulu/Desktop/MAGcode/function/keggAll.Rdata")
load(file = "/function/keggAll.Rdata")
#GOannotation <- read.delim("/Users/caiyulu/Desktop/MAGcode/function/GOannotation.tsv", stringsAsFactors=FALSE)
GOinfo <- read.delim("/function/go.tb", stringsAsFactors=FALSE)
# GO enrichment
#setwd("/function/du72")
gocsv<-list.files("/function/egg_result",pattern = "*tsv$")
for(i in gocsv){
GOannotation <- read.delim(i, stringsAsFactors=FALSE)
#GOannotation<-GOannotation[-1,] 
GOgene_list<-GOannotation$gene
## go
gene_go<-enricher(GOgene_list,
                  TERM2GENE=GOannotation[c(2,1)],
                  TERM2NAME=GOinfo[1:2])
write.csv(as.data.frame(gene_go@result),file = paste(substr(i,1,4),"_result.csv"), sep="\t", row.names=F)}

# KEGG enrichment
setwd("/kegg")
kcsv<-list.files("/kegg",pattern = "*.tsv$")
for(i in kcsv){
KOannotation <- read.delim(i, stringsAsFactors=FALSE)
KOgene_list<-KOannotation$gene
kegg_result<-enricher(KOgene_list,
                      TERM2GENE=KOannotation[c(2,1)],
                      TERM2NAME=KOannotation[c(2,4)])
write.csv(as.data.frame(kegg_result@result),file ="du_kegg_result.csv")}


## KOannotation.tsv
## GOannotation.tsv


#combine with the kegg_level 
library(tidyr)
#setwd("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/sedi_function/")
kr<-read.delim("KOannotation.tsv",stringsAsFactors=FALSE)
kr1<-separate(kr,col = gene,into = c("mag",NA),sep = "k")
sedi_ko<-kr1[,c(1,2)]
colnames(sedi_ko)[2]<-"ko"
kegg<-read.csv("kegg_level_table.csv")
sedi_ko_f<-merge(sedi_ko,kegg,by = "ko",all = TRUE)
    