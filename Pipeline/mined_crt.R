#crt
library("Biostrings")
library(tidyr)
files <- list.files("./",pattern=".txt")

list_string_diff = function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE, only.position = TRUE){
    if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ")
    if(ignore.case){
        a = toupper(a)
        b = toupper(b)
    }

    split_seqs = strsplit(c(a, b), split = "")
    only.diff = split_seqs[[1]] != split_seqs[[2]]
    only.diff[
        (split_seqs[[1]] %in% exclude) |
        (split_seqs[[2]] %in% exclude)
    ] = NA

    diff.info = data.frame(which(is.na(only.diff)|only.diff),
                                 split_seqs[[1]][only.diff], split_seqs[[2]][only.diff])
    names(diff.info) = c("position", "seq.a", "seq.b")

    if(!show.excluded) diff.info = na.omit(diff.info)
    if(only.position){
        diff.info$position
    }else diff.info
}

for (file in files){
seqs <- readLines(file)
pos <- grep("CRISPR",seqs)
ratios <- NULL
names <- NULL
spacers <- NULL
for (i in 1:(length(pos)-1)){
name <- seqs[pos[i]]
name <- gsub("\\s+Range:\\s+\\S+\\s+\\S+\\s+\\S+","",name)
name <- gsub(" ","_",name)
target <- seqs[(pos[i]+3):(pos[i+1]-6)]
target <- strsplit(target,split="\t")
repeats <- unlist(lapply(target,`[`,3))
num <- table(repeats)
max <- names(num[num==max(num)])[1]
sumdiff <- 0
sumlength <- nchar(max)*length(repeats)
for (rep in repeats){
diff <- length(list_string_diff(rep,max))
sumdiff <- sumdiff+diff
spacer <- pos[i+1]-pos[i]-9
}
ratio <- sumdiff/sumlength
names <- c(names,name)
ratios <- c(ratios,ratio)
spacers <- c(spacers,spacer)
}
name <- seqs[pos[length(pos)]]
name <- gsub("\\s+Range:\\s+\\S+\\s+\\S+\\s+\\S+","",name)
name <- gsub(" ","_",name)
target <- seqs[(pos[length(pos)]+3):(length(seqs)-9)]
target <- strsplit(target,split="\t")
repeats <- unlist(lapply(target,`[`,3))
num <- table(repeats)
max <- names(num[num==max(table(repeats))])[1]
sumdiff <- 0
sumlength <- nchar(max)*length(repeats)
for (rep in repeats){
diff <- length(list_string_diff(rep,max))
sumdiff <- sumdiff+diff
spacer <- length(seqs)-pos[length(pos)]-8
}
ratio <- sumdiff/sumlength
names <- c(names,name)
ratios <- c(ratios,ratio)
spacers <- c(spacers,spacer)
data <- data.frame(names=names,ratios=ratios,spacers=spacers)
data <- data[(data$ratios<=0.03 & data$spacers>=3),]
spacerdatas <- data.frame(name=NULL,seq=NULL)
for (i in 1:dim(data)[1]){
name <- paste0(gsub("_"," ",data$names[i])," ")
posspacer<-grep(name,seqs)
targetspacer <- seqs[(posspacer+3):(posspacer+2+data$spacers[i])]
targetspacer <- strsplit(targetspacer,split="\t")
spacerseq <- unlist(lapply(targetspacer,`[`,4))
spacername <- paste0(">",gsub(".txt","_",file),data$names[i],"_",c(1:data$spacers[i]))
spacerdata <- data.frame(name=spacername,seq=spacerseq)
spacerdatas <- rbind(spacerdatas,spacerdata)
}
spacerdatas <- unite(spacerdatas,"nameseq",name,seq,sep="\n")
writeLines(spacerdatas$nameseq,gsub(".txt","_spacers.fa",file))
}


#minced
library("Biostrings")
setwd("./MAGnoderup_10k.fasta.split")
files <- list.files("./",pattern=".txt")

list_string_diff = function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE, only.position = TRUE){
    if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ")
    if(ignore.case){
        a = toupper(a)
        b = toupper(b)
    }

    split_seqs = strsplit(c(a, b), split = "")
    only.diff = split_seqs[[1]] != split_seqs[[2]]
    only.diff[
        (split_seqs[[1]] %in% exclude) |
        (split_seqs[[2]] %in% exclude)
    ] = NA

    diff.info = data.frame(which(is.na(only.diff)|only.diff),
                                 split_seqs[[1]][only.diff], split_seqs[[2]][only.diff])
    names(diff.info) = c("position", "seq.a", "seq.b")

    if(!show.excluded) diff.info = na.omit(diff.info)
    if(only.position){
        diff.info$position
    }else diff.info
}

for (file in files){
seqs <- readLines(file)
pos <- grep("CRISPR",seqs)
ratios <- NULL
names <- NULL
spacers <- NULL
for (i in 1:(length(pos)-1)){
name <- seqs[pos[i]]
name <- gsub("\\s+Range:\\s+\\S+\\s+\\S+\\s+\\S+","",name)
name <- gsub(" ","_",name)
target <- seqs[(pos[i]+3):(pos[i+1]-4)]
target <- strsplit(target,split="\t")
repeats <- unlist(lapply(target,`[`,3))
num <- table(repeats)
max <- names(num[num==max(num)])[1]
sumdiff <- 0
sumlength <- nchar(max)*length(repeats)
for (rep in repeats){
diff <- length(list_string_diff(rep,max))
sumdiff <- sumdiff+diff
spacer <- pos[i+1]-pos[i]-7
}
ratio <- sumdiff/sumlength
names <- c(names,name)
ratios <- c(ratios,ratio)
spacers <- c(spacers,spacer)
}
name <- seqs[pos[length(pos)]]
name <- gsub("\\s+Range:\\s+\\S+\\s+\\S+\\s+\\S+","",name)
name <- gsub(" ","_",name)
target <- seqs[(pos[length(pos)]+3):(length(seqs)-3)]
target <- strsplit(target,split="\t")
repeats <- unlist(lapply(target,`[`,3))
num <- table(repeats)
max <- names(num[num==max(table(repeats))])[1]
sumdiff <- 0
sumlength <- nchar(max)*length(repeats)
for (rep in repeats){
diff <- length(list_string_diff(rep,max))
sumdiff <- sumdiff+diff
spacer <- length(seqs)-pos[length(pos)]-6
}
ratio <- sumdiff/sumlength
names <- c(names,name)
ratios <- c(ratios,ratio)
spacers <- c(spacers,spacer)
data <- data.frame(names=names,ratios=ratios,spacers=spacers)
data <- data[(data$ratios<=0.03 & data$spacers>=3),]
data$names <- paste0(data$names,"_")
seqfilen <- gsub(".txt","_spacers.fa",file)
seqfile <- readDNAStringSet(seqfilen)
extractname <- names(seqfile)[grep(paste(as.character(data$names),collapse = "|"),names(seqfile))]
grepseq <-  extractList(seqfile,extractname)
writeXStringSet(grepseq@unlistData,paste0(gsub(".txt","",file),"_spacers_filter.fa"))
}
}

#spacer
makeblastdb -in MAG_cas10k_minced_spacers_grep_norep.fa -dbtype nucl
blastn -query GSV_2022.fasta -db MAG_cas10k_minced_spacers_grep_norep.fa -out MAG_cas10k_minced_spacers_grep_norep.out -outfmt 6 -num_threads 100 &
awk '$5 <= 1 && $3>=95' MAG_cas10k_spacers_grepall_25.out > MAG_cas10k_spacers_grepall_25_demis.out

#prophage
makeblastdb -in GSV_2022.fasta -dbtype nucl
blastn -query MAG_renamecatall.fa -db GSV_2022.fasta -out prophage.out -outfmt 6 -num_threads 200 
awk '$3 >90 && $4 > 500' prophage.out > prophage_500_90.out
seqkit fx2tab -n -l GSV_2022.fasta > GSV_2022_nl.txt
seqkit fx2tab -n -l MAG_renamecatall.fa > MAG_renamecatall_nl.txt



