#!/bin/sh

#Public data conversion
fastq-dump /rawdata/$sampleid --split-3 --gzip -O /rawdata/$sampleid,(Double-ended data)

#Rawdata to CleanData
java -jar /public/home/bma/application/Trimmomatic-0-2.39/trimmomatic-0.39.jar PE -threads 10 \
/rawdata/$sampleid_1.fastq.gz /rawdata/$sampleid2.fastq.gz \
/02.paired_data/$sampleid_paired_1.clean.fastq.gz \
/02.unpaired_data/$sampleid_unpaired_1.clean.fastq.gz \
/02.paired_data/$sampleid_paired_2.clean.fastq.gz \
/02.unpaired_data/$sampleid_unpaired_2.clean.fastq.gz \
ILLUMINACLIP:/adapter-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

Attention:
(1)the pathway to adapter.fa
(2)the order: p1 unpaired1 p2 unpaired2

#Assembly
megahit -t 10 -m 0.5 --min-contig-len 500 --k-step 10 --k-min 27 \
-1 /02.paired_data/$sampleid_paired_1.clean.fastq.gz \\
-2 /02.paired_data/$sampleid_paired_2.clean.fastq.gz \
-o /assembly/$sampleid

#Binning
metawrap binning  -t 28 -m 1000 -l 1500 \
-o /binning/$sampleid_binning \
 -a /assembly/$sampleid.fa 
 --metabat2 --maxbin2 --concoct \
 /02.paired_data/$sampleid_*.fastq \


#refine
metawrap bin_refinement  \
-o /refine/$sampleid_refine \
 -t 10 \
 -A concoct_bins/ -B metabat2_bins/ -C maxbin2_bins/ \
 -c 50 -x 10 \

 metawrap_50_10_bins is the bins after refining

 #dereplication and clustering

##dereplication 40,039 soil mag
 dRep dereplicate SMAG_drep/ \
  -g /SMAGALL/*.fa \
  --S_algorithm fastANI -comp 50 -con 10 -pa 0.9 -sa 0.95 -nc 0.1 \
  -cm larger --clusterAlg centroid --multiround_primary_clustering --primary_chunksize 5000

##clustering
dRep dereplicate /genomeall_drep/ \
-p 56 -g /genomeall/*.fa \
--S_algorithm fastANI -comp 50 -con 10 -l 50000 -pa 0.9 -sa 0.95 -nc 0.1 
--multiround_primary_clustering --primary_chunksize 12000(when the number of genomes > 20000)

#Taxonomic analyses
gtdbtk classify_wf --cpus 20 \
 --genome_dir /SMAG/ --out_dir /GTDB_SMAG/ -x fa \

#tRNA & rRNA identification

##trna
for bacteria: for i in ${ID[@]};do tRNAscan-SE -B -Q -o /trna/${i}trna.out /bac_mag/${i}.fa;done
for archae: for i in ${ID[@]};do tRNAscan-SE -A -Q -o /trna/${i}trna.out /arc_mag/${i}.fa;done

##rrna

for bacteria: for i in ${ID[@]};do barrnap --quiet -k bac /bac_mag/${i}.fa --outseq /rrna_hit/${i}hit.fa > /rrna_out/${i}rrna.out;done
for archae: for i in ${ID[@]};do barrnap --quiet -k arc /arc_mag/${i}.fa --outseq /rrna_hit/${i}hit.fa > /rrna_out/${i}rrna.out;done


#Phylogenetic analyse

##construct your own config_file (set the parameters)
ex:phylophlan_write_config_file  \
    -d a \
    -o newphylo.cfg \
  #--force_nucleotides \(based on the input)
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 iqtree \
    --tree2 raxml \
    --overwrite \
    --verbose 2>&1

##tree
phylophlan -i  SMAG_drep/ -d phylophlan --diversity high --accurate \
-f /config/newphylo.cfg \
-o /tree/SMAG_drep_tree --nproc 40 --verbose 2>&1 | tee tree.log


#Abundance of MAGs

##rename the contig name in the mag
sed -i 's/>k141/>$samplein.11k141/g' $samplebin.1.fa

##cat the mag(after rename)
cat ./*.fa > catmag.fa

##index mag
bowtie2-build ./catmag.fa /catindex/catmag

##mapping the read to magindex
ID=(sample)
for i in ${ID[@]};do
bowtie2 -p 28 -x /catindex/catmag \
-1 /02.paired_data/${i}_paired_1.clean.fastq.gz  \
-2 /02.paired_data/${i}paired_2.clean.fastq.gz | samtools sort -O bam -@ 28 -o - > /bam/${i}.bam;done

##filter the read-percent-identity and read-aligned-percent
ID=(sample)
for j in ${ID[@]}; do coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 \
-b /bam/${i}.bam -o /bam/${j}_filter.bam -t 32;done

##Depth of coverage
ID=(sample)
for j in ${ID[@]}; do coverm contig --trim-max 90 --trim-min 10 --methods trimmed_mean \
 --bam-files /bam/${j}_filter.bam > /csv/${j}_filter.bam_coverage.csv;done


##weight by the contig length in the mag

###rscript
#Rscript
library(tidyr)
library(data.table)
#read the mapping coverage csv
csv<-list.files("~/data1/LUCY/mapping/catdone",pattern="*coverage.csv$") 
#read the contig length
contig_length<-fread("/media/contig.txt")
colnames(contig_length)[1]<-"mag"
#colnames(num_contig)[2]<-"mag"
for (i in csv){
coverage<-fread(i,sep="\t",header=F)
coverage1<-coverage[-1,]
colnames(coverage1)[1]<-"mag"
catcoverage<-merge(coverage1,contig_length,by="mag",all=TRUE)
catcoverage<-separate(catcoverage,mag,into=c("mag",NA),sep="k")
catcoverage$V2<-as.numeric(catcoverage$V2)
catcoverage$length<-as.numeric(catcoverage$length)
catcoverage$cat_abundance<-(catcoverage$V2)*(catcoverage$length)
catcoverage1<-aggregate(cbind(catcoverage$length,catcoverage$cat_abundance),by=list(catcoverage$mag),sum)
#weight by the contig length
abundance<-catcoverage1$V2/catcoverage1$V1
catcoverage1$ave_abundance<-abundance
write.csv(catcoverage1,file=paste("/media/length/cat_",i,sep=""))


#MAG function annotation

##Protein clustering
mmseqs easy-linclust /catmag/mag.faa clusterRes tmp --min-seq-id 0.95 -c 0.9 --threads 96 --cluster-mode 2 --cov-mode 1
~clusterRes_rep_seq.fasta is representative amino acids

##Functiona annotation
python3 /eggnog-mapper-2.1.7/emapper.py  \
-i /catmag/mag_repre.fasta \
--report_orthologs --output /function/mag_repre -m diamond --cpu 36 


#BGC identification
antismash -c 96 --cb-general --cb-knownclusters --cb-subclusters --cc-mibig --genefinding-tool prodigal --allow-long-headers /catmag/mag_repre.fasta

##BGC classfication

python /BiG-SCAPE/bigscape.py \
-i /all_gbk/mag_gbk \
-o /BGC/BGC_class \
-c 80 --cutoffs 0.3 --clan_cutoff 0.2 0.8 --mibig 

