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

#SNV catalog construction

##01.snv-calling & shared-snp.list

SPECIES=112.sp001002135

GENOMES_DIR=/03.snv/101-150/genome
PAIRS_IN=/03.snv/101-150/$SPECIES"_CATALOG"/input.tsv
OUTPUT_DIR=/03.snv/101-150/$SPECIES"_CATALOG"
OUTPUT_DIR2=/03.snv/101-150/$SPECIES"_CATALOG"/snps


while IFS= read -r pair
do
            genome1=$(echo $pair | cut -d' ' -f1)
            genome2=$(echo $pair | cut -d' ' -f2)

                GB_PATH_1=$GENOMES_DIR/$genome1.fa
                        GB_PATH_2=$GENOMES_DIR/$genome2.fa

                            nucmer $GENOMES_DIR/$genome1.fa $GENOMES_DIR/$genome2.fa --prefix $OUTPUT_DIR/$genome1-$genome2
                                delta-filter -q -r $OUTPUT_DIR/$genome1-$genome2.delta > $OUTPUT_DIR/$genome1-$genome2.filter.delta
                                    show-coords $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.coords
                                        show-snps -Clr $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR2/$genome1-$genome2.snps
                                            show-diff $OUTPUT_DIR/$genome1-$genome2.filter.delta > $OUTPUT_DIR/$genome1-$genome2.diff

                                                echo "done whole genome alignment for $genome1 and $genome2"
                                        done < $PAIRS_IN
#get share snp List files
ls $OUTPUT_DIR2 | grep snps | cut -d'.' -f1 | xargs -I[] bash -c 'sed "1,5d" '$OUTPUT_DIR2'/[].snps | awk "$0"' '$2 != "." && $3 != "." {printf "%s\t%s||%s||%s||%s\n", "[]", $14, $1, $2, $3}' | cut -f2 | LC_ALL=C sort -k2,2 -S20G --parallel=4 | uniq -c | awk '$1 > 1 {print $2}' > $SPECIES.snps.list

##02.Snv catalog
###method 1
python3 generate_catalog_revised.py --shared 112.sp001002135.snps.list --in-list <(ls /03.snv/101-150/112.sp001002135_CATALOG/snps/*.snps) --out 112.sp001002135.catalog.tsv

### method 2
SNV_CATALOG.R

#Pangenome construction
## genomes selestion
completeness>=80%
contamination<=5%
HQ genomes >=10 
## fasta to gff (prokka v1.14.5)
prokka ./mag.fa ./mag.gff

## create a pangenome with a minimum amino acid identity at 90% (‘-i 90’) and a core gene defined at 90% presence (‘-cd 90’)
roary -e --mafft -i 90 -cd 90 -f  output_dir  *.gff


#Crispr-cas system information

##crispr array
###contig filter
conda activate pangenome
for i in $(cat mag.list); do perl 0.casfind-1.pl /${i}.fa /2.5eur/${i}.3kb.fa /cas/${i}.temp;done

###crispr array and protein prediction
for i in $(cat mag.list);do pilercr -in /cas/${i}.temp -out /cas/${i}.spacer;done
for i in $(cat mag.list);do prodigal -a /cas/$id.pep -i /cas/${i}.temp  -p single -f gff -o /cas/${i}.temp/$id.temp.gff;done

### protein filter 10kb region of crispr array
for i in $(cat mag.list);do perl 0.casfind-2.pl /cas/${i}.spacer /${i}.fa;done
for i in $(cat mag.list);do perl 0.casfind-3.pl /cas/${i}.temp.gff /${i}.fa;done
for i in $(cat mag.list);do bedtools intersect -wo -a /cas/${i}.temp.gff.loc -b /cas/${i}.spacer.loc > /cas/${i}.bed;done
for i in $(cat mag.list);do perl 0.casfind-4.pl /cas/$id.temp.bed /cas/$id.pep /cas/$id.temp.pep.fasta /cas/$id.pep.cas.fasta /$id.fa;done

##cas protein identification and alignment

###du-replicated
cd-hit -i /casdb/Casdb.fa -o /casdb/NR/Casdb.0.95.fa -c 0.95 -aS 0.9 -n 5 

###filter protein length (200-1000aa)
open(AA,"Casdb.fa");
open(BB,">Casdb_filter.fa");
while($line=<AA>)
{  chomp $line;
   $line2=<AA>;  chomp $line2;
   $ll=$line2;
   $len=length($ll);
   if ($len>=200 & $len<=1000){
   print BB $line,"\_L",length($line2),"\n",$line2,"\n";}
}
close AA; close BB;

###(cas novelty) align to NR database
/diamond blastp -k 1 -e 1e-10  --query-cover 40 -d /uniref_100.dmnd -q /casdb/NR/Casdb.0.95.fa -o /casdb/NR/CAS


#Connecting MAGs to viruses
#Linkages from Crispr-cas system information

##Remove MAGs containing fewer than four CRISPR-associated proteins.
##https://github.com/sandialabs/CasCollect/tree/master/ref
##Download all the hmm files of cas proteins
prodigal -i mag10k.fa -a mag10k.faa -d mag10k.ffn -p meta -f gff > 
hmmsearch -Z 1 --noali --cpu 10 --tblout Cas3_0_ID.out Cmr3_1_IIIB.hmm MAG_protein.faa

##filter CRISPR-associated proteins >=4
###R
library(data.table)
casout <- read.csv("cas.out",comment.char="#",sep="",header=F)
casout <- casout[casout$V5<1e-5,]
proteins <- unique(casout$V1)
bins <- gsub("_\\S+","",proteins)
count <- data.frame(table(bins))
count <- count[count$Freq>=4,]
data <- data.frame(bins,proteins)
data <- data[data$bins%in%count$bins,]
seqs <- unique(gsub("_\\d+$","",data$proteins))
writeLines(seqs,"cas_bins_name.txt")

##sequence grep
/public/ylwang/seqkit grep -f cas_bins_name.txt MAG_renamecatall.fa > MAG_cas.fa

##Contigs longer than 10 kb
seqkit seq -m 10000 MAG_cas.fa > MAG_cas10k.fa

##minced(0.4.2)
/application/minced-0.4.2/minced -gffFull -spacers MAGnoderup_10k.part_001.fasta MAGnoderup_10k.part_001.txt MAGnoderup_10k.part_001.gff

##CRT crispr
java -cp /application/CRT1.2-CLI.jar crt MAGnoderup_10k.part_030.fasta MAGnoderup_10k.part_030.out
## Prediction
virsorter2.sif run -w test -i mag.fa -j 4 all

## Quality Control
checkv end_to_end input_file.fna output_directory -t 28


##prophage
makeblastdb -in GSV_2022.fasta -dbtype nucl
nohup blastn -query MAG_renamecatall.fa -db GSV_2022.fasta -out prophage.out -outfmt 6 -num_threads 200 &
awk '$3 >90 && $4 > 500' prophage.out > prophage_500_90.out
seqkit fx2tab -n -l GSV_2022.fasta > GSV_2022_nl.txt
seqkit fx2tab -n -l MAG_renamecatall.fa > MAG_renamecatall_nl.txt



 



