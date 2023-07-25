
# Soil microbial dark matter explored from genome-resolved metagenomics

This directory contains scripts related to the manuscript "Soil microbial dark matter explored from genome-resolved metagenomics". 

Before running, you must ensure that all required softwares and databases are installed successfully. 

The installation method refer to the manual of each software. The name, version and availablity of the software are as follows:  
|Software|Availability|
|:-----|:---------|
|fastq-dump (v2.9.6)|https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk|
|Trimmomatic (v2.39)|https://github.com/usadellab/Trimmomatic|
|MEGAHIT (v1.2.9)|https://github.com/voutcn/megahit|
|Metabat (v2.12.1)|https://bitbucket.org/berkeleylab/metabat|
|MaxBin (v2.2.6)|https://sourceforge.net/projects/maxbin/|
|CONCOCT|https://github.com/BinPro/CONCOCT|
|metaWRAP(v1.2.1)|https://github.com/bxlab/metaWRAP|
|CheckM (v1.0.11)|https://github.com/Ecogenomics/CheckM|
|Barrnap (v.0.9)|https://github.com/tseemann/barrnap|
|tRNAscan-SE (v.2.0.9)|https://github.com/UCSC-LoweLab/tRNAscan-SE|
|dRep (v2.2.4)|https://github.com/MrOlm/drep|
|Mash (v2.3)|https://github.com/marbl/mash|
|Mummer (v4.0.0)|https://github.com/gmarcais/mummer|
|FastANI|https://github.com/ParBLiSS/FastANI|
|GTDB-TK (v.1.6.0)|https://github.com/Ecogenomics/GtdbTk|
|PhyloPhlAn (v3.0.60)|https://github.com/biobakery/phylophlan|
|Diamond (v0.9.14.115)|https://github.com/bbuchfink/diamond|
|mafft (version v7.310)|https://github.com/The-Bioinformatics-Group/Albiorix/wiki/mafft|
|trimal (version 1.4rev15)|https://github.com/scapella/trimal|
|RAxML (version 8.1.12)|https://github.com/stamatak/standard-RAxML|
|FastTree (version 2.1.10)|https://github.com/PavelTorgashov/FastTree|
|HMM|https://github.com/guyz/HMM|
|Prodigal (v2.6.3)|https://github.com/hyattpd/prodigal/wiki|
|BWA (v0.7.17)|https://github.com/lh3/bwa|
|Samtools (v1.10)|https://github.com/samtools/|
|IQ-TREE (v1.6.6)|http://www.cibiv.at/software/iqtree|
|antiSMASH(v6.1)|https://github.com/antismash/antismash|
|BiG-SCAPE|https://github.com/medema-group/BiG-SCAPE|
|Clustal (v2.1)|http://www.clustal.org/|
|BlastKOALA (v.2.21)|https://www.kegg.jp/blastkoala/|
|Roary (v3.12.0)|https://github.com/sanger-pathogens/Roary|
|MMseqs2|https://github.com/soedinglab/MMseqs2|
|eggNOG-mapper (v2.1.6)|https://github.com/eggnogdb/eggnog-mapper|
|PILER-CR|https://github.com/widdowquinn/pilercrpy|
|VirSorter2 (v2.0 alpha)|https://github.com/jiarong/VirSorter2|
|CheckV (v1.0)|https://bitbucket.org/berkeleylab/CheckV|

Note: Make all needed command of software availabled in the "~/bin" directory or in system environment variables.The version is only the version used in the paper and does not have to be the same,  and some softwares are included in other software, so you don't have to install it repeatedly. For example, bwa, bowtie2, Samtools and MEGAHIT are included in metawrap. 


# OVERVIEW OF PIPELINE

The scripts of metagenomic analysis are placed in "[Pipeline](https://github.com/Caiyulu-818/SMAG/blob/main/Pipeline/)" 

#  download the public sra

~/Applications/Aspera\ Connect.app/Contents/Resources/ascp -QT -l 300m -P33001 -i /Users/caiyulu/Applications/Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/SRR695/000/SRR6956730 /Volumes/sra/\
## $sampleid=SRR6956730

#  Public data conversion

##  install
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz
tar zxvf sratoolkit.2.10.8-centos_linux64.tar.gz
cd sratoolkit.2.10.8-centos_linux64/bin
vdb-config -i # position
vi ~/.bashrc
i
export PATH=$PATH:/home/urname/local/app/sratoolkit/bin
ESC, :wq  
source ~/.bashrc # source 
##  use
fastq-dump /rawdata/$sampleid --split-3 --gzip -O /rawdata/$sampleid,(Double-ended data)
##   example output
/rawdata/$sampleid_1.fastq.gz 
/rawdata/$sampleid2.fastq.gz 


# Rawdata to CleanData

##  install
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip
unzip Trimmomatic-0-2.39.zip
cd Trimmomatic-0-2.39
ls
├── adapters
│ ├── NexteraPE-PE.fa
│ ├── TruSeq2-PE.fa
│ ├── TruSeq2-SE.fa
│ ├── TruSeq3-PE-2.fa
│ ├── TruSeq3-PE.fa
│ └── TruSeq3-SE.fa
├── LICENSE
└── trimmomatic-0.39.jar
1 directory, 8 files

##  use
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

# Assembly

##  install
git clone https://github.com/voutcn/megahit.git
cd megahit
make

##  use
megahit -t 10 -m 0.5 --min-contig-len 500 --k-step 10 --k-min 27 \
-1 /02.paired_data/$sampleid_paired_1.clean.fastq.gz \\
-2 /02.paired_data/$sampleid_paired_2.clean.fastq.gz \
-o /assembly/$sampleid
time 8h/sample 96cpu/196GB

# Binning

##  install-conda
conda create -y -n metawrap-env python=2.7
conda activate metawrap-env
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky

##  Unix/Linux only
conda install --only-deps -c ursky metawrap-mg

##   show config
metaWRAP --show-config

##  install-docker in aws
see Dockerfile (metawrapdocker) 

##  use
source ~/miniconda3/bin/activate
conda activate metawrap-env
metawrap binning  -t 28 -m 1000 -l 1500 \
-o /binning/$sampleid_binning \
 -a /assembly/$sampleid.fa 
 --metabat2 --maxbin2 --concoct \
 /02.paired_data/$sampleid_*.fastq 
## attention to the format of the reads(name_*.fastq) 
## -l 1500 minimum contig length to bin (default=1000bp). 
## Note: metaBAT will default to 1500bp minimum


# refine

##  install same to above
metawrap bin_refinement  \
-o /refine/$sampleid_refine \
 -t 10 \
 -A concoct_bins/ -B metabat2_bins/ -C maxbin2_bins/ \
 -c 50 -x 10 \

 metawrap_50_10_bins is the bins after bin_refinement


 # dereplication and clustering

 ##  install
conda config --add channels bioconda
conda install drep
dRep check_dependencies

##  dereplication 40,039 soil mag
conda config --add channels bioconda
conda install drep
dRep check_dependencies
dRep dereplicate SMAG_drep/ \
  -g /SMAGALL/*.fa \
  --S_algorithm fastANI -comp 50 -con 10 -pa 0.9 -sa 0.95 -nc 0.1 \
  -cm larger --clusterAlg centroid --multiround_primary_clustering --primary_chunksize 5000
##  clustering
dRep dereplicate /genomeall_drep/ \
-p 56 -g /genomeall/*.fa \
--S_algorithm fastANI -comp 50 -con 10 -l 50000 -pa 0.9 -sa 0.95 -nc 0.1\
--multiround_primary_clustering --primary_chunksize 12000 (when the number of genomes > 20000) 

# Taxonomic analyses


##  install
conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1

##   wget the gtdb DATABASE
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz

## use
conda activate gtdbtk-2.1.0
conda env config vars set GTDBTK_DATA_PATH=home/ec2-user/data1/DB/gtdb/release202
gtdbtk classify_wf --cpus 20 \
 --genome_dir /SMAG/ --out_dir /GTDB_SMAG/ -x fa \

# tRNA & rRNA identification

##  install
conda install tRNAscan-SE
conda install barrnap

##  use
##  trna
for bacteria: for i in ${ID[@]};do tRNAscan-SE -B -Q -o /trna/${i}trna.out /bac_mag/${i}.fa;done
for archae: for i in ${ID[@]};do tRNAscan-SE -A -Q -o /trna/${i}trna.out /arc_mag/${i}.fa;done

##  rrna

for bacteria: for i in ${ID[@]};do barrnap --quiet -k bac /bac_mag/${i}.fa --outseq /rrna_hit/${i}hit.fa > /rrna_out/${i}rrna.out;done
for archae: for i in ${ID[@]};do barrnap --quiet -k arc /arc_mag/${i}.fa --outseq /rrna_hit/${i}hit.fa > /rrna_out/${i}rrna.out;done


# Phylogenetic analyse

##  install
conda install -c bioconda phylophlan
##  use
##  construct your own config_file (set the parameters)
ex:phylophlan_write_config_file  \
    -d a \
    -o newphylo.cfg \
    --force_nucleotides \(based on the input)
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 iqtree \
    --tree2 raxml \
    --overwrite \
    --verbose 2>&1

##  constructtree
phylophlan -i  SMAG_drep/ -d phylophlan --diversity high --accurate \
-f /config/newphylo.cfg \
-o /tree/SMAG_drep_tree --nproc 40 --verbose 2>&1 | tee tree.log


# Abundance of MAGs
##  install
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar -jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
bwa -h
## use
##  rename the contig name in the mag
sed -i 's/>k141/>$samplein.11k141/g' $samplebin.1.fa
##  cat the mag(after rename)
cat ./*.fa > catmag.fa
##  index mag
bwa index ./catmag.fa 
##  mapping the read to magindex
ID=(sample)
for i in ${ID[@]};do
    /bwa-0.7.17/bwa mem -t 96 /catindex/catmag.fa \
/02.paired_data/${i}_paired_1.clean.fastq.gz \
/02.paired_data/${i}_paired_2.clean.fastq.gz| samtools sort -O bam -@ 96 \
-o - > /$result_dir/$3.bam;done
##  filter the read-percent-identity and read-aligned-percent

##  install
conda install coverm
##  use
ID=(sample)
for j in ${ID[@]}; do coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 \
-b /bam/${i}.bam -o /bam/${j}_filter.bam -t 32;done
##  Depth of coverage
ID=(sample)
for j in ${ID[@]}; do coverm contig --trim-max 90 --trim-min 10 --methods trimmed_mean \
 --bam-files /bam/${j}_filter.bam > /csv/${j}_filter.bam_coverage.csv;done
##  weight by the contig length in the mag

##  # rscript
 Rscript
library(tidyr)
library(data.table)
## read the mapping coverage csv
csv<-list.files("~/data1/LUCY/mapping/catdone",pattern="*coverage.csv$") 
## read the contig length
contig_length<-fread("/media/contig.txt")
colnames(contig_length)[1]<-"mag"
 colnames(num_contig)[2]<-"mag"
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
### weight by the contig length
abundance<-catcoverage1$V2/catcoverage1$V1
catcoverage1$ave_abundance<-abundance
write.csv(catcoverage1,file=paste("/media/length/cat_",i,sep=""))


# MAG function annotation
##  install
conda install -c bioconda mmseqs2

##  Protein clustering
mmseqs easy-linclust /catmag/mag.faa clusterRes tmp --min-seq-id 0.95 -c 0.9 --threads 96 --cluster-mode 2 --cov-mode 1
~clusterRes_rep_seq.fasta is representative amino acids
##  Functiona annotation
##  install
python>=3.7
git clone https://github.com/jhcepas/eggnog-mapper.git
cd eggnog-mapper
./download_eggnog_data.py 
##  use
python3 /eggnog-mapper-2.1.7/emapper.py  \
-i /catmag/mag_repre.fasta \
--report_orthologs --output /function/mag_repre -m diamond --cpu 36 


# BGC identification

##  install
conda create -n antismash antismash
conda activate antismash
download-antismash-databases
conda deactivate
##  use
antismash -c 96 --cb-general --cb-knownclusters --cb-subclusters --cc-mibig --genefinding-tool prodigal --allow-long-headers /catmag/mag_repre.fasta

##  BGC classfication
##  install
conda env create -f environment.yml
source activate bigscape
conda install numpy scipy scikit-learn
conda install -c bioconda hmmer biopython fasttree
conda install -c anaconda networkx
mdkir BiG-SCAPE
cd BiG-SCAPE
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz 
hmmpress Pfam-A.hmm
python bigscape.py --version
may error

##  download the antiSMASH databases 
download_antismash_databases

##  attention
vim environment.yml
insert
conda create --name bigscape_fix python=3.7 scikit-learn=0.19.2
conda install biopython networkx
conda install -c bioconda hmmer fasttree

##  use
python /BiG-SCAPE/bigscape.py \
-i /all_gbk/mag_gbk \
-o /BGC/BGC_class \
-c 80 --cutoffs 0.3 --clan_cutoff 0.2 0.8 --mibig 

# SNV catalog construction

##  01.snv-calling & shared-snp.list

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
## get share snp List files
ls $OUTPUT_DIR2 | grep snps | cut -d'.' -f1 | xargs -I[] bash -c 'sed "1,5d" '$OUTPUT_DIR2'/[].snps | awk "$0"' '$2 != "." && $3 != "." {printf "%s\t%s||%s||%s||%s\n", "[]", $14, $1, $2, $3}' | cut -f2 | LC_ALL=C sort -k2,2 -S20G --parallel=4 | uniq -c | awk '$1 > 1 {print $2}' > $SPECIES.snps.list

##  02.Snv catalog
##  # method 1
python3 generate_catalog_revised.py --shared 112.sp001002135.snps.list --in-list <(ls /03.snv/101-150/112.sp001002135_CATALOG/snps/*.snps) --out 112.sp001002135.catalog.tsv
###  method 2
SNV_CATALOG.R

# Pangenome construction
##   genomes selestion
completeness>=80%
contamination<=5%
HQ genomes >=10 
##   fasta to gff (prokka v1.14.5)

##  install
conda create -n prokka
conda install -c bioconda prokka
##  use
prokka ./mag.fa ./mag.gff

##   create a pangenome with a minimum amino acid identity at 90% (‘-i 90’) and a core gene defined at 90% presence (‘-cd 90’)
##  install
conda create -n roary roary
##  use
roary -e --mafft -i 90 -cd 90 -f  output_dir  *.gff

# Crispr-cas system information

##  crispr array
# #  # contig filter
conda activate pangenome
for i in $(cat sample.list); do perl 0.casfind-1.pl /${i}.fa /2.5eur/${i}.3kb.fa /cas/${i}.temp;done

### crispr array and protein prediction

####  pilercr install
git clone git@github.com:widdowquinn/pilercrpy.git
cd pilercrpy
python setup.py install

### use
for i in $(cat sample.list);do pilercr -in /cas/${i}.temp -out /cas/${i}.spacer;done
for i in $(cat sample.list);do prodigal -a /cas/${i}.pep -i /cas/${i}.temp  -p single -f gff -o /cas/${i}.temp/${i}.temp.gff;done

###  protein filter 10kb region of crispr array
for i in $(cat sample.list);do perl 0.casfind-2.pl /cas/${i}.spacer /${i}.fa;done
for i in $(cat sample.list);do perl 0.casfind-3.pl /cas/${i}.temp.gff /${i}.fa;done
for i in $(cat sample.list);do bedtools intersect -wo -a /cas/${i}.temp.gff.loc -b /cas/${i}.spacer.loc > /cas/${i}.bed;done
for i in $(cat sample.list);do perl 0.casfind-4.pl /cas/${i}.temp.bed /cas/${i}.pep /cas/${i}.temp.pep.fasta /cas/${i}.pep.cas.fasta /${i}.fa;done

## overlap region predication

####bedtools install
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz 
tar zxvf bedtools-2.30.0.tar.gz
cd bedtools2 
make 
### use
bedtools intersect -wo -a ~/cas/output/${i}.temp.gff.loc -b ~/cas/output/${i}.spacer.loc > ~/cas/output/${i}.bed
### use
perl 0.casfind-4.pl ~/cas/output/${i}.bed ~/cas/output/${i}.temp.pep ~/cas/output/${i}.temp.pep.fasta ~/cas/output/${i}.pep.cas.fasta ~/cas/genome/${i}.fa

##  cas protein identification and alignment

### du-replicated
####  cdhit install
wget https://github.com/weizhongli/cdhit/releases/download/V4.6.7/cd-hit-v4.6.7-2017-0501-Linux-binary.tar.gz
tar -zxvf cd-hit-v4.6.7-2017-0501-Linux-binary.tar.gz

####  use
cat /cas/*.pep.cas.fasta > /casdb/Casdb.fa
cd-hit -i /casdb/Casdb.fa -o /casdb/NR/Casdb.0.95.fa -c 0.95 -aS 0.9 -n 5 
### filter protein length (200-1000aa)
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

### (cas novelty) align to NR database
####  install
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.31/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
####  use
/diamond blastp -k 1 -e 1e-10  --query-cover 40 -d /uniref_100.dmnd -q /casdb/NR/Casdb.0.95.fa -o /casdb/NR/CAS


# Connecting MAGs to viruses
## Linkages from Crispr-cas system information

##  Remove MAGs containing fewer than four CRISPR-associated proteins.
##  https://github.com/sandialabs/CasCollect/tree/master/ref
##  Download all the hmm files of cas proteins
prodigal -i mag10k.fa -a mag10k.faa -d mag10k.ffn -p meta -f gff > 
hmmsearch -Z 1 --noali --cpu 10 --tblout Cas3_0_ID.out Cmr3_1_IIIB.hmm MAG_protein.faa

##  filter CRISPR-associated proteins >=4
## R
library(data.table)
casout <- read.csv("cas.out",comment.char="# ",sep="",header=F)
casout <- casout[casout$V5<1e-5,]
proteins <- unique(casout$V1)
bins <- gsub("_\\S+","",proteins)
count <- data.frame(table(bins))
count <- count[count$Freq>=4,]
data <- data.frame(bins,proteins)
data <- data[data$bins%in%count$bins,]
seqs <- unique(gsub("_\\d+$","",data$proteins))
writeLines(seqs,"cas_bins_name.txt")

##  sequence grep
conda install seqkit
/public/ylwang/seqkit grep -f cas_bins_name.txt MAG_renamecatall.fa > MAG_cas.fa

##  Contigs longer than 10 kb
seqkit seq -m 10000 MAG_cas.fa > MAG_cas10k.fa

##  minced(0.4.2)
### install
wget https://github.com/ctSkennerton/minced/archive/0.4.2.tar.gz
tar -zxvf 0.4.2.tar.gz
rm  0.4.2.tar.gz
cd minced-0.4.2/
make
##  use
/minced-0.4.2/minced -gffFull -spacers MAGnoderup_10k.part_001.fasta MAGnoderup_10k.part_001.txt MAGnoderup_10k.part_001.gff

##  CRT crispr
### install
http://www.room220.com/crt/
wget http://www.room220.com/crt/CRT1.2-CLI.jar.zip
unzip CRT1.2-CLI.jar.zip
### use
java -cp /application/CRT1.2-CLI.jar crt MAGnoderup_10k.part_030.fasta MAGnoderup_10k.part_030.out

##   virus Prediction
##  install
conda create -n vs2 virsorter=2
conda activate vs2
##  use
virsorter2.sif run -w test -i mag.fa -j 4 all

##   Quality Control
### install
conda install -c conda-forge -c bioconda checkv
### download the database
checkv download_database .
### download the database manually
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
export CHECKVDB=/path/to/checkv-db-v1.0
### use
checkv end_to_end input_file.fna output_directory -t 28

##  prophage
makeblastdb -in GSV_2022.fasta -dbtype nucl
nohup blastn -query MAG_renamecatall.fa -db GSV_2022.fasta -out prophage.out -outfmt 6 -num_threads 200 &
awk '$3 >90 && $4 > 500' prophage.out > prophage_500_90.out
seqkit fx2tab -n -l GSV_2022.fasta > GSV_2022_nl.txt
seqkit fx2tab -n -l MAG_renamecatall.fa > MAG_renamecatall_nl.txt


# Statistical analysis and visualization

Some processing steps in the pipeline, statistical analysis were handled by scripting with R, Shell, Perl or Python languages. These scripts were placed in "[Scripts](https://github.com/Caiyulu-818/SMAG/blob/main/Pipeline)" directory. All related input data for visualization are in ["figure"](https://github.com/Caiyulu-818/SMAG/scripts) directory.




 


