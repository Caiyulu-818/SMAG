#SMAG documation

All the information related the study can be viewed at https://microbma.github.io/project/SMAG.html

Soil microbial dark matter mining based on genome-resolved metagenomics

We constructed the SMAG from 2,990 soil metagenomes from large-scale genome-resolved metagenomics to expand the genomic catalog of soil microbiomes!

After download, all datasets can be unpacked using: cat ./mag.tar.gz* > mag.tar.gz; tar -xjvf <mag.tar.gz>

##SMAG catalog
Soil MAGs were assembled from 2,990 soil metagenomes from 9 different ecosystems across globe. All the MAGs were recovered for individual metagenomic assemblies using three different tools with default options: Metabat (v2.12.1), MaxBin (v2.2.6), CONCOCT (v1.0.0), which all on the basis of tetranucleotide frequencies(TNF) and coverage information. The resulting MAGs were refined using the module ‘bin_refinement’ from metaWRAP (v1.2.1).

##SMAG MAGs (N=40,039, MAG.tar.gz)
To upload the large gz file, we split it into smaller file with the prefix "mag.tar.gz", downloaders can use the cat ./mag.tar.gz* > mag.tar.gz; tar -xjvf <mag.tar.gz> to process the MAGs.

##MAG MAGs after dereplication (N=21,077,MAGdrep.tar.gz) All the resource was deposited at cyverse [/iplant/home/lucyzju/Caiyu/]

All MAGs are estimated to be >= 50% complete and < 10% contaminated

The MAGs after dereplication meet or exceed the medium-quality level of the minimum information about a metagenome-assembled genome (MIMAG).

If you need to download please contact lucy20@zju.edu.cn for permission.

##RAxML_bestTree.mag20177prodigal_refined.tre

The dereplicated soil MAGs phylogenetic tree built by phylophlan.

##magvirus.fa

The virus sequence annotated from the SMAG.

##SNV_CATALOG.tar.gz

The SNV catalog constructed from the SMAG.

##Extended Data Table 1-6.xlsx

The supplement tables related to the study.