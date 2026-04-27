# SMAG
Soil microbial dark matter mining based on genome-resolved metagenomics 

We constructed the SMAG from 3,304 soil metagenomes from large-scale genome-resolved metagenomics to expand the genomic catalog of soil microbiomes!

After download, all datasets can be unpacked using: `cat ./mag.tar.gz* > mag.tar.gz; tar -xjvf <mag.tar.gz>`

# Community-maintained download guide

To facilitate systematic discovery and convenient access to high-quality MAG resources, Dr. Rui Li from the Beijing Institute of Genomics has created an open-source index project, [Awesome MAG] (https://github.com/inspirewind/awesome-mag)

SMAG has been included in this project, together with a detailed download guide:(https://github.com/inspirewind/awesome-mag/blob/main/sources/smag/download.md)

We thank Dr. Rui Li (lirui@big.ac.cn) for organizing and maintaining this helpful community resource.

## SMAG catalog

Soil MAGs were assembled from 3,304 soil metagenomes from 9 different ecosystems across globe. All the MAGs were recovered for individual metagenomic assemblies using three different tools with default options: [Metabat (v2.12.1)](https://github.com/bioboxes/metaBAT), [MaxBin (v2.2.6)](https://github.com/movingpictures83/MaxBin), [CONCOCT (v1.0.0)](https://github.com/ConcoctLang/concoct), which all on the basis of tetranucleotide frequencies(TNF) and coverage information. The resulting MAGs were refined using the module ‘bin_refinement’ from [metaWRAP (v1.2.1)](https://github.com/bxlab/metaWRAP).

<b>SMAG MAGs (N=40,039)</b>   
The SMAG catalog of the soil metagenomes, SNV catalogs and viruses predicted from SMAG for this publication are available at [Zenodo](https://doi.org/10.5281/zenodo.7941562). To upload the large gz file, we split it into smaller file with the prefix "mag.tar.gz",
downloaders can use the `cat ./mag.tar.gz* > mag.tar.gz; tar -xjvf <mag.tar.gz>` to process the MAGs.

<b>SMAG MAGs after dereplication (N=21,077)</b>
All the 21,077 MAGs were deposited at [cyverse] (https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAGdrep.tar.gz),
All the 40,039 MAGs were deposited at [cyverse] (https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu_SMAG_catalog_2023/MAG.tar.gz)


* All MAGs are estimated to be >= 50% complete and < 10% contaminated
* The MAGs after dereplication meet or exceed the medium-quality level of the minimum information about a metagenome-assembled genome (MIMAG).
* If you need more information please contact lucy20@zju.edu.cn.
