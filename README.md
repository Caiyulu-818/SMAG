# SMAG
Soil microbial dark matter mining based on genome-resolved metagenomics 

We constructed the SMAG from 3,463 soil metagenomes from large-scale genome-resolved metagenomics to expand the genomic catalog of soil microbiomes !

After download, all datasets can be unpacked using: `tar -xjvf <MAG.tar.gz>`

## SMAG catalog

Soil MAGs were assembled from 3,463 soil metagenomes from 9 different ecosystems across globe. all the MAGs were recovered for individual metagenomic assemblies using three different tools with default options: [Metabat (v2.12.1)](https://github.com/bioboxes/metaBAT), [MaxBin (v2.2.6)](https://github.com/movingpictures83/MaxBin), [CONCOCT (v1.0.0)](https://github.com/ConcoctLang/concoct), which all on the basis of tetranucleotide frequencies(TNF) and coverage information. The resulting MAGs were refined using the module ‘bin_refinement’ from [metaWRAP (v1.2.1)](https://github.com/bxlab/metaWRAP).

<b>SMAG MAGs (N=40039)</b>   
download from your browser: [https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu/MAG.tar.gz]

<b>SMAG MAGs after dereplication (N=21077)</b>
download from your browser:[https://data.cyverse.org/dav-anon/iplant/home/lucyzju/Caiyu/MAGdrep.tar.gz]

* All MAGs are estimated to be >50% complete and <10% contaminated
* The MAGs after dereplication meet or exceed the medium-quality level of the minimum information about a metagenome-assembled genome (MIMAG).
