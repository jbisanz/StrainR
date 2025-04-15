# StrainR: quantification of highly related strains in shotgun metagenomics

***

# Depreciation

Note: StrainR has been depreciated in favor of StrainR2 which can be found [here](https://github.com/BisanzLab/StrainR2).

***

## The problem

Experiments on complex defined microbial communities are empowering a new wave of hypothesis-driven research; however, the quantification of strain abundances is not always a trivial task. If selective and differential media is available and micro-organisms are robust, then culture is preferable, but what if the organisms are obligate anaerobes and/or difficult to selectively enumerate by culture? Culture-free DNA-based methods are preferable in this case, and the availability of genomes for the members of the communities provides a number of opportunities. qPCR is a valid option for smaller communities; however the time and cost involved in assay developement and run-time scale horribly with increasing community complexity. For many, sequencing becomes preferable. Amplicon sequencing of the 16S rRNA gene would be powerful tool, but requires that all members of the community have distinct 16S rRNA variable regions which will often not be possible where closely related strains are being studied. Where this is true, shotgun metagenomic sequencing becomes necessary; however, quantification is not as simple as counting mapped reads. We developed this tool to quantify strains of Eggerthella lenta within a 22-member community to overcome biases due to sequence variation within this strain. Intuitively, the more divergent the strain within the species, the higher its abundance is reported in the community wherein their is a 28-fold difference in abundance reported between strains which are present at equal concentrations in a community (**Figure 1AB**).

![Figure 1](https://github.com/jbisanz/StrainR/blob/master/figures/non_normalized.jpg)
**Figure 1. Abundance profiles as a function of sequence variation in even composition communities. (A)** A dendrogram based on the presence of unique conical Kmers present in genomes which mirrors the phylogenetic relationship of strains (UPGMA clustering of binary distance). **(B)** Some strains (W1BHI6 and CC75D52) are reported at 28-fold higher levels as a function of sequence divergence not true elevated abundance. (C) Corrected abundances using fragments per thousand unique kmers per million reads mapped stabilizes abundance estimates and reports a 1.8-fold skew from highest to lowest organism.

## The solution

There are three main considerations that must be accounted for when quantifying strains:
1. The majority of reads in the sample are not informative as they map to conserved regions in the species
2. Strains that have more divergent sequence have more mapping sites
3. Not all reference genomes are of equal assembly quality
4. Multi-copy elements, ex plasmids, skew the abundance of the organism
5. Certain library protocols, ex Nextera, may create coverage skews over the genome

To overcome these issues, we apply the following methology in the StrainR pipeline:
1. Fragment all draft genomes in silico to ~N50 of worst assembly to avoid influence of assembly quality
2. Generate a set of canonical Kmers of size approximately equivalent to read length for each sub-fragment
3. Generate sparse matrix on a per-fragment basis retaining only unique kmers on a per-fragment basis
4. Strictly filter reads to minimize expected errors
5. Map reads to fragments maintaining only perfect unambiguous matches
6. Calculate abundance on a per fragment basis using **FKM: Fragments per thousand unique Kmers per million reads mapped**
7. Calculate the abundance of the strain based on the median FKM of all of its sub-fragments

## Dependencies
Unix-type operating system (tested on OS X Mojave and Ubuntu 14.04.5)

R Packages:
* R >3.5.0
* tidyverse
* Biostrings
* doParallel
* foreach
* data.table
* Matrix
* Matrix.utils
* vegan
* openssl
* dada2
* ShortRead

Available to run via system call:
* BBmap >=37.97
* jellyfish 2
* samtools >=1.9

## Use

### Build database

This setup needs only be done once per community and can be recycled for all samples of the community. This is created using the `PreProcessR` script. Example usage:
```
Rscript PreProcessR.R \
 --indir /labmainshare/qb3share/jbisanz/ElComp/contigs/ \ # directory containing your assemblies
 --outdir MockDB \
 --threads 4 \
 --jobs 4 \
 --kmersize 120
 ```
 
### Map Sample

```
Rscript StrainR.R \
 --forward reads/Sample_R1.fastq.gz \
 --reverse reads/Sample_R2.fastq.gz  \
 --reference MockDB \
 --bbmap /labmainshare/qb3share/shared_resources/sftwrshare/bbmap/bbmap.sh \
 --threads 12 \
 --mem 24 \
 --outdir Mock_StrainR/ \
 --outprefix Sample
```

### Example Output

The reported about gives the abundance on a per-fragment basis. Using the median of these abundances is suggested to give a stable abundance estimation.

`head -n 2 Mock_StrainR/Sample.abundances`

StrainID | ContigID | Start_Stop | Unique_Kmers | Length_Contig | Bases | Coverage | Mapped_Reads | Mapped_Frags | Total_Mapped_Reads_In_Sample | FKM
---------|----------|------------|--------------|---------------|-------|----------|--------------|--------------|------------------------------|-----
Eggerthella_lenta_1356FAA | Eggerthella_lenta_1356FAA\|GL622582.1 | 3332860_3378358 | 2758 | 45498 | 184080 | 4.0459 | 1416 | 708 | 278548 | 921.5925414859568


The following example was generated in silico using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) using reads simulated based on a real 150x150 Illumina NovaSeq Run with community compositions based on predetermined proportions:

![Figure 2](https://github.com/jbisanz/StrainR/blob/master/figures/mock_barplots.jpg)
**Figure 2. Composition of in silico generated communities across sequencing depth.** Inset figures show the proportion of reads generated as a function of input community.

