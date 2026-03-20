# human-gut-metagenomics
Taxonomic classification, community comparison and differential abundance analysis for human gut microbiomes with different diets

**Author:** Rebekah Hest  
**Course:** BINF6110 - Genomic Methods for Bioinformatics

## Assignment 3

This repository contains Assignment 3, a study to explore and compare taxonomic communities between microbial samples from shotgun metagenomics data of human gut microbiomes with different habitual diets including vegans, vegetarians and omnivores.

### Introduction

The human gut microbiome is a complex community that contributes to host metabolism and immune function by synthesizing vitamins and producing necessary metabolites that impact overall health (Toma et al., 2026). These diverse microbiota profiles consist primarily of bacteria but also include archaea, viruses, and other microorganisms where no single preferred composition has been determined because optimal functioning for disease prevention is specific to each individual (Rinninella et al., 2019). However, long- and short-term gut microbial ecology has been found to be largely modulated by diet and thus diet selection can be studied as prevention and treatment strategies for common diseases (Huang et al., 2024).
<br><br>
Plant-based diets tend to produce fiber-degrading microbes which provide anti-inflammatory effects and cardiovascular protection (Tomova et al., 2019). In contrast, a diet rich in animal foods can lead to protein fermentation which can be detrimental for the host’s health and contribute to can bowel diseases (Fackelmann et al., 2025). Characterization of microbial communities can be conducted via targeted 16S rRNA gene sequencing or shotgun metagenomics. Compared with 16S rRNA gene sequencing, shotgun metagenomics provides greater taxonomic resolution and microbial coverage and can reveal more biologically informative taxa (Durazzi et al., 2021). Furthermore, it enables direct characterization of the functional potential of the microbiome (Joval et al., 2016). Therefore, the objective of this study is to compare gut microbial diversity and taxonomic composition between individuals consuming vegan and omnivorous diets using shotgun metagenomic sequencing data.
<br><br>
Taxonomic classification and species-level abundance estimation in metagenomic datasets can be performed using several computational approaches including alignment-based, marker gene-based, and k-mer based methods. Alignment-based approaches like BLAST offer high accuracy and specificity for detecting distant homologous relationships but require substantial computational resources for large reference databases (Altschul et al., 1990). Marker-based gene methods, such as metaPhlAn, instead classify reads using a curated set of clade-specific marker genes, reducing computational requirements but limiting detection to organisms represented by these markers (Blanco-Minguez et al., 2023). In contrast, k-mer based approaches assign reads based on exact matches of short sequence fragments, enabling rapid classification of large metagenomic datasets (Wood & Salzberg, 2014). Although k-mer–based classifiers such as Kraken2 may produce false positive assignments due to shared k-mers among related taxa, KrakenUniq attempts to mitigate this issue by tracking unique k-mers, but its substantially higher computational requirements often limit its practical use (Breitwieser et al., 2018). Kraken2 classifications can also be refined using Bracken to improve species-level abundance estimation through Bayesian re-estimation (Lu et al., 2017). Benchmarking studies have demonstrated that k-mer–based classifiers such as Kraken2 achieve high classification accuracy while maintaining substantially faster runtimes than other approaches (Ye et al., 2019).
### Methods

#### Data Acquisition

Shotgun metagenomic data from the gut microbiome of healthy Italian individuals was downloaded from the Sequence Read Archive (NCBI-SRA: SRP126540). Six samples were selected based on diet (omnivore n = 3; vegan n = 3). Raw sequencing reads were downloaded using the SRA Toolkit (v3.2.1) with commands `prefetch` and `fasterq-dump` to retrieve and convert SRA files to FASTQ	 format. Taxonomic classification was performed using the Kraken2 Core_nt reference database (`k2_core_nt_20251015.tar.gz`), which includes a comprehensive collection of nucleotide sequences derived from GenBank, RefSeq, Third Party Annotation (TPA), and Protein Data Bank (PDB) entries.

#### Quality Control

Initial quality assessment of raw sequencing reads was performed using FastQC (v0.12.1) and summarized with MultiQC (v1.13) using the commands `fastqc` and `multiqc`. Quality metrics including per-base sequence quality, GC content, and sequence length distribution were evaluated. Multiple peaks were observed in the Per Sequence GC Content plots, suggesting the presence of adapter contamination or mixed sequence composition. Therefore, adapter trimming and quality filtering were performed using Trimmomatic (v0.39) in paired-end mode (`PE`). Nextera adapter sequences were removed using `ILLUMINACLIP:NexteraPE-PE.fa:2:40:15`, compatible with the Nextera DNA library preparation kit used for NextSeq 500 sequencing. Additional quality filtering was performed using default parameters`LEADING:3`, `TRAILING:3`, `MINLEN:36`, and `SLIDINGWINDOW:4:20` to remove low-quality bases and short reads (Bolger, Lohse, & Usadel, 2014). Post-trimming sequence quality was reassessed using `fastqc` and `multiqc` to confirm improvement in read quality metrics.

#### Taxonomic Classification

Quality-filtered paired-end reads (`*_paired.fastq`) were taxonomically classified using Kraken2 (v 2.1.6) with the `--paired` option against the Core_nt reference database (Wood, Lu, & Langmead, 2019). A confidence threshold of `--confidence 0.15` was applied to reduce false positive taxonomic assignments. The option `--use-names` was included to report taxonomic classifications using scientific names rather than taxonomic identifiers. Memory mapping (`--memory-mapping`) was not used because the database was loaded directly into memory on a high-memory compute node.

#### Abundance Estimation

Species-level abundance estimates were generated using Bracken (v3.0), which refines Kraken2 taxonomic assignments through Bayesian re-estimation (Lu et al., 2017). Kraken2 report files were used as input, and abundance estimates were calculated using a read length parameter of `-r 150` bp and taxonomic level `-l S` to produce species-level abundance profiles.

#### Data Pre-Processing

Bracken-adjusted taxonomic reports were converted into BIOM format using `kraken-biom`, which aggregates taxonomic abundance estimates across samples into a standardized BIOM table (Daboub, 2016). The BIOM file was imported into R (v.2026.01.1+403) using the `biomformat` package and integrated with sample metadata using the phyloseq package (v.1.52.0) to construct a phyloseq object containing taxonomic abundance data and associated sample metadata (McMurdie & Holmes, 2012). The metadata file (`metadata.tsv`) was adapted from the SRA metadata and included dietary classification (omnivore or vegan) and location (city) for each sample. This phyloseq object was used for downstream taxonomic composition summaries and diversity analyses. Rarefaction curves were generated to assess sampling effort. 

#### Relative Abundance Analysis

Relative taxonomic composition was evaluated at the phylum level by aggregating taxa using the `tax_glom` function in the phyloseq package (v.1.52.0). Relative abundances were visualized as stack bar plots per sample and grouped by dietary group to compare community composition. The top ten phylum by total abundance were retained and the remaining phylum were grouped into “Other” to improve interpretability. 

#### Alpha Diversity Analysis (Within-samples)

Alpha diversity was assessed using the Shannon and Simpson diversity indices, which capture both species richness and evenness within samples. Diversity metrics were calculated using the `estimate_richness` function in the phyloseq package (v.1.52.0). Diversity measures focused on rare taxa were excluded from this analysis since these estimates rely heavily on singletons but were absent following abundance re-estimation (Deng, Umbach, & Neufeld, 2024). Within-group differences were tested using a Wilcox rank-sum test for both Shannon and Simpson indices.


#### Beta Diversity Analysis (Between-samples)

Beta diversity was evaluated using Bray–Curtis dissimilarity, which incorporates differences in taxon abundances between samples. Ordination was performed using Principal Coordinates Analysis (PCoA) based on Bray–Curtis distances, implemented with the `ordinate` function in phyloseq. Differences in community composition between dietary groups were statistically assessed using permutational multivariate analysis of variance (PERMANOVA) with 999 permutations via the `adonis2` function in the vegan (2.7.3) package (Oksanen et al., 2026).

#### Differential Abundance Analysis 

To identify statistically significant differences in abundance of taxa between dietary groups, Analysis of Compositions of Microbes with Bias Correction (ANCOM-BC2) was conducted. This methodology accounts for sampling bias common in microbe studies while controlling the False Discovery Rate (FDR) (Lin & Peddada, 2020). The analysis was performed at the genus level (`tax_level = "Genus"`) using the `ancombc2` function from the ANCOMBC package (v.2.10.1), with died specified as the fixed effect (fix_formula = “Diet) and p-values were adjusted using the Holm Method (`p_adj_method = "holm"`). Taxa lacking genus-level classification were assigned to the lowest available taxonomic rank (e.g., family) to improve interpretability. Results were visualized for the top twenty taxa by absolute log-fold change. 

### Results

### Discussion

## References
