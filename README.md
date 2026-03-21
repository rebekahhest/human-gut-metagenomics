# human-gut-metagenomics
Taxonomic classification, community comparison and differential abundance analysis for human gut microbiomes with different diets

**Author:** Rebekah Hest  
**Course:** BINF6110 - Genomic Methods for Bioinformatics

## Assignment 3

This repository contains Assignment 3, a study to explore and compare taxonomic communities between microbial samples from shotgun metagenomics data of human gut microbiomes with different habitual diets including vegans and omnivores.

### Introduction

The human gut microbiome is a complex community that contributes to host metabolism and immune function by synthesizing vitamins and producing necessary metabolites that impact overall health (Toma et al., 2026). These diverse microbiota profiles consist primarily of bacteria but also include archaea, viruses, and other microorganisms where no single preferred composition has been determined because optimal functioning for disease prevention is specific to each individual (Rinninella et al., 2019). However, long- and short-term gut microbial ecology has been found to be largely modulated by diet and thus diet selection can be studied as prevention and treatment strategies for common diseases (Huang et al., 2024).
<br><br>
Plant-based diets tend to produce fiber-degrading microbes which provide anti-inflammatory effects and cardiovascular protection (Tomova et al., 2019). In contrast, a diet rich in animal foods can lead to protein fermentation which can be detrimental for the host’s health and contribute to can bowel diseases (Fackelmann et al., 2025). Characterization of microbial communities can be conducted via targeted 16S rRNA gene sequencing or shotgun metagenomics. Compared with 16S rRNA gene sequencing, shotgun metagenomics provides greater taxonomic resolution and microbial coverage and can reveal more biologically informative taxa (Durazzi et al., 2021). Furthermore, it enables direct characterization of the functional potential of the microbiome (Joval et al., 2016). Therefore, the objective of this study is to compare gut microbial diversity and taxonomic composition between individuals consuming vegan and omnivorous diets using shotgun metagenomic sequencing data.
<br><br>
Taxonomic classification and species-level abundance estimation in metagenomic datasets can be performed using several computational approaches including alignment-based, marker gene-based, and k-mer based methods. Alignment-based approaches like BLAST offer high accuracy and specificity for detecting distant homologous relationships but require substantial computational resources for large reference databases (Altschul et al., 1990). Marker-based gene methods, such as metaPhlAn, instead classify reads using a curated set of clade-specific marker genes, reducing computational requirements but limiting detection to organisms represented by these markers (Blanco-Minguez et al., 2023). In contrast, k-mer based approaches assign reads based on exact matches of short sequence fragments, enabling rapid classification of large metagenomic datasets (Wood & Salzberg, 2014). Although k-mer–based classifiers such as Kraken2 may produce false positive assignments due to shared k-mers among related taxa, KrakenUniq attempts to mitigate this issue by tracking unique k-mers, but its substantially higher computational requirements often limit its practical use (Breitwieser et al., 2018). Kraken2 classifications can also be refined using Bracken to improve species-level abundance estimation through Bayesian re-estimation (Lu et al., 2017). Benchmarking studies have demonstrated that k-mer–based classifiers such as Kraken2 achieve high classification accuracy while maintaining substantially faster runtimes than other approaches (Ye et al., 2019).
<br><br>
To conduct Differential Abundance (DA) analysis for microbiome studies, several methods have been developed to address the challenge of handling the compositional nature of sequencing data. ANOVA-Like Differential Expression tool (ALDEx2) estimates relative abundance using a Bayesian estimation by Monte Carlo sampling (Fernandes et al., 2014). 
Analysis of composition of microbes (ANCOM-BC) applies bias correction to estimate log fold changes while accounting for unequal sampling fractions (Li and Peddada, 2023). ALDEx-2 has been shown to exceed the nominal level of False Detection Rate (FDR) (5%) and result in substantially smaller compared to other DA methods, whereas ANCOM-BC performs reasonably well for controlling FDR and maintains high recall power (Cappellato, Baruzzo, and Di Camillo, 2022; Li and Peddada, 2020). 

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

#### Taxonomic Composition

Across all samples, the gut microbiome was dominated by Bacteriodota and Bacillota, which account for the majority of relative abundance together in both dietary groups (Figure 1). Omnivore samples showed greater variability in phylum-level composition and a more balanced distribution of taxa, with one sample (SRR8146936) exhibiting a particularly high proportion of Bacillota. In contrast, vegan samples were more consistent in relative abundances among samples, except for one sample (SRR814951) that was strongly represented by Bactillota. Minor contributions from other phylum, including Uroviricota and Mycoplasmatota in lower relative abundance were observed across all samples.
<br><br>
However, at the genus level, the differential abundance analysis using ANCOM-BC2 identified no taxa that were significantly different between dietary groups after multiple testing correction via the Holm-method (q > 0.05). Despite the lack of statistically significant results, several genera exhibited notable differences in log fold change between groups (Figure 2). Genera such as _Romboutsia_, _Parvimonas_, and _Anaerobutyricum_ showed higher relative abundance in vegan samples, whereas taxa including _Fusobacterium_ and _Campylobacter_ were more abundant in omnivore samples. The top twenty genera exhibiting the largest log-fold change primarily represented the phylum Bactillota followed by Bacteroida. 

#### Microbial Community Structure
Alpha diversity was assessed using Shannon diversity in which there was no significant different between omnivore samples (mean ~ 3.18) and vegan samples (mean ~ 2.92) when evaluated using a Wilcox rank sum test (p = 0.7) (Figure 3). Similarly, Simpson diversity found no significant difference between dietary groups, omnivores (mean ~ 0.85) and vegans (mean ~ 0.79) using a Wilcox test (p = 0.7) (Figure 4). Both metrics displayed overlapping distributions, indicating comparable within-group diversity across diets (Figure 5). Between both diversity measures, the omnivore samples contained both high and low alpha diversity estimates, whereas the vegan samples were more consistent across samples.
<br><br>
The beta diversity measured using Bray-Curtis dissimilarity and PCoA revealed no distinct clustering of samples by group despite diet explaining 76.8% of the variation (PC1 = 51% and PC2 = 25.8%) (Figure 6). Along the first principal component, the samples were interspersed but the second principal component showed slight separation of the dietary groups. The PERMANOVA supported this showing that diet explained only 18.5% of the observed variation with no statistically significant difference between groups (R<sup>2</sup> = 0.185, p = 0.5). Together, these results suggest substantial overlap in composition exists between the groups and that diet was not a strong determinant of microbial community structure.

### Discussion

This study investigated the taxonomic compositional differences in human gut microbiome between individuals on an omnivore and vegan diet using shotgun metagenomic sequencing. The results indicate that diet was not a strong predictor of microbial variation in the dataset. The relative abundance profiles were similar across all samples, with both groups dominated by only a few phyla. In particular, Bactilla and Bacteriodota together comprised roughly 75% of the classified phyla for all samples, with varying ratios of the two taxa observed among the samples.
<br><br>
The alpha diversity analysis further supported this observation as neither Shannon nor Simpson indices showed no significant differences between dietary groups. This suggests that within-sample diversity and evenness were comparable between vegan and omnivore microbes. High alpha diversity has been observed to increase microbiome productivity and reduced disease risk (Bell et al., 2005). However, several diseases have been characterized by elevated diversity, indicating that alpha diversity alone is not a reliable indicator of host health (Williams et al., 2024). Therefore, the lack of significant differences in alpha diversity between dietary groups should be interpreted cautiously and does not necessarily reflect differences in health between diets.
<br><br>
Beta diversity metrics are often considered more sensitive than alpha diversity for detecting differences between microbial communities, particularly when using abundance-based measures such as Bray–Curtis dissimilarity (Kers and Saccenti, 2022). Despite this, no clear clustering of samples by dietary group was observed, suggesting that differences in community composition may have been subtle relative to within-group variability. However, beta diversity measures have been known to violate mathematical properties and can compromise downstream analyses including PCoA and PERMANOVA tests (Zhu et al., 2026). Coupled with the small sample size, these limitations may reduce the ability to accurately detect finer-scale differences in microbial community structure between dietary groups and may require a framework that further refines dissimilarity metrics.
<br><br>
Limitations in differential abundance detection may arise from the taxonomic classification using Kraken2 and abundance re-estimation with Bracken where ambiguous reads of higher taxonomic rank may be preferentially assigned to more common taxa. Consequently, rare taxa may be underrepresented in downstream analyses (Xu, Rajeev, and Salvador, 2023). ANCOM-BC2 may further reduce sensitivity to rare taxa, as it performs more reliably on common taxa and may have reduced power for low-abundance groups (Shi et al., 2024). Together, these factors bias the analysis toward dominant taxa and may not reliably capture the entire microbial community.
<br><br>
Although no taxa were significantly different between dietary groups, several genera exhibited notable differences in log fold change. For example, _Fusobacterium_ and _Campylobacter_ were more abundant in omnivore samples, while genera such as _Romboutsia_ and _Anaerobutyricum_ were relatively enriched in vegan samples. Previous studies have associated certain genera with dietary patterns, including increased abundance of short-chain fatty acid–producing taxa in plant-based diets (Egas-Montenegro et al., 2026; Soldan et al., 2024). While these observations are consistent with broader trends reported in the literature, they should be interpreted cautiously given the lack of statistical significance and limited sample size.
<br><br>
Overall, these results highlight the importance of sample size, tool selection and metric interpretation for microbiome studies. Additional research on microbiota composition and diversity is required to concretely determine the role and mechanisms by which diet influences host health.

## References
