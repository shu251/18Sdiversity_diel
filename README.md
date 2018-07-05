## A Hard Day's Night: Diel shifts in microbial eukaryotic activity in the North Pacific Subtropical Gyre

Current citation: Sarah K. Hu*, Paige E. Connell, Lisa Y. Mesrop, & David A. Caron. (In prep). A Hard Day’s Night: Diel shifts in microbial eukaryotic activity in the North Pacific Subtropical Gyre.

### Summary
We investigated diel shifts in community composition (based on rDNA sequences) and relative activity (based on rRNA sequences) of naturally occurring protistan assemblages (18S V4 tag sequencing of rDNA and rRNA, respectively) in the North Pacific Subtropical Gyre (NPSG) at 4 hour intervals over a period of 3 days (Lagrangian sampling).
This study is a part of the Simon’s Collaboration on Ocean Processes and Ecology (SCOPE). [More information here](http://scope.soest.hawaii.edu/)

#### Methods & Dataset
* Sample collection
* [DNA/RNA extraction protocol](https://www.protocols.io/view/rna-and-optional-dna-extraction-from-environmental-hk3b4yn)
* [Library preparation](https://www.protocols.io/view/18s-v4-tag-sequencing-pcr-amplification-and-librar-hdmb246)
* 18S rRNA gene [tag-sequencing pipeline]((https://github.com/shu251/https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1)
* Raw sequence data are available under SRA BioProject PRJNA393172
* Filtered/QC-ed sequences available here:  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.846380.svg)](https://doi.org/10.5281/zenodo.1243295)
* Original OTU table (output from QIIME) also available [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.846380.svg)](https://doi.org/10.5281/zenodo.1243295)

#### Sequence QC
Script for processing raw fastq files from QC to OTU clustering: diel_QC_04282018.pl

*Summary of steps:*
1. Merge with fastq-join (QIIME1)
2. Filter low quality reads in QIIME1, Q>30
3. Remove V4 primers (cutadapt)
4. Filter sequences that are too short or too long
5. Merged individual sample .fasta files into one
6. Chimera check (reference based, PR2) using vsearch
7. OTU clustering in QIIME1, open-reference
8. Taxonomy assignment using uclust in QIIME1 and PR2 database (v.4.7)

#### Import into R for analysis of data
See [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.846380.svg)](https://doi.org/10.5281/zenodo.1243295) to download full OTU table.
Rscript included here, 'Rscript_Diel18S_SHu_05082018.r' imports this OTU table and generates all figures found in manuscript.

#### R script annotations
R script is broadly sectioned out for each part of data analysis


*Set-up:*
Import raw OTU table, get sequence run information, update time of day label, manually curate taxonomic group names, and generate Figure S2 (sequence run summary and distribution of OTUs).
Finally, randomly subsample (rarefy) samples in both RNA and DNA libraries so they have the same number of sequences (per sample).

*Diversity analysis*
* For Figure 2 in the main text, extract the total number of RNA or DNA OTUs from the subsampled dataset.
* To generate Figure 3 bar plot, sum the total number of OTUs belonging to each taxonomic group (OTU richness for each taxonomic group)
* Generate Figure 4 area plot by summing all sequences within a taxonomic group. Ahead of plotting, remove Unassigned, NA, and metazoan groups.
* Calculate and plot RNA:DNA read ratios. First, calculate RNA:DNA ratio for each OTU. Remove OTUs which did not have RNA or DNA sequences (meaning, only keep OTUs that had both RNA and DNA sequences). Calculate the mean and standard mean error of RNA:DNA ratios within each taxonomic group.

*Statistical analysis*
* RAIN analysis: center log-ratio transform data and perform RAIN analysis on RNA and DNA library results separately. Compile results so that significantly diel OTUs have p values < 0.05.
* Prepare input table for extended Local Similarity analysis. First funnel OTUs, so only OTUs found in every sample and which had more than 10 sequences are retained. Perform centered log-ratio transformation. Export .txt file to conduct eLSA [See tutorial here](https://stamps.mbl.edu/index.php/Local_Similarity_Analysis_(LSA)_Tutorial). The following section of the R script imports results from eLSA to extract significantly co-occurring OTUs and summarize findings by taxonomic group.

*Additional supplementary plots*
* Plot environmental data
* Calculate Bray-Curtis dissimilarity matrix and generate a two-dimensional MDS plot to visualize results.

##### Last updated by S. Hu 06-04-2018
