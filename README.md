<div align="justify">
  
# Investigation of Empty Droplets detection in SUMseq
  
This repo is an annual project made during studies in Bioinformatics Institute by Shakir Suleimanov, Maria Lukina and Vladimir Grigoriants.

This project focuses on identifying and visualizing dependencies that lead to erroneous droplet classification in SUMseq data, particularly due to lower overall coverage in some samples and uneven UMI distribution among cell nuclei. We developed methods to detect and visualize these differences using various analytical approaches.

## Introduction to SUMseq technology

SUM-seq (Single-cell Ultra-high-throughput Multiomic sequencing) is introduced as a cost-effective and scalable sequencing technique designed for multiplexed multiomics profiling. This method enables the simultaneous profiling of chromatin accessibility and gene expression in single nuclei at an ultra-high-throughput scale, accommodating up to millions of cells and hundreds of samples. The technique refines the two-step combinatorial indexing approach, initially introduced by Datlinger et al. for snRNA-seq, adapting it specifically for a multiomic context.

In the SUM-seq process, accessible chromatin and nuclear RNA are first tagged with a sample-specific index through transposition and reverse transcription, respectively. Following this, a second index is added using a droplet-based microfluidic platform, such as 10X Chromium. This dual indexing approach allows the overloading of microfluidic droplets while maintaining the ability to assign matched RNA- and ATAC-seq reads to the same individual cell.
<div align="center">
  <img src="https://drive.google.com/uc?export=view&id=1F8vMIbyUR42zOxa378jNyo4ll95jZBCW" alt="image.jpg" />
  <p><i>Introduction to combinatorial indexing</i></p>
</div>

## Key Methods and Findings

### Identification and Visualization of Droplet Differences

- **Knee Plots:** 
  - Created knee plots with log scales on both axes to depict UMI/barcode counts.
  - **Knee Point:** Intersection of linear approximations separating non-empty droplets.
  - **Inflection Point:** Marks the change in the rate of decrease, also used for droplet selection.
  - Each point represented total coverage per barcode, colored based on the proportion of UMIs from a specific sample.
  - For comparison of two samples observed higher UMI counts for UMI enriched barcodes predominantly from the first sample with higher coverage, indicating discrepancies compared to unfiltered and empty droplets.

<div align="center">
  <img src="https://drive.google.com/uc?export=view&id=1q9LQIqvAJZS_kUgajG5j5SJknzYie8YB" alt="image.jpg" />
  <p><i>UMI uneven distribution impacts two-sample based EmptyDrops classification of non-empty/empty droplets</i></p>
</div>

### Analysis with Eight Samples

- **Droplet Categorization:** 
  - Categorized droplets based on UMI counts into three groups: above the knee, between the knee and inflection point, and below the minimum threshold of 100 UMI.
  - Found approximately 70% of theoretically non-empty droplets passed the EmptyDrops (ED) test without significant differences between samples.
  - Noted significant variability in the intermediate zone, with some samples having few non-empty droplets, while others had few unselected by ED.

<div align="center">
  <img src="https://drive.google.com/uc?export=view&id=1seXhuik9YlcFv7jlQI8iG_GMGFkiC2B7" alt="image.jpg" />
  <p><i>Inconsistency in EmptyDrops filtration in different samples occur mostly between inflection and knee points</i></p>
</div>

### Fisher’s Test Analysis

- **Odds Ratio Examination:** 
  - Analyzed the odds ratio for droplets being selected or not across samples using Fisher’s test.
  - Found minor differences in overall droplet selection for low coverage droplets and dramatic differences in the intermediate range, suggesting substantial, yet unexplained, sample-specific variations.

<div align="center">
  <img src="https://drive.google.com/uc?export=view&id=1Bmmt_ZMkSWp6xurwOFrWj_V-ic72r1t5" alt="image.jpg" />
  <p><i>Number of barcodes within inflection-knee range varies dramatically between samples</i></p>
</div>

### Empty Droplet Expression Profile Analysis

- **Similarity Metrics:**
  - Used Spearman correlation, Euclidean distances, and LogProbabilities from EmptyDrops to assess similarity between barcodes from different samples and ambient RNA expression profiles.
  - Addressed the challenge of gene basis comparison using the *highly_variable_genes* function in Scanpy.
  - Found that EmptyDrops' confidence in classifying droplets varied across samples, strongly correlating with sample coverage.
  - Even within a single coverage range, persistent differences between samples were observed, attributed to factors beyond coverage.

<div align="center">
  <img src="https://drive.google.com/uc?export=view&id=1qgrL_-U3kwMQjIg1hz5CzQUF6RH4CATe" alt="image.jpg" />
  <p><i>Number of barcodes within inflection-knee range varies dramatically between samples</i></p>
</div>

## Conclusions

Our findings highlight the complexity of droplet classification in SUMseq data, revealing significant sample-specific variations and the need for more refined analytical approaches to enhance detection accuracy and mitigate erroneous classifications. Our future work will focus on several tasks:

- Analysis of the differences between samples excluding the impact of the sample coverage 
- Analysis of the distribution of the barcodes from several samples regarding inflection and knee point of the samples
- Tuning the parameters of the EmptyDrops analysis regarding the individual sample metrics
 
