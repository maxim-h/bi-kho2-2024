 
# This is the repository for cell calling analysis in SUMseq.


1. <ins>*Count_barcodes.R*</ins> - small script to create barplots with the number of barcodes in different samples before and after filtration.

2. <ins>*EmptyDrops_filtering.R*</ins> - functions to filter empty droplets from a single-cell matrix with custom parameters. The output of the function is a list of barcodes that passed the filtration.

3. <ins>*Basic_metrics.R*</ins> - functions to create pictures with basic metrics for individual sample (UMI per barcodes, gene per barcode, knee-plor, % of mitochondrial genes) and UMI distribution between samples for SUM-seq. Uses EmptyDrops_filtering.R to use filtered datasets for calculations.

4. <ins>*Analysis.R*</ins> - final script to implement previously described functions and create pictures.


