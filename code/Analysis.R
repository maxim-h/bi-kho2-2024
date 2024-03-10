source('Basic_metrics.R')

#-----------------------------------------------------Knee plot and one sample statistic-------------------------------------------------------------------------------------------------------------------------


counts_path <- 'PATH_TO_SPARSE_MATRIX'
barcodes_path <- 'PATH_TO_BARCODES'
features_path <- 'PATH_TO_GENES'
file_name <- 'FILENAME'

one_sample_plot <- umi_gene_dist(counts_path = counts_path, barcodes_path = barcodes_path,
                               genes_path = features_path, file_name = file_name, 
                               emptydrops_run = TRUE, fdr_threshold = 0.05)

#-----------------------------------------------------UMI distribution between samples-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- 'PATH_TO_SPARSE_MATRIX_S1'
genes_path_s1 <- 'PATH_TO_BARCODES_S1'
barcodes_path_s1 <- 'PATH_TO_GENES_S1'

counts_path_s2 <- 'PATH_TO_SPARSE_MATRIX_S2'
genes_path_s2 <- 'PATH_TO_BARCODES_S2'
barcodes_path_s2 <- 'PATH_TO_GENES_S2'

file_name = 'FILENAME'

genefull_rawdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                                counts_path_s2, barcodes_path_s2, genes_path_s2,
                                file_name, emptydrops_run = TRUE, fdr_threshold = 0.01)

genefull_rawdata_plots <- plot_distributions(genefull_rawdata) 



