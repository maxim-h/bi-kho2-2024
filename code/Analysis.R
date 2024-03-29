source('Basic_metrics.R')

#-----------------------------------------------------Knee plot and one sample statistic-------------------------------------------------------------------------------------------------------------------------


counts_path <- counts[1]
barcodes_path <- barcodes[1]
genes_path <- genes[1]
file_name <- filenames[1]

one_sample_plot <- umi_gene_dist(counts_path = counts_path, barcodes_path = barcodes_path,
                               genes_path = genes_path, file_name = file_name, 
                               emptydrops_run = TRUE, fdr_threshold = 0.05)

#-----------------------------------------------------UMI distribution between samples-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- '../Data/Solo.out_1/Gene/raw/matrix.mtx.gz'
genes_path_s1 <- '../Data/Solo.out_1/Gene/raw/features.tsv.gz'
barcodes_path_s1 <- '../Data/Solo.out_1/Gene/raw/barcodes.tsv.gz'

counts_path_s2 <- '../Data/Solo.out_2/Gene/raw/matrix.mtx.gz'
genes_path_s2 <- '../Data/Solo.out_2/Gene/raw/features.tsv.gz'
barcodes_path_s2 <- '../Data/Solo.out_2/Gene/raw/barcodes.tsv.gz'

file_name = 'FILENAME'

sample_distribution_dataset <- sample_dist(counts, barcodes, genes, 
                                filenames, emptydrops_run = TRUE, fdr_threshold = 0.01)

sample_distribution_plot <- plot_distributions(sample_distribution_dataset) 


counts <- c('../Data/3.STAR/Fixed_fresh_K562_noPEG_R1/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_fresh_K562_noPEG_R2/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_fresh_K562_PEG_R1/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_fresh_K562_PEG_R2/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_frozen_K562_noPEG_R1/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_frozen_K562_noPEG_R2/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_frozen_K562_PEG_R1/Solo.out/GeneFull/raw/matrix.mtx',
               '../Data/3.STAR/Fixed_frozen_K562_PEG_R2/Solo.out/GeneFull/raw/matrix.mtx'
               )

genes <- c('/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_noPEG_R1/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_noPEG_R2/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_PEG_R1/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_PEG_R2/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_noPEG_R1/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_noPEG_R2/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_PEG_R1/Solo.out/GeneFull/raw/features.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_PEG_R2/Solo.out/GeneFull/raw/features.tsv'
              )

barcodes <- c('/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_noPEG_R1/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_noPEG_R2/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_PEG_R1/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_fresh_K562_PEG_R2/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_noPEG_R1/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_noPEG_R2/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_PEG_R1/Solo.out/GeneFull/raw/barcodes.tsv',
              '/Users/suleymanov-ef/Desktop/Additional education/Bioinformatics Institute/Bioinformatics/Projects/Empty_drops/bi-kho2-2024/Data/3.STAR/Fixed_frozen_K562_PEG_R2/Solo.out/GeneFull/raw/barcodes.tsv'
              )

filenames <- c('Fixed_fresh_K562_noPEG_R1',
               'Fixed_fresh_K562_noPEG_R2',
               'Fixed_fresh_K562_PEG_R1',
               'Fixed_fresh_K562_PEG_R2',
               'Fixed_frozen_K562_noPEG_R1',
               'Fixed_frozen_K562_noPEG_R2',
               'Fixed_frozen_K562_PEG_R1',
               'Fixed_frozen_K562_PEG_R2'
               )