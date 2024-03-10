library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
source('EmptyDrops_filtering.R')

#-----------------------------------------------------Function to calculate basic metrics
umi_gene_dist <- function(counts_path, barcodes_path, genes_path, file_name,
                          emptydrops_run = FALSE, lower = 100, 
                          test.ambient = TRUE, fdr_threshold = 0.05, 
                          return_indices = FALSE){
  
  #Read files
  counts <- Matrix::readMM(counts_path)
  genes <- readr::read_tsv(genes_path, col_names = FALSE)
  gene_ids <- genes$X1
  cell_ids <- readr::read_tsv(barcodes_path, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts) <- gene_ids
  colnames(counts) <- cell_ids

  #Calculate UMI per cell barcode distribution
  umi_per_CB <- as.data.frame(Matrix::colSums(counts)[Matrix::colSums(counts) > 0])
  colnames(umi_per_CB) <- c('UMIcount')
  
  #Calculate gene per cell barcode distribution
  gene_per_CB <- as.data.frame(Matrix::colSums(counts > 0)[Matrix::colSums(counts) > 0])
  colnames(gene_per_CB) <- c('Genecount')
  
  #Run emptydrops IF needed
  if (emptydrops_run == TRUE) {
    
    #Data filtering with empty drops
    filtered_counts <- filter_barcodes_with_emptyDrops(matrix = counts, lower = lower, test.ambient = test.ambient)
    
    #Filtered barcodes
    filtered_barcodes_list <- return_filtered_barcodes_or_indices(filtered_counts, fdr_threshold = fdr_threshold, 
                                                                     return_indices = return_indices)
    
    umi_per_CB <- filter(umi_per_CB, rownames(umi_per_CB) %in% filtered_barcodes_list)
    gene_per_CB <- filter(gene_per_CB, rownames(gene_per_CB) %in% filtered_barcodes_list)
  }

  
  #Create plot for UMI per cell barcode distribution
  umi_per_CB_plot <- ggplot(umi_per_CB, aes(x = log10(UMIcount))) +
    geom_histogram(bins=15, fill = 'orange', col = 'darkgreen') +
    labs(y = 'Frequency', x = 'UMI count, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Create plot for gene per cell barcode distribution
  gene_per_CB_plot <- ggplot(gene_per_CB, aes(x = log10(Genecount))) +
    geom_histogram(bins=15, fill = 'darkgreen', col = 'orange') +
    labs(y = 'Frequency', x = 'gene count, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Create DF of two stats
  umi_gene_data <- data.frame(umi_per_CB, gene_per_CB)
  
  #Create scatterplot for UMI per CB ~ gene per CB
  gene_umi <- ggplot(umi_gene_data, aes(x = Genecount, y = UMIcount)) +
    geom_point(color = 'darkgreen') +
    labs(y = 'UMI count', x = 'gene count') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))

  #Calculate mitochondrial genes content
  mito_genes <- c(grep(pattern = c("mt-"), genes$X2, value = TRUE), grep(pattern = c("MT-"), genes$X2, value = TRUE))
  mito_genes_ens <- (filter(genes, X2 %in% mito_genes))$X1
  percent_vector <- c(round(sum(counts[rownames(counts) %in% mito_genes_ens, ]) / sum(counts), 4), 1 - round(sum(counts[rownames(counts) %in% mito_genes_ens, ]) / sum(counts), 4))
  label_vector <- c('MT', 'NonMT')
  
  #Create DF for piechart
  pie_data <- data.frame(
    group = label_vector,
    value = percent_vector,
    labels = percent_vector*100
  )
  
  #Create piechart
  pie_chart <- ggplot(pie_data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    geom_text(aes(label = labels),
              position = position_stack(vjust = 0.5)) +
    coord_polar("y", start=0) +
    theme_void()
  
  #Create sorted UMI_count DF
  umi_per_CB_sorted <- umi_per_CB %>%
    arrange(desc(umi_per_CB))
  
  #Create knee plot
  knee_plot <- ggplot(umi_per_CB_sorted, aes(x = c(1:length(umi_per_CB_sorted$UMIcount)), y = UMIcount)) +
    geom_line(color = 'darkgreen', size = 1) +
    labs(y = 'UMI count, log10', x = 'cell barcodes, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")) +
    scale_x_log10() +
    scale_y_log10()
  
  #Arrange 3 plots into 1 figure
  final_plot <- grid.arrange(
    umi_per_CB_plot, 
    gene_per_CB_plot,
    gene_umi,
    pie_chart,
    knee_plot,
    nrow = 3,
    ncol = 2,
    top = textGrob(file_name, gp=gpar(fontsize=18,font=8))
  )
  
  return (final_plot)
}


#-----------------------------------------------------Function to calculate UMI distribution between samples

sample_dist <- function(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                        counts_path_s2, barcodes_path_s2, genes_path_s2,
                        file_name, emptydrops_run = FALSE, lower = 100, 
                        test.ambient = TRUE, fdr_threshold = 0.05, 
                        return_indices = FALSE){
  
  #----------------------------------------Import Solo1-----------------------------------------------------------------------------------------------
  
  #Read files
  counts_s1 <- Matrix::readMM(counts_path_s1)
  genes_s1 <- readr::read_tsv(genes_path_s1, col_names = FALSE)
  gene_ids_s1 <- genes_s1$X1
  cell_ids_s1 <- readr::read_tsv(barcodes_path_s1, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts_s1) <- gene_ids_s1
  colnames(counts_s1) <- cell_ids_s1
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB_s1_mtx <- Matrix::colSums(counts_s1)[Matrix::colSums(counts_s1) > 0]
  
  #Convert to DF
  DF_umi_per_CB_s1 <- as.data.frame(umi_per_CB_s1_mtx)
  
  #Run emptydrops IF needed
  if (emptydrops_run == TRUE) {
    
    #Data filtering with empty drops
    filtered_counts_s1 <- filter_barcodes_with_emptyDrops(matrix = counts_s1, lower = lower, test.ambient = test.ambient)
    
    #Filtered barcodes
    filtered_barcodes_list_s1 <- return_filtered_barcodes_or_indices(filtered_counts_s1, fdr_threshold = fdr_threshold, 
                                                                  return_indices = return_indices)
    
    DF_umi_per_CB_s1 <- filter(DF_umi_per_CB_s1, rownames(DF_umi_per_CB_s1) %in% filtered_barcodes_list_s1)
  }
  
  #----------------------------------------Import Solo2-----------------------------------------------------------------------------------------------
  
  #Read files
  counts_s2 <- Matrix::readMM(counts_path_s2)
  genes_s2 <- readr::read_tsv(genes_path_s2, col_names = FALSE)
  gene_ids_s2 <- genes_s2$X1
  cell_ids_s2 <- readr::read_tsv(barcodes_path_s2, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts_s2) <- gene_ids_s2
  colnames(counts_s2) <- cell_ids_s2
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB_s2_mtx <- Matrix::colSums(counts_s2)[Matrix::colSums(counts_s2) > 0]
  
  #Convert to DF
  DF_umi_per_CB_s2 <- as.data.frame(umi_per_CB_s2_mtx)
  
  #Run emptydrops IF needed
  if (emptydrops_run == TRUE) {
    
    #Data filtering with empty drops
    filtered_counts_s2 <- filter_barcodes_with_emptyDrops(matrix = counts_s2, lower = lower, test.ambient = test.ambient)
    
    #Filtered barcodes
    filtered_barcodes_list_s2 <- return_filtered_barcodes_or_indices(filtered_counts_s2, fdr_threshold = fdr_threshold, 
                                                                  return_indices = return_indices)
    
    DF_umi_per_CB_s2 <- filter(DF_umi_per_CB_s2, rownames(DF_umi_per_CB_s2) %in% filtered_barcodes_list_s2)
  }
  
  #----------------------------------------Calculate distributions----------------------------------------------------------------------------------------------
  
  #Merge dataset
  DF_umi_per_CB_s1_s2 <- merge(DF_umi_per_CB_s1, DF_umi_per_CB_s2, 
                               by = 'row.names', all = FALSE)
  
  colnames(DF_umi_per_CB_s1_s2) <- c('Barcodes', 'umi_per_CB_s1', 'umi_per_CB_s2')
  
  #Replcase NA's with 0
  DF_umi_per_CB_s1_s2[is.na(DF_umi_per_CB_s1_s2)] <- 0
  
  #Calculate distribution
  DF_umi_per_CB_s1_s2 <- DF_umi_per_CB_s1_s2 %>%
    mutate(dist_s1 = round(umi_per_CB_s1/(umi_per_CB_s1 + umi_per_CB_s2),5)) %>%
    mutate(dist_s2 = round(umi_per_CB_s2/(umi_per_CB_s1 + umi_per_CB_s2),5)) %>%
    mutate(sum = umi_per_CB_s1 + umi_per_CB_s2)
  
  #----------------------------------------Binomial testing----------------------------------------------------------------------------------------------
  
  binomial.test <- function(s1, sum, p = 0.5) {binom.test(s1, sum, p)$p.value}
   
  DF_umi_per_CB_s1_s2$p_val_binom <- mapply(binomial.test, DF_umi_per_CB_s1_s2$umi_per_CB_s1, DF_umi_per_CB_s1_s2$sum)

  DF_umi_per_CB_s1_s2$p_val_binom2 <- dbinom(DF_umi_per_CB_s1_s2$umi_per_CB_s1, DF_umi_per_CB_s1_s2$sum, 0.5)
  
  return(DF_umi_per_CB_s1_s2)
}

#-----------------------------------------------------Function to plot UMI distribution between samples

plot_distributions <- function(distributions_dataset){
  #Plot distribution with 0 and 1 values
  dist_plot <- ggplot(distributions_dataset) +
    geom_density(aes(x = dist_s1), fill = '#28B463', alpha = 0.7) +
    geom_density(aes(x = dist_s2), fill = '#FF7F50', alpha = 0.7) +
    labs(y = 'Counts', x = 'Distribution') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Create function to calculate number of UMI in sample1 for each step of distribution 
  mean_median_max_umi_count_dist_s1 <- function(k){
    mean <- mean(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    median <- median(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    max <- max(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    return(c(mean, median, max))
  }
  
  #Create function to calculate number of UMI in sample2 for each step of distribution 
  mean_median_max_umi_count_dist_s2 <- function(k){
    mean <- mean(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    median <- median(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    max <- max(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    return(c(mean, median, max))
  }
  
  #Initialize steps
  k <- seq(0, 1, by = 0.05)
  
  #Implement above function and calculate number of UMI in sample1 
  mean_median_max_umi_count_with_k_s1 <- as.data.frame(sapply(k, mean_median_max_umi_count_dist_s1))
  colnames(mean_median_max_umi_count_with_k_s1) <- k
  
  #Create dataset for the plot
  mean_median_max_umi_count_with_k_s1 <- mean_median_max_umi_count_with_k_s1 %>%
    #select(where(~!any(is.na(.)))) %>%
    mutate(Param = c('Mean', 'Median', 'Max')) %>%
    tidyr::pivot_longer(!Param)
  
  #Plot with UMI counts sample 1 for each step of distribution
  mean_max_umi_count_dist_plot_s1 <- ggplot(mean_median_max_umi_count_with_k_s1, aes(x = name, y = value, color = Param)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    ylab('UMI count sample1') +
    xlab('Distribution steps') +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  
  #Implement above function and calculate number of UMI in sample2 
  mean_median_max_umi_count_with_k_s2 <- as.data.frame(sapply(k, mean_median_max_umi_count_dist_s2))
  colnames(mean_median_max_umi_count_with_k_s2) <- k
  
  #Create dataset for the plot
  mean_median_max_umi_count_with_k_s2 <- mean_median_max_umi_count_with_k_s2 %>%
    #select(where(~!any(is.na(.)))) %>%
    mutate(Param = c('Mean', 'Median', 'Max')) %>%
    tidyr::pivot_longer(!Param)
  
  #Plot with UMI counts sample2 for each step of distribution
  mean_max_umi_count_dist_plot_s2 <- ggplot(mean_median_max_umi_count_with_k_s2, aes(x = name, y = value, color = Param)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    ylab('UMI count sample2') +
    xlab('Distribution steps') +
    theme(legend.position = "bottom") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  
  #Combine plots
  final_plot <- grid.arrange(
    mean_max_umi_count_dist_plot_s1,
    mean_max_umi_count_dist_plot_s2,
    dist_plot,
    nrow = 2,
    ncol = 2,
    top = textGrob(file_name, gp=gpar(fontsize=18,font=8))
  )

  return(final_plot)
}


#-----------------------------------------------------Solo1 gene raw data-------------------------------------------------------------------------------------------------------------------------


counts_path_s1_g_raw = '../Data/Solo.out_1/Gene/raw/matrix.mtx.gz'
barcodes_path_s1_g_raw = '../Data/Solo.out_1/Gene/raw/barcodes.tsv.gz'
features_path_s1_g_raw = '../Data/Solo.out_1/Gene/raw/features.tsv.gz'
file_name_s1_g_raw = 'Solo1_gene_raw'

plot_s1_g_raw <- umi_gene_dist(counts_path = counts_path_s1_g_raw, barcodes_path = barcodes_path_s1_g_raw,
                               genes_path = features_path_s1_g_raw, file_name = file_name_s1_g_raw, 
                               emptydrops_run = TRUE, fdr_threshold = 0.01)

#---------------------------------------------------- Solo1 gene filtered data-----------------------------------------------------------------------------------------------------------------


counts_path_s1_g_filt = '../Data/Solo.out_1/Gene/filtered/matrix.mtx.gz'
barcodes_path_s1_g_filt = '../Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz'
features_path_s1_g_filt = '../Data/Solo.out_1/Gene/filtered/features.tsv.gz'
file_name_s1_g_filt = 'Solo1_gene_filt'

plot_s1_g_filt <- umi_gene_dist(counts_path = counts_path_s1_g_filt, barcodes_path = barcodes_path_s1_g_filt,
                               genes_path = features_path_s1_g_filt, file_name = file_name_s1_g_filt)

#-----------------------------------------------------Solo1 full gene raw data----------------------------------------------------------------------------------------------------------------------------


counts_path_s1_fg_raw = '../Data/Solo.out_1/GeneFull/raw/matrix.mtx.gz'
barcodes_path_s1_fg_raw = '../Data/Solo.out_1/GeneFull/raw/barcodes.tsv.gz'
features_path_s1_fg_raw = '../Data/Solo.out_1/GeneFull/raw/features.tsv.gz'
file_name_s1_fg_raw = 'Solo1_geneFull_raw'

plot_s1_fg_raw <- umi_gene_dist(counts_path = counts_path_s1_fg_raw, barcodes_path = barcodes_path_s1_fg_raw,
              genes_path = features_path_s1_fg_raw, file_name = file_name_s1_fg_raw)

#-----------------------------------------------------Solo1 full gene filtered data---------------------------------------------------------------------------------------------------------------------------


counts_path_s1_fg_filt = '../Data/Solo.out_1/GeneFull/filtered/matrix.mtx.gz'
barcodes_path_s1_fg_filt = '../Data/Solo.out_1/GeneFull/filtered/barcodes.tsv.gz'
features_path_s1_fg_filt = '../Data/Solo.out_1/GeneFull/filtered/features.tsv.gz'
file_name_s1_fg_filt = 'Solo1_geneFull_filt'

plot_s1_fg_filt <- umi_gene_dist(counts_path = counts_path_s1_fg_filt, barcodes_path = barcodes_path_s1_fg_filt,
                                genes_path = features_path_s1_fg_filt, file_name = file_name_s1_fg_filt)


solo1_out <- grid.arrange(plot_s1_g_raw, 
                          plot_s1_g_filt,
                          plot_s1_fg_raw,
                          plot_s1_fg_filt)


solo1_out


#-----------------------------------------------------Solo2 gene raw data-------------------------------------------------------------------------------------------------------------------------


counts_path_s2_g_raw = '../Data/Solo.out_2/Gene/raw/matrix.mtx.gz'
barcodes_path_s2_g_raw = '../Data/Solo.out_2/Gene/raw/barcodes.tsv.gz'
features_path_s2_g_raw = '../Data/Solo.out_2/Gene/raw/features.tsv.gz'
file_name_s2_g_raw = 'Solo2_gene_raw'

plot_s2_g_raw <- umi_gene_dist(counts_path = counts_path_s2_g_raw, barcodes_path = barcodes_path_s2_g_raw,
                               genes_path = features_path_s2_g_raw, file_name = file_name_s2_g_raw)

#---------------------------------------------------- Solo2 gene filtered data-----------------------------------------------------------------------------------------------------------------


counts_path_s2_g_filt = '../Data/Solo.out_2/Gene/filtered/matrix.mtx.gz'
barcodes_path_s2_g_filt = '../Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz'
features_path_s2_g_filt = '../Data/Solo.out_2/Gene/filtered/features.tsv.gz'
file_name_s2_g_filt = 'Solo2_gene_filt'

plot_s2_g_filt <- umi_gene_dist(counts_path = counts_path_s2_g_filt, barcodes_path = barcodes_path_s2_g_filt,
                                genes_path = features_path_s2_g_filt, file_name = file_name_s2_g_filt)

#-----------------------------------------------------Solo2 full gene raw data----------------------------------------------------------------------------------------------------------------------------


counts_path_s2_fg_raw = '../Data/Solo.out_2/GeneFull/raw/matrix.mtx.gz'
barcodes_path_s2_fg_raw = '../Data/Solo.out_2/GeneFull/raw/barcodes.tsv.gz'
features_path_s2_fg_raw = '../Data/Solo.out_2/GeneFull/raw/features.tsv.gz'
file_name_s2_fg_raw = 'Solo2_geneFull_raw'

plot_s2_fg_raw <- umi_gene_dist(counts_path = counts_path_s2_fg_raw, barcodes_path = barcodes_path_s2_fg_raw,
                                genes_path = features_path_s2_fg_raw, file_name = file_name_s2_fg_raw)

#-----------------------------------------------------Solo2 full gene filtered data---------------------------------------------------------------------------------------------------------------------------


counts_path_s2_fg_filt = '../Data/Solo.out_2/GeneFull/filtered/matrix.mtx.gz'
barcodes_path_s2_fg_filt = '../Data/Solo.out_2/GeneFull/filtered/barcodes.tsv.gz'
features_path_s2_fg_filt = '../Data/Solo.out_2/GeneFull/filtered/features.tsv.gz'
file_name_s2_fg_filt = 'Solo2_geneFull_filt'

plot_s2_fg_filt <- umi_gene_dist(counts_path = counts_path_s2_fg_filt, barcodes_path = barcodes_path_s2_fg_filt,
                                 genes_path = features_path_s2_fg_filt, file_name = file_name_s2_fg_filt)


solo2_out <- grid.arrange(plot_s2_g_raw, 
                          plot_s2_g_filt,
                          plot_s2_fg_raw,
                          plot_s2_fg_filt)


solo2_out

#-----------------------------------------------------geneFull raw data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("../Data/Solo.out_1/GeneFull/raw/matrix.mtx.gz")
genes_path_s1 <- ("../Data/Solo.out_1/GeneFull/raw/features.tsv.gz")
barcodes_path_s1 <- ("../Data/Solo.out_1/GeneFull/raw/barcodes.tsv.gz")

counts_path_s2 <- ("../Data/Solo.out_2/GeneFull/raw/matrix.mtx.gz")
genes_path_s2 <- ("../Data/Solo.out_2/GeneFull/raw/features.tsv.gz")
barcodes_path_s2 <- ("../Data/Solo.out_2/GeneFull/raw/barcodes.tsv.gz")

file_name = 'Gene full Raw data'

genefull_rawdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                                counts_path_s2, barcodes_path_s2, genes_path_s2,
                                file_name, emptydrops_run = TRUE, fdr_threshold = 0.01)

genefull_rawdata_plots <- plot_distributions(genefull_rawdata) 

#-----------------------------------------------------geneFull filt data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("../Data/Solo.out_1/GeneFull/filtered/matrix.mtx.gz")
genes_path_s1 <- ("../Data/Solo.out_1/GeneFull/filtered/features.tsv.gz")
barcodes_path_s1 <- ("../Data/Solo.out_1/GeneFull/filtered/barcodes.tsv.gz")

counts_path_s2 <- ("../Data/Solo.out_2/GeneFull/filtered/matrix.mtx.gz")
genes_path_s2 <- ("../Data/Solo.out_2/GeneFull/filtered/features.tsv.gz")
barcodes_path_s2 <- ("../Data/Solo.out_2/GeneFull/filtered/barcodes.tsv.gz")

file_name = 'Gene full Filt data'

genefull_filtdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                                 counts_path_s2, barcodes_path_s2, genes_path_s2,
                                 file_name)

genefull_filtdata_plots <- plot_distributions(genefull_filtdata) 

#-----------------------------------------------------gene raw data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("../Data/Solo.out_1/Gene/raw/matrix.mtx.gz")
genes_path_s1 <- ("../Data/Solo.out_1/Gene/raw/features.tsv.gz")
barcodes_path_s1 <- ("../Data/Solo.out_1/Gene/raw/barcodes.tsv.gz")

counts_path_s2 <- ("../Data/Solo.out_2/Gene/raw/matrix.mtx.gz")
genes_path_s2 <- ("../Data/Solo.out_2/Gene/raw/features.tsv.gz")
barcodes_path_s2 <- ("../Data/Solo.out_2/Gene/raw/barcodes.tsv.gz")

file_name = 'Gene Raw data'

gene_rawdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                            counts_path_s2, barcodes_path_s2, genes_path_s2,
                            file_name)

gene_rawdata_plots <- plot_distributions(gene_rawdata) 

#-----------------------------------------------------gene filt data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("../Data/Solo.out_1/Gene/filtered/matrix.mtx.gz")
genes_path_s1 <- ("../Data/Solo.out_1/Gene/filtered/features.tsv.gz")
barcodes_path_s1 <- ("../Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz")

counts_path_s2 <- ("../Data/Solo.out_2/Gene/filtered/matrix.mtx.gz")
genes_path_s2 <- ("../Data/Solo.out_2/Gene/filtered/features.tsv.gz")
barcodes_path_s2 <- ("../Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz")

file_name = 'Gene Filt data'

gene_filtdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                             counts_path_s2, barcodes_path_s2, genes_path_s2,
                             file_name)

gene_filtdata_plots <- plot_distributions(gene_filtdata) 

combined_plot <- grid.arrange(
  genefull_rawdata_plots,
  genefull_filtdata_plots,
  gene_rowdata_plots,
  gene_filtdata_plots,
  nrow = 2,
  ncol = 2
)




