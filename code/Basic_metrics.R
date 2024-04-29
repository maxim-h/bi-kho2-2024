library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ComplexUpset)
source('EmptyDrops_filtering.R')

#----------------------------------------Function to calculate basic metrics

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
    labs(y = 'UMI count, log10', x = 'Cell barcodes, ranked') +
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


#----------------------------------------Function to calculate UMI distribution between samples

read_10X_samples <- function(counts_path, barcodes_path, genes_path, 
                        file_name, emptydrops_run = FALSE, lower = 100, 
                        test.ambient = TRUE, fdr_threshold = 0.05, 
                        return_indices = FALSE){
  
  #----------------------------------------Import data-----------------------------------------------------------------------------------------------
  
  number_of_samples <- length(counts)
  
  # Create empty list
  experiments <- list(data = list(), names = list(), counts_full = list(), counts = list(), distributions = list(), non_filt_pvals = list(), umi_sum = list())

  
  for (i in 1:number_of_samples){
    count_tables <- Matrix::readMM(counts[i])
    rownames(count_tables) <- readr::read_tsv(genes[i], col_names = FALSE)$X1
    colnames(count_tables) <- readr::read_tsv(barcodes[i], col_names = FALSE)$X1
    counts_data <- as.data.frame(Matrix::colSums(count_tables)[Matrix::colSums(count_tables) > 0])
    
    #Run emptydrops IF needed
    if (emptydrops_run == TRUE) {
      
      #Data filtering with empty drops
      filtered_counts <- filter_barcodes_with_emptyDrops(matrix = count_tables, lower = lower, test.ambient = test.ambient)
      
      #Filtered barcodes
      filtered_barcodes_list <- return_filtered_barcodes_or_indices(filtered_counts, fdr_threshold = fdr_threshold, 
                                                                       return_indices = return_indices)
      
      counts_data <- filter(counts_data, rownames(counts_data) %in% filtered_barcodes_list)
    }
    colnames(counts_data) <- glue::glue('UMI_count_sample_{i}')
    experiments$data[[i]] <- counts_data
    experiments$names[[i]] <- filenames[i]
    experiments$umi_sum[[i]] <- sum(counts_data)
  }
  return(experiments)
}

#----------------------------------------Calculate distributions----------------------------------------------------------------------------------------------
  
calculate_distributions <- function(experiments, number_of_samples){
  
  merged_dataset <- experiments$data[[1]] %>%
    tibble::rownames_to_column(var = 'Barcodes')
  
  for (i in 2:number_of_samples){
    merged_dataset <- merge(merged_dataset, tibble::rownames_to_column(experiments$data[[i]], var = 'Barcodes'), by = 'Barcodes', all = FALSE)
  }

  merged_dataset_full <- experiments$data[[1]] %>%
    tibble::rownames_to_column(var = 'Barcodes')
  
  for (i in 2:number_of_samples){
    merged_dataset_full <- merge(merged_dataset_full, tibble::rownames_to_column(experiments$data[[i]], var = 'Barcodes'), by = 'Barcodes', all = TRUE)
  }
  
  merged_dataset_full[is.na(merged_dataset_full)] <- 0
  
  #Calculate distribution
  merged_dataset_distributions <- merged_dataset %>%
    dplyr::mutate(sum_UMI = rowSums(select(., starts_with('UMI')))) %>%
    dplyr::mutate_at(vars(starts_with("UMI")), function(x) (x / .$sum_UMI))

  experiments$counts_full <- merged_dataset_full
  experiments$counts <- merged_dataset
  experiments$distributions <- merged_dataset_distributions
  
  return(experiments)
}

#----------------------------------------Binomial testing----------------------------------------------------------------------------------------------
  
multinomial_testing <- function(experiments, number_of_samples,  prob = 0.5){
  
  #library(EMT)
    
  #perform_multinomial_test <- function(row_data) {
  #  p_value <- EMT::multinomial.test(row_data, rep(1/number_of_samples, number_of_samples))$p.value
  #  return(p_value)
  #}
  
  # Apply the function row-wise to the selected columns and store p-values
  #merged_dataset_distributions$p_val_multinom <- mapply(perform_multinomial_test, select(merged_dataset, starts_with('UMI_')))

  if (number_of_samples == 2){
    binomial.test <- function(s1, sum, p = 0.5) {binom.test(s1, sum, p)$p.value}
    
    # Binomial test with custom function
    experiments$counts$p_val <- mapply(binomial.test, experiments$counts$UMI_count_sample_1, experiments$counts$UMI_count_sample_1 + experiments$counts$UMI_count_sample_2)
    
    # P-value adjustment for custom binomial test
    experiments$counts$p_val_adj <- p.adjust(experiments$counts$p_val, method='BH', n = length(experiments$counts$UMI_count_sample_1))
  } else {
    probs <- rep(1/number_of_samples, number_of_samples)
    for (i in 1:length(experiments$counts$Barcodes)){
      values <- as.numeric(as.vector(select(experiments$counts, starts_with('UMI_'))[1, ]))
      experiments$counts$p_val[i] <- EMT::multinomial.test(values, probs, MonteCarlo = TRUE)$p.value
    }
    experiments$counts$p_val_adj <- p.adjust(experiments$counts$p_val, method='BH', n = length(experiments$counts$UMI_count_sample_1))
  }
  
  pvals <- select(experiments$counts, c(Barcodes, p_val_adj))
  
  experiments$non_filt_pvals <- pvals
  # Binomial test with dbinom
  #DF_umi_per_CB_s1_s2$p_val_binom2 <- dbinom(DF_umi_per_CB_s1_s2$umi_per_CB_s1, DF_umi_per_CB_s1_s2$sum, 0.5)
  
  # P-value adjustment for binomial test with dbinom
  #DF_umi_per_CB_s1_s2$p_val_adj_binom2 <- p.adjust(DF_umi_per_CB_s1_s2$p_val_binom2, method='BH', n = length(DF_umi_per_CB_s1_s2$p_val_binom2))
  
  return(experiments)
}

#----------------------------------------Assign p-values----------------------------------------------------------------------------------------------

assign_non_filt_pvals <- function(experiments, pvals){
  experiments$counts <- dplyr::inner_join(experiments$counts, pvals, by='Barcodes')
  return(experiments)
}

#----------------------------------------Function to plot UMI distribution between samples

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
          axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5))

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


plot_distributions_density <- function(distributions_dataset){
  
  distributions_densityplot <- distributions_dataset %>%
    dplyr::select(!starts_with('sum')) %>%
    tidyr::pivot_longer(!Barcodes, names_to = 'samples', values_to = 'distribution') %>%
    ggplot(., aes(x = distribution, fill = samples)) +
    geom_density(alpha = 0.7) +
    scale_x_log10() +
    xlab('Distributions') +
    xlab('Density') +
    theme_bw()
  
  return(distributions_densityplot)
}


plot_distributions_boxplots <- function(distributions_dataset, number_of_samples, alpha=0.25){
  
  distributions_boxplot <- distributions_dataset %>%
    tibble::column_to_rownames(var = 'Barcodes') %>%
    select(starts_with('UMI_')) %>%
    tibble::rownames_to_column(var = 'Barcodes') %>%
    tidyr::pivot_longer(!Barcodes, names_to = 'samples', values_to = 'distribution') %>%
    mutate(samples = as.factor(samples)) %>%
    ggplot(aes(x = samples, y = distribution, fill = samples)) +
    geom_point(position = 'jitter', alpha = alpha) +
    geom_boxplot(outlier.alpha = 0) +
    geom_hline(yintercept=1/number_of_samples, linetype = "dashed", color = "darkgreen", size = 1) +
    theme_classic() +
    theme(axis.text.x=element_blank())
  
  return(distributions_boxplot)
}


plot_intersections <- function(counts_dataset_full, min_size=100){
  
  test_upset_data <- counts_dataset_full %>%
    tibble::column_to_rownames(var = 'Barcodes') %>%
    select(starts_with('UMI_')) %>%
    mutate_if(is.numeric, ~1 * (. > 0)) 

  intersection_plot <- ComplexUpset::upset(test_upset_data, colnames(test_upset_data), name = "Intersection", width_ratio=0.1, min_size=min_size)
  
  return(intersection_plot)
}


knee_plot_sample_distribution <- function(experiments, sample, threshold = 0.05){
  
  sample_number <- glue::glue('UMI_count_sample_{sample}')
  
  knee_plot_data <- cbind(experiments$counts[sample_number], experiments$distributions[sample_number], experiments$counts$p_val_adj)
  
  colnames(knee_plot_data) <- c('UMI_count', 'UMI_count_distribution', 'P_val_adj')
  
  knee_plot_data <- knee_plot_data %>%
    arrange(desc(UMI_count)) %>%
    mutate(FDR_sign = ifelse(P_val_adj < threshold, 1, -1))
           
  knee_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = UMI_count, color = UMI_count_distribution)) +
    geom_point(size=3) +
    theme_bw() +
    ylab('UMI distribution') +
    #scale_color_continuous(type = "viridis") +
    scale_color_gradientn(colors = c('#fde725','#5ec962', '#21918c', '#3b528b')) +
    scale_x_log10(breaks=c(0, 10**round(log10(nrow(knee_plot_data))))) +
    scale_y_log10() +
    ggtitle(sample_number) +
    theme(legend.position = "None",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  bar_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = FDR_sign, color = UMI_count_distribution)) +
    geom_bar(stat = 'identity') +
    geom_hline(yintercept = 0) + 
    theme_bw() +
    xlab('Barcodes, ranked by sum') +
    scale_color_gradientn(colors = c('#fde725','#5ec962', '#21918c', '#3b528b')) +
    #scale_color_continuous(type='viridis') +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    scale_x_log10(breaks=c(0,10**round(log10(nrow(knee_plot_data)))))
  
  final_plot <- ggarrange(knee_plot,
            bar_plot, 
            nrow=2,
            align='v',
            heights = c(2,1))

  return(final_plot)
}



