library(glue)
library(dplyr)
#source('emptyDrops_filtering.R')

#----------------------------------------Function to calculate UMI distribution between samples

read_10X_samples <- function(counts_path_list, barcodes_path_list, genes_path_list, 
                             file_name_list, emptydrops_run = FALSE, lower = 100, 
                             test.ambient = TRUE, fdr_threshold = 0.05, 
                             return_indices = FALSE){
  
  #----------------------------------------Import data-----------------------------------------------------------------------------------------------
  
  number_of_samples <- length(counts_path_list)
  
  # Create empty list
  experiments <- list(data = list(), names = list(), counts_full = list(), counts = list(), distributions = list(), non_filt_pvals = list(), umi_sum = list())
  
  
  for (i in 1:number_of_samples){
    print(glue::glue('Start reading sample {i}'))
    count_tables <- Matrix::readMM(counts_path_list[i])
    rownames(count_tables) <- readr::read_tsv(genes_path_list[i], col_names = FALSE)$X1
    colnames(count_tables) <- readr::read_tsv(barcodes_path_list[i], col_names = FALSE)$X1
    counts_data <- as.data.frame(Matrix::colSums(count_tables)[Matrix::colSums(count_tables) > 0])
    
    #Run emptydrops IF needed
    if (emptydrops_run == TRUE) {
      print(glue::glue('Start filtering sample {i} with EmptyDrops'))
      #Data filtering with empty drops
      filtered_counts <- filter_barcodes_with_emptyDrops(matrix = count_tables, lower = lower, test.ambient = test.ambient)
      
      #Filtered barcodes
      filtered_barcodes_list <- return_filtered_barcodes_or_indices(filtered_counts, fdr_threshold = fdr_threshold, 
                                                                    return_indices = return_indices)
      
      counts_data <- filter(counts_data, rownames(counts_data) %in% filtered_barcodes_list)
      print(glue::glue('Sample {i} has been filtered with EmptyDrops'))
    }
    colnames(counts_data) <- glue::glue('UMI_count_sample_{i}')
    experiments$data[[i]] <- counts_data
    experiments$names[[i]] <- filenames[i]
    experiments$umi_sum[[i]] <- sum(counts_data)
    print(glue::glue('Sample {i} finished'))
    print(glue::glue('====================='))
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

multinomial_testing <- function(experiments, number_of_samples,  probs = rep(1/number_of_samples, number_of_samples)){
  
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
    for (i in 1:length(experiments$counts$Barcodes)){
      values <- as.numeric(as.vector(select(experiments$counts, starts_with('UMI_'))[1, ]))
      #experiments$counts$p_val[i] <- EMT::multinomial.test(values, probs, MonteCarlo = TRUE)$p.value
      experiments$counts$p_val[i] <- chisq.test(values, probs)$p.value
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

#----------------------------------------UMI proportions by samples----------------------------------------------------------------------------------------------

calculate_umi_proportion_by_samples <- function(experiments){
  return(unlist(experiments$umi_sum) / sum(unlist(experiments$umi_sum)))
}