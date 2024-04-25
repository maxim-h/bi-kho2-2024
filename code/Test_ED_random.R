source('Basic_metrics.R')

counts <- c('../Data/Solo.out_1/Gene/raw/matrix.mtx.gz',
            '../Data/Solo.out_2/Gene/raw/matrix.mtx.gz')

genes <- c('../Data/Solo.out_1/Gene/raw/features.tsv.gz',
           '../Data/Solo.out_2/Gene/raw/features.tsv.gz')

barcodes <- c('../Data/Solo.out_1/Gene/raw/barcodes.tsv.gz',
              '../Data/Solo.out_2/Gene/raw/barcodes.tsv.gz')

filenames <- c('Solo1', 'Solo2')

common_list = list()

for (i in 1:10){
  experiment <- read_10X_samples(counts, barcodes, genes, 
                                     filenames, emptydrops_run = TRUE, fdr_threshold = 0.05)
  common_list[[i]] <- experiment
}

bc1 <- as.data.frame(common_list[[1]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc2 <- as.data.frame(common_list[[2]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc3 <- as.data.frame(common_list[[3]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc4 <- as.data.frame(common_list[[4]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc5 <- as.data.frame(common_list[[5]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc6 <- as.data.frame(common_list[[6]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc7 <- as.data.frame(common_list[[7]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc8 <- as.data.frame(common_list[[8]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc9 <- as.data.frame(common_list[[9]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')
bc10 <- as.data.frame(common_list[[10]]$data[[1]]) %>%
  tibble::rownames_to_column(var = 'Barcodes')

merged_dataset <- experiments$data[[1]] %>%
  tibble::rownames_to_column(var = 'Barcodes')

for (i in 2:number_of_samples){
  merged_dataset <- merge(merged_dataset, tibble::rownames_to_column(experiments$data[[i]], var = 'Barcodes'), by = 'Barcodes', all = FALSE)
} %>%
  tibble::rownames_to_column(var = 'Barcodes')

bc_merged <- merge(bc1, bc2, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc3, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc4, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc5, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc6, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc7, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc8, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc9, by='Barcodes', all = TRUE)
bc_merged <- merge(bc_merged, bc10, by='Barcodes', all = TRUE)

colnames(bc_merged) <- c('Barcodes', 's1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10')

bc_merged[is.na(bc_merged)] <- 0

test_upset_bc <- bc_merged %>%
  tibble::column_to_rownames(var = 'Barcodes') %>%
  mutate_if(is.numeric, ~1 * (. > 0)) 

intersection_plot <- ComplexUpset::upset(test_upset_bc, colnames(test_upset_bc), name = "Intersection", width_ratio=0.1)
