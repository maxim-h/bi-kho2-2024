library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)

#-----------------------------------------------------Function to calculate basic metrics
umi_gene_dist <- function(counts_path, barcodes_path, genes_path, file_name){
  
  #Read files
  counts <- Matrix::readMM(counts_path)
  genes <- readr::read_tsv(genes_path, col_names = FALSE)
  gene_ids <- genes$X1
  cell_ids <- readr::read_tsv(barcodes_path, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts) <- gene_ids
  colnames(counts) <- cell_ids
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB <- Matrix::colSums(counts)[Matrix::colSums(counts) > 0]
  
  #Calculate gene per cell barcode distribution
  gene_per_CB <- Matrix::colSums(counts > 0)[Matrix::colSums(counts) > 0]
  
  #Create plot for UMI per cell barcode distribution
  umi_per_CB_plot <- ggplot(as.data.frame(umi_per_CB), aes(x = log10(umi_per_CB))) +
    geom_histogram(bins=15, fill = 'orange', col = 'darkgreen') +
    labs(y = 'Frequency', x = 'UMI count, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Create plot for gene per cell barcode distribution
  gene_per_CB_plot <- ggplot(as.data.frame(gene_per_CB), aes(x = log10(gene_per_CB))) +
    geom_histogram(bins=15, fill = 'darkgreen', col = 'orange') +
    labs(y = 'Frequency', x = 'gene count, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Create DF of two stats
  umi_gene_data <- data.frame(umi_per_CB, gene_per_CB)
  
  #Create scatterplot for UMI per CB ~ gene per CB
  gene_umi <- ggplot(umi_gene_data, aes(x = gene_per_CB, y = umi_per_CB)) +
    geom_point(color = 'darkgreen') +
    labs(y = 'UMI count', x = 'gene count') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))

  #Calculate mitochondrial genes content
  mito_genes <- grep(pattern = "MT-", genes$X2, value = TRUE)
  percent_vector <- c(sum(counts[rownames(counts) %in% mito_genes, ]) / sum(counts), 1 - (sum(counts[rownames(counts) %in% mito_genes, ]) / sum(counts)))
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
  DF_umi_per_CB <- as.data.frame(umi_per_CB) %>%
    arrange(desc(umi_per_CB)) %>%
    mutate(cb = as.numeric(rownames(as.data.frame(umi_per_CB))))
  
  #Create knee plot
  knee_plot <- ggplot(DF_umi_per_CB, aes(x = cb, y = umi_per_CB)) +
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
                        file_name){
  
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
  umi_per_CB_s1 <- Matrix::colSums(counts_s1)[Matrix::colSums(counts_s1) > 0]
  
  #Convert to DF
  DF_umi_per_CB_s1 <- as.data.frame(umi_per_CB_s1)
  
  #----------------------------------------Import Solo2-----------------------------------------------------------------------------------------------
  
  #Read files
  counts_s2 <- Matrix::readMM("Data/Solo.out_2/GeneFull/raw/matrix.mtx.gz")
  genes_s2 <- readr::read_tsv("Data/Solo.out_2/GeneFull/raw/features.tsv.gz", col_names = FALSE)
  gene_ids_s2 <- genes_s2$X1
  cell_ids_s2 <- readr::read_tsv("Data/Solo.out_2/GeneFull/raw/barcodes.tsv.gz", col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts_s2) <- gene_ids_s2
  colnames(counts_s2) <- cell_ids_s2
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB_s2 <- Matrix::colSums(counts_s2)[Matrix::colSums(counts_s2) > 0]
  
  #Convert to DF
  DF_umi_per_CB_s2 <- as.data.frame(umi_per_CB_s2)
  
  #----------------------------------------Create distributions----------------------------------------------------------------------------------------------
  
  #Merge dataset
  DF_umi_per_CB_s1_s2 <- merge(DF_umi_per_CB_s1, DF_umi_per_CB_s2, 
                               by = 'row.names', all = TRUE)
  
  #Replcase NA's with 0
  DF_umi_per_CB_s1_s2[is.na(DF_umi_per_CB_s1_s2)] <- 0
  
  #Calculate distribution
  DF_umi_per_CB_s1_s2 <- DF_umi_per_CB_s1_s2 %>%
    mutate(dist_s1 = round(umi_per_CB_s1/(umi_per_CB_s1 + umi_per_CB_s2),4)) %>%
    mutate(dist_s2 = round(umi_per_CB_s2/(umi_per_CB_s1 + umi_per_CB_s2),4))
  
  #Plot distribution with 0 and 1 values
  dist_plot_with01 <- ggplot(DF_umi_per_CB_s1_s2) +
    geom_density(aes(x = dist_s1), fill = '#28B463', alpha = 0.7) +
    geom_density(aes(x = dist_s2), fill = '#FF7F50', alpha = 0.7) +
    labs(y = 'Counts', x = 'Distribution') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Plot distribution without 0 and 1 values
  dist_plot<- ggplot(DF_umi_per_CB_s1_s2[(DF_umi_per_CB_s1_s2$dist_s1 > 0) & (DF_umi_per_CB_s1_s2$dist_s2 > 0), ]) +
    geom_density(aes(x = dist_s1), fill = '#28B463', alpha = 0.7) +
    geom_density(aes(x = dist_s2), fill = '#FF7F50', alpha = 0.7) +
    labs(y = 'Counts', x = 'Distribution') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold"))
  
  #Combine plots
  final_plot <- grid.arrange(
    dist_plot, 
    dist_plot_with01,
    nrow = 2,
    ncol = 1,
    top = textGrob(file_name, gp=gpar(fontsize=18,font=8))
  )
  
  return(final_plot)
}

#-----------------------------------------------------Solo1 gene raw data-------------------------------------------------------------------------------------------------------------------------


counts_path_s1_g_raw = 'Data/Solo.out_1/Gene/raw/matrix.mtx.gz'
barcodes_path_s1_g_raw = 'Data/Solo.out_1/Gene/raw/barcodes.tsv.gz'
features_path_s1_g_raw = 'Data/Solo.out_1/Gene/raw/features.tsv.gz'
file_name_s1_g_raw = 'Solo1_gene_raw'

plot_s1_g_raw <- umi_gene_dist(counts_path = counts_path_s1_g_raw, barcodes_path = barcodes_path_s1_g_raw,
                                genes_path = features_path_s1_g_raw, file_name = file_name_s1_g_raw)

#---------------------------------------------------- Solo1 gene filtered data-----------------------------------------------------------------------------------------------------------------


counts_path_s1_g_filt = 'Data/Solo.out_1/Gene/filtered/matrix.mtx.gz'
barcodes_path_s1_g_filt = 'Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz'
features_path_s1_g_filt = 'Data/Solo.out_1/Gene/filtered/features.tsv.gz'
file_name_s1_g_filt = 'Solo1_gene_filt'

plot_s1_g_filt <- umi_gene_dist(counts_path = counts_path_s1_g_filt, barcodes_path = barcodes_path_s1_g_filt,
                               genes_path = features_path_s1_g_filt, file_name = file_name_s1_g_filt)

#-----------------------------------------------------Solo1 full gene raw data----------------------------------------------------------------------------------------------------------------------------


counts_path_s1_fg_raw = 'Data/Solo.out_1/GeneFull/raw/matrix.mtx.gz'
barcodes_path_s1_fg_raw = 'Data/Solo.out_1/GeneFull/raw/barcodes.tsv.gz'
features_path_s1_fg_raw = 'Data/Solo.out_1/GeneFull/raw/features.tsv.gz'
file_name_s1_fg_raw = 'Solo1_geneFull_raw'

plot_s1_fg_raw <- umi_gene_dist(counts_path = counts_path_s1_fg_raw, barcodes_path = barcodes_path_s1_fg_raw,
              genes_path = features_path_s1_fg_raw, file_name = file_name_s1_fg_raw)

#-----------------------------------------------------Solo1 full gene filtered data---------------------------------------------------------------------------------------------------------------------------


counts_path_s1_fg_filt = 'Data/Solo.out_1/GeneFull/filtered/matrix.mtx.gz'
barcodes_path_s1_fg_filt = 'Data/Solo.out_1/GeneFull/filtered/barcodes.tsv.gz'
features_path_s1_fg_filt = 'Data/Solo.out_1/GeneFull/filtered/features.tsv.gz'
file_name_s1_fg_filt = 'Solo1_geneFull_filt'

plot_s1_fg_filt <- umi_gene_dist(counts_path = counts_path_s1_fg_filt, barcodes_path = barcodes_path_s1_fg_filt,
                                genes_path = features_path_s1_fg_filt, file_name = file_name_s1_fg_filt)


solo1_out <- grid.arrange(plot_s1_g_raw, 
                          plot_s1_g_filt,
                          plot_s1_fg_raw,
                          plot_s1_fg_filt)


solo1_out


#-----------------------------------------------------Solo2 gene raw data-------------------------------------------------------------------------------------------------------------------------


counts_path_s2_g_raw = 'Data/Solo.out_2/Gene/raw/matrix.mtx.gz'
barcodes_path_s2_g_raw = 'Data/Solo.out_2/Gene/raw/barcodes.tsv.gz'
features_path_s2_g_raw = 'Data/Solo.out_2/Gene/raw/features.tsv.gz'
file_name_s2_g_raw = 'Solo2_gene_raw'

plot_s2_g_raw <- umi_gene_dist(counts_path = counts_path_s2_g_raw, barcodes_path = barcodes_path_s2_g_raw,
                               genes_path = features_path_s2_g_raw, file_name = file_name_s2_g_raw)

#---------------------------------------------------- Solo2 gene filtered data-----------------------------------------------------------------------------------------------------------------


counts_path_s2_g_filt = 'Data/Solo.out_2/Gene/filtered/matrix.mtx.gz'
barcodes_path_s2_g_filt = 'Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz'
features_path_s2_g_filt = 'Data/Solo.out_2/Gene/filtered/features.tsv.gz'
file_name_s2_g_filt = 'Solo2_gene_filt'

plot_s2_g_filt <- umi_gene_dist(counts_path = counts_path_s2_g_filt, barcodes_path = barcodes_path_s2_g_filt,
                                genes_path = features_path_s2_g_filt, file_name = file_name_s2_g_filt)

#-----------------------------------------------------Solo2 full gene raw data----------------------------------------------------------------------------------------------------------------------------


counts_path_s2_fg_raw = 'Data/Solo.out_2/GeneFull/raw/matrix.mtx.gz'
barcodes_path_s2_fg_raw = 'Data/Solo.out_2/GeneFull/raw/barcodes.tsv.gz'
features_path_s2_fg_raw = 'Data/Solo.out_2/GeneFull/raw/features.tsv.gz'
file_name_s2_fg_raw = 'Solo2_geneFull_raw'

plot_s2_fg_raw <- umi_gene_dist(counts_path = counts_path_s2_fg_raw, barcodes_path = barcodes_path_s2_fg_raw,
                                genes_path = features_path_s2_fg_raw, file_name = file_name_s2_fg_raw)

#-----------------------------------------------------Solo2 full gene filtered data---------------------------------------------------------------------------------------------------------------------------


counts_path_s2_fg_filt = 'Data/Solo.out_2/GeneFull/filtered/matrix.mtx.gz'
barcodes_path_s2_fg_filt = 'Data/Solo.out_2/GeneFull/filtered/barcodes.tsv.gz'
features_path_s2_fg_filt = 'Data/Solo.out_2/GeneFull/filtered/features.tsv.gz'
file_name_s2_fg_filt = 'Solo2_geneFull_filt'

plot_s2_fg_filt <- umi_gene_dist(counts_path = counts_path_s2_fg_filt, barcodes_path = barcodes_path_s2_fg_filt,
                                 genes_path = features_path_s2_fg_filt, file_name = file_name_s2_fg_filt)


solo2_out <- grid.arrange(plot_s2_g_raw, 
                          plot_s2_g_filt,
                          plot_s2_fg_raw,
                          plot_s2_fg_filt)


solo2_out

#-----------------------------------------------------Solo1 geneFull raw data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("Data/Solo.out_1/GeneFull/raw/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/GeneFull/raw/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/GeneFull/raw/barcodes.tsv.gz")

counts_path_s2 <- ("Data/Solo.out_2/GeneFull/raw/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/GeneFull/raw/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/GeneFull/raw/barcodes.tsv.gz")

file_name = 'Gene full Raw data'

genefull_rawdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                                counts_path_s2, barcodes_path_s2, genes_path_s2,
                                file_name)


#-----------------------------------------------------Solo1 geneFull filt data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("Data/Solo.out_1/GeneFull/filtered/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/GeneFull/filtered/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/GeneFull/filtered/barcodes.tsv.gz")

counts_path_s2 <- ("Data/Solo.out_2/GeneFull/filtered/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/GeneFull/filtered/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/GeneFull/filtered/barcodes.tsv.gz")

file_name = 'Gene full Filt data'

genefull_filtdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                                 counts_path_s2, barcodes_path_s2, genes_path_s2,
                                 file_name)


#-----------------------------------------------------Solo1 gene raw data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("Data/Solo.out_1/Gene/raw/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/Gene/raw/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/Gene/raw/barcodes.tsv.gz")

counts_path_s2 <- ("Data/Solo.out_2/Gene/raw/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/Gene/raw/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/Gene/raw/barcodes.tsv.gz")

file_name = 'Gene Raw data'

gene_rawdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                            counts_path_s2, barcodes_path_s2, genes_path_s2,
                            file_name)


#-----------------------------------------------------Solo1 gene raw data-------------------------------------------------------------------------------------------------------------------------

counts_path_s1 <- ("Data/Solo.out_1/Gene/filtered/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz")

counts_path_s2 <- ("Data/Solo.out_2/Gene/filtered/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz")

file_name = 'Gene Filt data'

gene_filtdata <- sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, 
                             counts_path_s2, barcodes_path_s2, genes_path_s2,
                             file_name)


combined_plot <- grid.arrange(
  gene_rawdata,
  gene_filtdata,
  genefull_rawdata,
  genefull_filtdata,
  nrow = 2,
  ncol = 2
)

