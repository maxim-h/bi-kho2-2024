library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(DropletUtils)
library(glue)
#функция возвращает DF, полученный из отсортированного по  количеству UMI в двух(!) образцах DF,
# из которого были выбраны только СВ, встречающиеся в двух образцах

UMI_per_CB_sample_dist <- function(counts_path_s1, barcodes_path_s1, genes_path_s1, counts_path_s2, barcodes_path_s2, genes_path_s2){

  #Read files
  counts <- Matrix::readMM(counts_path_s1)
  genes <- readr::read_tsv(genes_path_s1, col_names = FALSE)
  gene_ids <- genes$X1
  cell_ids <- readr::read_tsv(barcodes_path_s1, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts) <- gene_ids
  colnames(counts) <- cell_ids
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB <- Matrix::colSums(counts)[Matrix::colSums(counts) > 0]
  
  #Calculate gene per cell barcode distribution
  gene_per_CB <- Matrix::colSums(counts > 0)[Matrix::colSums(counts) > 0]
  
  
  #Create sorted UMI_count DF
  DF_umi_per_CB1 <- as.data.frame(umi_per_CB) %>%
    arrange(desc(umi_per_CB)) %>%
    mutate(cb = as.numeric(rownames(as.data.frame(umi_per_CB))))
  
  #Read files
  counts <- Matrix::readMM(counts_path_s2)
  genes <- readr::read_tsv(genes_path_s2, col_names = FALSE)
  gene_ids <- genes$X1
  cell_ids <- readr::read_tsv(barcodes_path_s2, col_names = FALSE)$X1
  
  #Assign row and colnames
  rownames(counts) <- gene_ids
  colnames(counts) <- cell_ids
  
  #Calculate UMI per cell barcode distribution
  umi_per_CB <- Matrix::colSums(counts)[Matrix::colSums(counts) > 0]
  
  #Calculate gene per cell barcode distribution
  gene_per_CB <- Matrix::colSums(counts > 0)[Matrix::colSums(counts) > 0]
  
  
  #Create sorted UMI_count DF
  DF_umi_per_CB2 <- as.data.frame(umi_per_CB) %>%
    arrange(desc(umi_per_CB)) %>%
    mutate(cb = as.numeric(rownames(as.data.frame(umi_per_CB))))
  
  
  DF_umi_per_CB_s1_s2_sort <- merge(DF_umi_per_CB1, DF_umi_per_CB2, 
                               by = 'row.names', all = TRUE)
  #DF_only_common <- merge(DF_umi_per_CB1, DF_umi_per_CB2, 
  #                by = 'row.names', all = FALSE)
  #Replcase NA's with 0
  DF_umi_per_CB_s1_s2_sort[is.na(DF_umi_per_CB_s1_s2_sort)] <- 0
  
  #Calculate distribution
  DF_umi_per_CB_s1_s2_sort <- DF_umi_per_CB_s1_s2_sort %>%
    mutate(dist_s1 = round(umi_per_CB.x/(umi_per_CB.x + umi_per_CB.y),3)) %>%
    mutate(dist_s2 = round(umi_per_CB.y/(umi_per_CB.x + umi_per_CB.y),3))
  
  #DF_for_knee_s1 <-subset(DF_umi_per_CB_s1_s2_sort, dist_s1 > 0)
  #DF_for_knee_s2 <-subset(DF_umi_per_CB_s1_s2_sort, dist_s2 > 0)
  
  DF_all <- DF_umi_per_CB_s1_s2_sort %>%
    mutate(dist_s1 = round(umi_per_CB.x/(umi_per_CB.x + umi_per_CB.y),3)) %>%
    mutate(dist_s2 = round(umi_per_CB.y/(umi_per_CB.x + umi_per_CB.y),3))%>%
    mutate(umi_per_CB_total = umi_per_CB.y + umi_per_CB.x)
  DF_all_sorted <- as.data.frame(DF_all) %>%
    arrange(desc(DF_all$umi_per_CB_total))
  DF_all_sorted <- as.data.frame(DF_all_sorted) %>%
    mutate(common_CB_number = as.numeric(rownames(as.data.frame(DF_all_sorted))))
  
  DF_only_common_sorted <-subset(DF_all_sorted, dist_s1 > 0)
  DF_only_common_sorted <-subset(DF_only_common_sorted, dist_s1 < 1)
  return (DF_only_common_sorted)
}

#принимает на вход df, c колонками
# common_CB_number номер баркода  из df отсортированного для knee_plot для общих CB
# umi_per_CB количество UMI для этго баркода в  ДВУХ образце 1
# dist_s1  - отношение количества UMI для этого баркода в образце 1 к количеству всех UMI для этого баркода
# и имя образца
knee_plot_UMI_sample_dist <- function(df, name){
  colors <- brewer.pal(n = 11, name = "PRGn")
  knee_plot_UMI_sample_dist <- ggplot(df, aes(x = df$common_CB_number, y = df$umi_per_CB_total, colour = df$dist_s1)) +
    geom_point(size = 1) +
    labs(y = 'UMI count, log10', x = 'cell barcodes, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_gradientn(colours = colors,
                           guide = guide_colorbar(title = "UMI in Sample1 presense"))+
    ggtitle(name)
  return (knee_plot_UMI_sample_dist)
}

#----------------------Empty_drops---------
## Функции для фильтрование баркодов из сырой матрицы при помощи EmptyDrops
filter_barcodes_with_emptyDrops <- function(sparse_matrix, lower = 100, test.ambient = TRUE) {
  emptyDrops_df = emptyDrops(sparse_matrix, lower = lower, test.ambient = test.ambient)
  emptyDrops_df$FDR[is.na(emptyDrops_df$FDR)] <- 1 # выставляю всем NA FDR единичку
  return(emptyDrops_df)
}

return_filtered_barcodes_or_indeces <- function(emptyDrops_df, fdr_threshold = 0.05, return_indeces = FALSE) {
  filtered_barcodes_indeces <- which(emptyDrops_df$FDR < fdr_threshold) #список со всеми отфильтр. баркодами
  
  if (return_indeces) {
    return(filtered_barcodes_indeces)
  }
  return(rownames(emptyDrops_df)[filtered_barcodes_indeces])
}

#Пример фильтрации

#Чтение данных
path_to_matrix <- "Data/Solo.out_1/Gene/raw/matrix.mtx.gz"
path_to_barcodes <- "Data/Solo.out_1/Gene/raw/barcodes.tsv.gz"
matrix_1_raw <- readMM(path_to_matrix)
cell_ids <- readr::read_tsv(path_to_barcodes, col_names = FALSE)$X1
colnames(matrix_1_raw) <- cell_ids

#Запуск emptyDrops
emptyDrops_df <- filter_barcodes_with_emptyDrops(matrix_1_raw)


#Вывод баркодов/индексов в отдельные списки
filtered_barcodes <- return_filtered_barcodes_or_indeces(emptyDrops_df, fdr_threshold = 0.05) #список со всеми отфильтр. баркодами
filtered_barcodes_indeces <- return_filtered_barcodes_or_indeces(emptyDrops_df, fdr_threshold = 0.05, return_indeces = TRUE) #список со всеми индексами отфильтр. баркодов


#------------------------------------DF with common CB for RAW and Filtered Data

#-----------------------------------------------------Solo1-filtered

counts_path_s1 <- ("Data/Solo.out_1/Gene/filtered/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz")

#-----------------------------------------------------Solo2-filtered
counts_path_s2 <- ("Data/Solo.out_2/Gene/filtered/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz")
file_name_s2 <- 'S2_gene_filtered_data'


dfGF <- UMI_per_CB_sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, counts_path_s2, barcodes_path_s2, genes_path_s2)


#-----------------------------------------------------Solo1_raw

counts_path_s1 <- ("Data/Solo.out_1/Gene/raw/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/Gene/raw/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/Gene/raw/barcodes.tsv.gz")

#-----------------------------------------------------Solo2_raw
counts_path_s2 <- ("Data/Solo.out_2/Gene/raw/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/Gene/raw/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/Gene/raw/barcodes.tsv.gz")

dfGR <- UMI_per_CB_sample_dist(counts_path_s1, barcodes_path_s1, genes_path_s1, counts_path_s2, barcodes_path_s2, genes_path_s2)

#-----------------------Merge DF with common BC with EmptyDrops filtration Data
ED_result <- subset(dfGR, Row.names %in% filtered_barcodes)
ED_result_filtered <- subset(dfGF, Row.names %in% filtered_barcodes)
dif_ED_result_filtered <- subset(dfGF, !(Row.names %in% filtered_barcodes))

#write.csv(ED_result, file = "common_gene_row_ED_filtered.csv", row.names = TRUE)
#write.csv(ED_result_filtered, file = "common_gene_filtered_ED_filtered.csv", row.names = TRUE)
#write.csv(dif_ED_result_filtered, file = "common_gene_filtered_not_in_ED_filered.csv", row.names = TRUE)

#------------------Knee_plots for Raw and two different filtrations
plotGR <- knee_plot_UMI_sample_dist(dfGR, 'Solo1(Gene_Raw)')
plotGF <- knee_plot_UMI_sample_dist(dfGF, 'Solo1(Gene_Filtered')

final_plot <- grid.arrange(
  plotGF, 
  plotGR,
  nrow = 2,
  ncol = 1,
  top = textGrob('UMI/CB distribution in 2 samples, colored by presence in Solo1 ', gp=gpar(fontsize=18,font=8))
)

final_plot


plotGR <- knee_plot_UMI_sample_dist(dfGR, 'Solo1(Gene_Raw)')
plotGF <- knee_plot_UMI_sample_dist(ED_result, 'Solo1(Gene_Filtered_byED)')
plotGF
final_plot <- grid.arrange(
  plotGF, 
  plotGR,
  nrow = 2,
  ncol = 1,
  top = textGrob('UMI/CB distribution in 2 samples, colored by presence in Solo1 ', gp=gpar(fontsize=18,font=8))
)

final_plot
