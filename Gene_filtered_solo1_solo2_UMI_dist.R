library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(viridis)
library(RColorBrewer)
#-----------------------------------------------------Solo1

counts_path_s1 <- ("Data/Solo.out_1/Gene/filtered/matrix.mtx.gz")
genes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/features.tsv.gz")
barcodes_path_s1 <- ("Data/Solo.out_1/Gene/filtered/barcodes.tsv.gz")


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

#-----------------------------------------------------Solo2
counts_path_s2 <- ("Data/Solo.out_2/Gene/filtered/matrix.mtx.gz")
genes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/features.tsv.gz")
barcodes_path_s2 <- ("Data/Solo.out_2/Gene/filtered/barcodes.tsv.gz")
file_name_s2 <- 'S2_gene_filtered_data'



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



#Create DF of two stats
umi_gene_data2 <- data.frame(umi_per_CB, gene_per_CB)


#Create sorted UMI_count DF
DF_umi_per_CB2 <- as.data.frame(umi_per_CB) %>%
  arrange(desc(umi_per_CB)) %>%
  mutate(cb = as.numeric(rownames(as.data.frame(umi_per_CB))))





#-----------------------Knee plot + S1/S2 distribution

DF_umi_per_CB_s1_s2_sort <- merge(DF_umi_per_CB1, DF_umi_per_CB2, 
                             by = 'row.names', all = TRUE)
#Replcase NA's with 0
DF_umi_per_CB_s1_s2_sort[is.na(DF_umi_per_CB_s1_s2_sort)] <- 0

#Calculate distribution
DF_umi_per_CB_s1_s2_sort <- DF_umi_per_CB_s1_s2_sort %>%
  mutate(dist_s1 = round(umi_per_CB.x/(umi_per_CB.x + umi_per_CB.y),2)) %>%
  mutate(dist_s2 = round(umi_per_CB.y/(umi_per_CB.x + umi_per_CB.y),2))

DF_for_knee_s1 <-subset(DF_umi_per_CB_s1_s2_sort, dist_s1 > 0)
DF_for_knee_s2 <-subset(DF_umi_per_CB_s1_s2_sort, dist_s2 > 0)
DF_for_knee_s1_without1 <-subset(DF_for_knee_s1, dist_s1 < 1)
DF_for_knee_s2_without1 <-subset(DF_for_knee_s2, dist_s2 < 1)
#DF_for_knee_s1$dist_s1 <- factor(DF_for_knee_s1$dist_s1) 
#DF_for_knee_s1$dist_s2 <- factor(DF_for_knee_s1$dist_s2)
#DF_for_knee_s2$dist_s1 <- factor(DF_for_knee_s2$dist_s1) 
#DF_for_knee_s2$dist_s2 <- factor(DF_for_knee_s2$dist_s2)

#---------------------------knee_plot_solo1
colors <- brewer.pal(n = 11, name = "PRGn") 

knee_plot_solo1 <- ggplot(DF_for_knee_s1_without1, aes(x = cb.x, y = umi_per_CB.x, colour = dist_s1)) +
  geom_point(size = 1) +
  labs(y = 'UMI count, log10', x = 'cell barcodes, log10') +
  theme_linedraw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_gradientn(colours = colors,
                         guide = guide_colorbar(title = "UMI in Sample1 presense"))+
  ggtitle('Solo1')

#---------------------------knee_plot_solo2

knee_plot_solo2 <- ggplot(DF_for_knee_s2_without1, aes(x = cb.y, y = umi_per_CB.y, colour = dist_s1)) +
  geom_point(size = 1) +
  labs(y = 'UMI count, log10', x = 'cell barcodes, log10') +
  theme_linedraw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_gradientn(colours = colors,
        guide = guide_colorbar(title = "UMI in Sample1 presense"))+
  ggtitle('Solo2')



#Combine plots
final_plot <- grid.arrange(
  knee_plot_solo1, 
  knee_plot_solo2,
  nrow = 2,
  ncol = 1,
  top = textGrob('Gene Filtered Data', gp=gpar(fontsize=18,font=8))
)

final_plot
