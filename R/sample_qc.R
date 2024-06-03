library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)


quality_metrics_one_sample <- function(sc_matrix, genes_path, file_name = "Sample", filter = FALSE, filtered_barcodes_list = NULL) {
  # Calculate UMI per cell barcode distribution
  umi_per_CB <- as.data.frame(Matrix::colSums(sc_matrix)[Matrix::colSums(sc_matrix) > 0])
  colnames(umi_per_CB) <- c("UMIcount")

  # Calculate gene per cell barcode distribution
  gene_per_CB <- as.data.frame(Matrix::colSums(sc_matrix > 0)[Matrix::colSums(sc_matrix) > 0])
  colnames(gene_per_CB) <- c("Genecount")

  if (filter == TRUE) {
    umi_per_CB <- filter(umi_per_CB, rownames(umi_per_CB) %in% filtered_barcodes_list)
    gene_per_CB <- filter(gene_per_CB, rownames(gene_per_CB) %in% filtered_barcodes_list)
  }

  genes <- readr::read_tsv(genes_path, col_names = FALSE)
  # Calculate mitochondrial genes content
  mito_genes <- c(grep(pattern = c("mt-"), genes$X2, value = TRUE), grep(pattern = c("MT-"), genes$X2, value = TRUE))
  mito_genes_ens <- (filter(genes, X2 %in% mito_genes))$X1
  percent_vector <- c(round(sum(sc_matrix[rownames(sc_matrix) %in% mito_genes_ens, ]) / sum(sc_matrix), 4), 1 - round(sum(sc_matrix[rownames(sc_matrix) %in% mito_genes_ens, ]) / sum(sc_matrix), 4))
  label_vector <- c("MT", "NonMT")

  # Create DF for piechart
  mito_data <- data.frame(
    group = label_vector,
    value = percent_vector,
    labels = percent_vector * 100
  )

  # Create plot for UMI per cell barcode distribution
  umi_per_CB_plot <- ggplot(umi_per_CB, aes(x = log10(UMIcount))) +
    geom_histogram(bins = 15, fill = "orange", col = "darkgreen") +
    labs(y = "Frequency", x = "UMI count, log10") +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )

  # Create plot for gene per cell barcode distribution
  gene_per_CB_plot <- ggplot(gene_per_CB, aes(x = log10(Genecount))) +
    geom_histogram(bins = 15, fill = "darkgreen", col = "orange") +
    labs(y = "Frequency", x = "gene count, log10") +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )

  # Create DF of two stats
  umi_gene_data <- data.frame(umi_per_CB, gene_per_CB)

  # Create scatterplot for UMI per CB ~ gene per CB
  gene_umi <- ggplot(umi_gene_data, aes(x = Genecount, y = UMIcount)) +
    geom_point(color = "darkgreen") +
    labs(y = "UMI count", x = "gene count") +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )

  # Create piechart
  pie_chart <- ggplot(mito_data, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    geom_text(aes(label = labels),
      position = position_stack(vjust = 0.5)
    ) +
    coord_polar("y", start = 0) +
    theme_void()

  # Create sorted UMI_count DF
  umi_per_CB_sorted <- umi_per_CB %>%
    arrange(desc(umi_per_CB))

  # Create knee plot
  knee_plot <- ggplot(umi_per_CB_sorted, aes(x = c(1:length(umi_per_CB_sorted$UMIcount)), y = UMIcount)) +
    geom_line(color = "darkgreen", size = 1) +
    labs(y = "UMI count, log10", x = "Cell barcodes, ranked") +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    ) +
    scale_x_log10() +
    scale_y_log10()

  # Arrange 3 plots into 1 figure
  final_plot <- grid.arrange(
    umi_per_CB_plot,
    gene_per_CB_plot,
    gene_umi,
    pie_chart,
    knee_plot,
    nrow = 3,
    ncol = 2,
    top = textGrob(file_name, gp = gpar(fontsize = 18, font = 8))
  )

  return(final_plot)
}
