library(Matrix)
#---------------------------------------Read 10X experiment

sc_matrix_10X <- function(counts_path, barcodes_path, genes_path, file_names) {
  count_tables <- Matrix::readMM(counts_path)
  rownames(count_tables) <- readr::read_tsv(genes_path, col_names = FALSE)$X1
  colnames(count_tables) <- readr::read_tsv(barcodes_path, col_names = FALSE)$X1

  return(count_tables)
}

calculate_umi_per_cb <- function(sc_matrix) {
  # Calculate UMI per cell barcode distribution
  umi_per_CB <- as.data.frame(Matrix::colSums(sc_matrix)[Matrix::colSums(sc_matrix) > 0])
  colnames(umi_per_CB) <- c("UMIcount")
  return(umi_per_CB)
}

calculate_gene_per_cb <- function(sc_matrix) {
  # Calculate gene per cell barcode distribution
  gene_per_CB <- as.data.frame(Matrix::colSums(sc_matrix > 0)[Matrix::colSums(sc_matrix) > 0])
  colnames(gene_per_CB) <- c("Genecount")
  return(gene_per_CB)
}
