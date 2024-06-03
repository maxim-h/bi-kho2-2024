library(DropletUtils)
library(glue)

filter_barcodes_with_emptyDrops <- function(matrix, lower = 100, test.ambient = TRUE) {
  emptyDrops_df <- emptyDrops(m = matrix, lower = lower, test.ambient = test.ambient)
  emptyDrops_df$FDR[is.na(emptyDrops_df$FDR)] <- 1 # выставляю всем NA FDR единичку
  return(emptyDrops_df)
}

return_filtered_barcodes_or_indices <- function(emptyDrops_df, fdr_threshold = 0.05, return_indices = FALSE) {
  filtered_barcodes_indices <- which(emptyDrops_df$FDR < fdr_threshold) # список со всеми отфильтр. баркодами
  if (return_indices) {
    return(filtered_barcodes_indices)
  }
  return(rownames(emptyDrops_df)[filtered_barcodes_indices])
}
