library(ggplot2)
library(dplyr)

#---------------------------------------Calculate distances between cells and emptyDrop profile

calculate_distances <- function(sc_matrix, emptyDrops_profile, cor_method = "spearman", euclidean = FALSE, returnDF = TRUE) {
  distances_intermediate_cell_vs_empty_drops <- list(barcodes = list(), distance_value = list())

  for (i in 1:sc_matrix@Dim[2]) {
    distances_intermediate_cell_vs_empty_drops$barcodes[i] <- sc_matrix@Dimnames[[2]][i]
    if (euclidean) {
      euclidean_dist <- function(a, b) sqrt(sum((a - b)^2))
      distances_intermediate_cell_vs_empty_drops$distance_value[i] <- euclidean_dist(as.vector(sc_matrix[, i]), emptyDrops_profile)
    } else {
      distances_intermediate_cell_vs_empty_drops$distance_value[i] <- (1 - abs(cor(as.vector(sc_matrix[, i]), emptyDrops_profile, method = cor_method)))
    }
  }
  distances_intermediate_cell_vs_empty_drops_df <- as.data.frame(x = unlist(distances_intermediate_cell_vs_empty_drops$distance_value), row.names = unlist(distances_intermediate_cell_vs_empty_drops$barcodes))

  if (returnDF) {
    return(distances_intermediate_cell_vs_empty_drops_df)
  } else {
    return(distances_intermediate_cell_vs_empty_drops)
  }
}


#---------------------------------------Knee plot with distances between cells and emptyDrop profile
knee_plot_with_distances <- function(dataset, barcodes_passed_ED) {
  counts_data_merged <- dataset %>%
    mutate(ED_empty = ifelse(Barcodes %in% barcodes_passed_ED, -1, 1)) %>%
    arrange(desc(UMI_count))

  knee_plot <- ggplot(counts_data_merged, aes(x = c(1:nrow(counts_data_merged)), y = UMI_count, color = Distance)) +
    geom_point(size = 3) +
    theme_bw() +
    ylab("UMI distribution") +
    # scale_color_continuous(type = "viridis") +
    scale_colour_gradientn(
      colours = c("#3b528b", "#21918c", "#5ec962", "#f6b092", "#fde725", "#f6b092", "#5ec962", "#21918c", "#3b528b"),
      # breaks = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0),
      # limits = c(-1.0, 1.0),
      values = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
    ) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(counts_data_merged))))) +
    scale_y_log10() +
    # ggtitle(sample_number) +
    theme(
      legend.position = "None",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  bar_plot1 <- ggplot(counts_data_merged, aes(x = c(1:nrow(counts_data_merged)), y = Distance, color = Distance)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = -1, linetype = "dashed") +
    ylim(c(-1, 1)) +
    theme_bw() +
    ylab("Distance") +
    scale_colour_gradientn(
      colours = c("#3b528b", "#21918c", "#5ec962", "#f6b092", "#fde725", "#f6b092", "#5ec962", "#21918c", "#3b528b"),
      # breaks = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0),
      # limits = c(-1.0, 1.0),
      values = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
    ) +
    # scale_color_continuous(type='viridis') +
    theme(
      legend.position = "None",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(counts_data_merged)))))

  bar_plot2 <- ggplot(counts_data_merged, aes(x = c(1:nrow(counts_data_merged)), y = ED_empty, color = Distance)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0) +
    ylim(c(-1, 1)) +
    theme_bw() +
    xlab("Barcodes, ranked by sum") +
    ylab("Empty FALSE/TRUE") +
    scale_colour_gradientn(
      colours = c("#3b528b", "#21918c", "#5ec962", "#f6b092", "#fde725", "#f6b092", "#5ec962", "#21918c", "#3b528b"),
      # breaks = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0),
      # limits = c(-1.0, 1.0),
      values = c(-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
    ) +
    # scale_color_continuous(type='viridis') +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(counts_data_merged)))))

  final_plot <- ggarrange(knee_plot,
    bar_plot1,
    bar_plot2,
    nrow = 3,
    align = "v",
    heights = c(3, 1, 1)
  )
  return(final_plot)
}

#---------------------------------------Merge distances with UMI_counts by Barcodes
merge_barcodes_counts_and_distances <- function(sc_matrix, distances) {
  counts_data <- as.data.frame(Matrix::colSums(sc_matrix))
  counts_data_merged <- merge(counts_data, distances, by = "row.names")
  colnames(counts_data_merged) <- c("Barcodes", "UMI_count", "Distance")
  return(counts_data_merged)
}
