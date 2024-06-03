library(ggplot2)
library(ggpubr)
library(dplyr)


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

plot_distributions <- function(distributions_dataset) {
  # Plot distribution with 0 and 1 values
  dist_plot <- ggplot(distributions_dataset) +
    geom_density(aes(x = dist_s1), fill = "#28B463", alpha = 0.7) +
    geom_density(aes(x = dist_s2), fill = "#FF7F50", alpha = 0.7) +
    labs(y = "Counts", x = "Distribution") +
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )

  # Create function to calculate number of UMI in sample1 for each step of distribution
  mean_median_max_umi_count_dist_s1 <- function(k) {
    mean <- mean(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    median <- median(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    max <- max(filter(distributions_dataset, dist_s1 == k)$umi_per_CB_s1)
    return(c(mean, median, max))
  }

  # Create function to calculate number of UMI in sample2 for each step of distribution
  mean_median_max_umi_count_dist_s2 <- function(k) {
    mean <- mean(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    median <- median(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    max <- max(filter(distributions_dataset, dist_s2 == k)$umi_per_CB_s2)
    return(c(mean, median, max))
  }

  # Initialize steps
  k <- seq(0, 1, by = 0.05)

  # Implement above function and calculate number of UMI in sample1
  mean_median_max_umi_count_with_k_s1 <- as.data.frame(sapply(k, mean_median_max_umi_count_dist_s1))
  colnames(mean_median_max_umi_count_with_k_s1) <- k

  # Create dataset for the plot
  mean_median_max_umi_count_with_k_s1 <- mean_median_max_umi_count_with_k_s1 %>%
    # select(where(~!any(is.na(.)))) %>%
    mutate(Param = c("Mean", "Median", "Max")) %>%
    tidyr::pivot_longer(!Param)

  # Plot with UMI counts sample 1 for each step of distribution
  mean_max_umi_count_dist_plot_s1 <- ggplot(mean_median_max_umi_count_with_k_s1, aes(x = name, y = value, color = Param)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    ylab("UMI count sample1") +
    xlab("Distribution steps") +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )

  # Implement above function and calculate number of UMI in sample2
  mean_median_max_umi_count_with_k_s2 <- as.data.frame(sapply(k, mean_median_max_umi_count_dist_s2))
  colnames(mean_median_max_umi_count_with_k_s2) <- k

  # Create dataset for the plot
  mean_median_max_umi_count_with_k_s2 <- mean_median_max_umi_count_with_k_s2 %>%
    # select(where(~!any(is.na(.)))) %>%
    mutate(Param = c("Mean", "Median", "Max")) %>%
    tidyr::pivot_longer(!Param)

  # Plot with UMI counts sample2 for each step of distribution
  mean_max_umi_count_dist_plot_s2 <- ggplot(mean_median_max_umi_count_with_k_s2, aes(x = name, y = value, color = Param)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    ylab("UMI count sample2") +
    xlab("Distribution steps") +
    theme(legend.position = "bottom") +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )

  # Combine plots
  final_plot <- grid.arrange(
    mean_max_umi_count_dist_plot_s1,
    mean_max_umi_count_dist_plot_s2,
    dist_plot,
    nrow = 2,
    ncol = 2,
    top = textGrob(file_name, gp = gpar(fontsize = 18, font = 8))
  )

  return(final_plot)
}

plot_distributions_density <- function(experiments) {
  distributions_dataset <- experiments$distributions

  distributions_densityplot <- distributions_dataset %>%
    dplyr::select(!starts_with("sum")) %>%
    tidyr::pivot_longer(!Barcodes, names_to = "samples", values_to = "distribution") %>%
    ggplot(., aes(x = distribution, fill = samples)) +
    geom_density(alpha = 0.7) +
    scale_x_log10() +
    xlab("Distributions") +
    xlab("Density") +
    theme_bw()

  return(distributions_densityplot)
}

plot_distributions_boxplots <- function(experiments, number_of_samples, alpha = 0.25) {
  distributions_dataset <- experiments$distributions

  distributions_boxplot <- distributions_dataset %>%
    tibble::column_to_rownames(var = "Barcodes") %>%
    select(starts_with("UMI_")) %>%
    tibble::rownames_to_column(var = "Barcodes") %>%
    tidyr::pivot_longer(!Barcodes, names_to = "samples", values_to = "distribution") %>%
    mutate(samples = as.factor(samples)) %>%
    ggplot(aes(x = samples, y = distribution, fill = samples)) +
    geom_point(position = "jitter", alpha = alpha) +
    geom_boxplot(outlier.alpha = 0) +
    geom_hline(yintercept = 1 / number_of_samples, linetype = "dashed", color = "darkgreen", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_blank())

  return(distributions_boxplot)
}

plot_intersections <- function(experiments, min_intersection_size = 100) {
  test_upset_data <- experiments$counts_full %>%
    tibble::column_to_rownames(var = "Barcodes") %>%
    select(starts_with("UMI_")) %>%
    mutate_if(is.numeric, ~ 1 * (. > 0))

  intersection_plot <- ComplexUpset::upset(test_upset_data, colnames(test_upset_data), name = "Intersection", width_ratio = 0.1, min_size = min_intersection_size)

  return(intersection_plot)
}

knee_plot_sample_distribution <- function(experiments, sample, threshold = 0.05) {
  sample_number <- glue::glue("UMI_count_sample_{sample}")

  knee_plot_data <- cbind(experiments$counts[sample_number], experiments$distributions[sample_number], experiments$counts$p_val_adj)

  colnames(knee_plot_data) <- c("UMI_count", "UMI_count_distribution", "P_val_adj")

  knee_plot_data <- knee_plot_data %>%
    arrange(desc(UMI_count)) %>%
    mutate(FDR_sign = ifelse(P_val_adj < threshold, 1, -1))

  knee_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = UMI_count, color = UMI_count_distribution)) +
    geom_point(size = 3) +
    theme_bw() +
    ylab(glue::glue("UMI counts in {sample} sample")) +
    # scale_color_continuous(type = "viridis") +
    scale_color_gradientn(colors = c("#fde725", "#5ec962", "#21918c", "#3b528b")) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(knee_plot_data))))) +
    scale_y_log10() +
    ggtitle(sample_number) +
    theme(
      legend.position = "None",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  bar_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = FDR_sign, color = UMI_count_distribution)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    xlab("Barcodes, ranked by sum") +
    ylab("Multinomial test\nfalse/true") +
    scale_color_gradientn(colors = c("#fde725", "#5ec962", "#21918c", "#3b528b")) +
    # scale_color_continuous(type='viridis') +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(knee_plot_data))))) +
    scale_y_continuous(breaks = c(-1, 1))

  final_plot <- ggarrange(knee_plot,
    bar_plot,
    nrow = 2,
    align = "v",
    heights = c(2, 1)
  )

  return(final_plot)
}

knee_plot_for_sum_of_sample_wth_sample_distribution <- function(experiments, sample, threshold = 0.05) {
  sample_number <- glue::glue("UMI_count_sample_{sample}")

  knee_plot_data <- cbind(as.data.frame(rowSums(select(experiments$counts, starts_with("UMI_count")))), experiments$distributions[sample_number], experiments$counts$p_val_adj)

  colnames(knee_plot_data) <- c("UMI_count", "UMI_count_distribution", "P_val_adj")

  knee_plot_data <- knee_plot_data %>%
    arrange(desc(UMI_count)) %>%
    mutate(FDR_sign = ifelse(P_val_adj < threshold, 1, -1))

  knee_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = UMI_count, color = UMI_count_distribution)) +
    geom_point(size = 3) +
    theme_bw() +
    ylab("UMI counts in total sample") +
    # scale_color_continuous(type = "viridis") +
    scale_color_gradientn(colors = c("#fde725", "#5ec962", "#21918c", "#3b528b")) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(knee_plot_data))))) +
    scale_y_log10() +
    ggtitle(sample_number) +
    theme(
      legend.position = "None",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  bar_plot <- ggplot(knee_plot_data, aes(x = c(1:nrow(knee_plot_data)), y = FDR_sign, color = UMI_count_distribution)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0) +
    theme_bw() +
    xlab("Barcodes, ranked by sum") +
    ylab("Multinomial test\nfalse/true") +
    scale_color_gradientn(colors = c("#fde725", "#5ec962", "#21918c", "#3b528b")) +
    # scale_color_continuous(type='viridis') +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_x_log10(breaks = c(0, 10**round(log10(nrow(knee_plot_data))))) +
    scale_y_continuous(breaks = c(-1, 1))

  final_plot <- ggarrange(knee_plot,
    bar_plot,
    nrow = 2,
    align = "v",
    heights = c(2, 1)
  )

  return(final_plot)
}
