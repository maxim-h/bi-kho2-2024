library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)


odds_ratio_analysis <- function(data, my_name) {
  results <- list()
  
  for (i in 1:(ncol(data))) {
    sample_data <- rbind.data.frame(data[, i], rowSums(data) - data[, i])
    fisher_result <- fisher.test(sample_data)
    results[[i]] <- data.frame(
      sample = as.character(i),
      odds_ratio = fisher_result$estimate,
      lower_ci = fisher_result$conf.int[1],
      upper_ci = fisher_result$conf.int[2]
    )
  }
  results_df <- do.call(rbind, results)
  
  fisher_plot <- ggplot(results_df, aes(x = sample, y = odds_ratio, color = factor(sample))) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), size = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
    labs(y = "Odds Ratio", title = my_name) +
    scale_y_continuous(oob = scales::squish_infinite, limits = c(0.01, 10000)) +
    scale_y_log10(limits = c(0.002, 16500)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      legend.position = "none")   
  
  return(fisher_plot)
}

# select data fo Fisher_test (filtered and not filtered by ED in different UMI counts interval)
sum_data <- read.csv('../Data/8_samples_summarized_data.csv')

# all barcodes with  UMI > inflection point
filtered_by_ED <- sum_data$count_between_knee_and_inflection_ED_yes + sum_data$count_more_than_knee
empty_by_ED <- sum_data$count_between_knee_and_inflection_ED_no
data_inf <- rbind.data.frame(filtered_by_ED, empty_by_ED)

# barcodes with inf < nUMI < knee 
filtered_by_ED <- sum_data$count_between_knee_and_inflection_ED_yes
empty_by_ED <- sum_data$count_between_knee_and_inflection_ED_no
data_kn_inf <- rbind.data.frame(filtered_by_ED, empty_by_ED)

# barcodes with 100 < nUMI < inflection
filtered_by_ED <- sum_data$count_between_100_and_inflection_ED_yes 
empty_by_ED <- sum_data$count_between_100_and_inflection_ED_no
data_inf_100 <- rbind.data.frame(filtered_by_ED, empty_by_ED)


# barcodes > 100 UMI
filtered_by_ED <- sum_data$count_between_100_and_inflection_ED_yes + sum_data$count_between_100_and_inflection_ED_yes
                                                                  + sum_data$count_more_than_knee 
empty_by_ED <- sum_data$count_between_100_and_inflection_ED_no + sum_data$count_between_100_and_inflection_ED_no 
data_100 <- rbind.data.frame(filtered_by_ED, empty_by_ED)

fisher_plot_100 <- odds_ratio_analysis(data_100, 'nUMIs > 100')
fisher_plot_inf <- odds_ratio_analysis(data_inf, 'nUMIs > inflection point')
fisher_plot_inf_100 <- odds_ratio_analysis(data_inf_100, '100 < nUMIs < inflection point' )
fisher_plot_kn_inf <- odds_ratio_analysis(data_kn_inf, 'inflection_point < nUMIs < knee point')


combined_plot <- (fisher_plot_100 |  fisher_plot_inf_100 |  fisher_plot_kn_inf | fisher_plot_inf) + 
  plot_annotation(
    title = "Odds ratio for each sample vs. summarized sample",
    theme = theme(plot.title = element_text(hjust = 0.5)) 
  )


print(combined_plot)



