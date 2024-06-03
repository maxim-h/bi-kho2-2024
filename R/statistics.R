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


