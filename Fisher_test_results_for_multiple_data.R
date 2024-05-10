library(patchwork)
library(ggplot2)

# пишем функцию для теста Фишера


fisher_test_results <- function(data) {
  results <- list()
  
  for (i in 1:(ncol(data) - 1)) {
    sample_data <- rbind.data.frame(data[, i], data$Summ_Sample)
    fisher_result <- fisher.test(sample_data)
    results[[i]] <- data.frame(
      sample = colnames(data)[i],
      odds_ratio = fisher_result$estimate,
      lower_ci = fisher_result$conf.int[1],
      upper_ci = fisher_result$conf.int[2],
      p_value = fisher_result$p.value
    )
  }
    
    return(results)
  }
  
  
  # создаем датафреймы для анализа
  
  # все капли выше c числом UMI > inflection point
  filtered_by_ED <- c(2737,	2265,	1848,	1146,	2484,	2037,	1644,	610,	14771)
  empty_by_ED <- c(431,	376,	1,	40,	563,	440,	250,	689,	2790)
  
  
  data_inf <- rbind.data.frame(filtered_by_ED, empty_by_ED)
  
  colnames(data_inf) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8", "Summ_Sample")
  rownames(data_inf) <- c('filtered_by_ED', 'empty_by_ED')
  
  # капли с inf < nUMI < knee 
  
  filtered_by_ED <- c(176,	68,	267,	326,	23,	13,	375,	4,	1252)
  empty_by_ED <-c(431,	376,	1,	40,	563,	440,	250,	689,	2790)
  data_kn_inf <- rbind.data.frame(filtered_by_ED, empty_by_ED)
  data
  colnames(data_kn_inf) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8", "Summ_Sample")
  rownames(data_kn_inf) <- c('filtered_by_ED', 'empty_by_ED')
  
  # капли с 100 < nUMI < inflection
  filtered_by_ED <- c(3342,	2408,	10692,	6061,	1929,	775,	5658,	1543,	32408)
  empty_by_ED <-c(1976,	1743,	3272,	2166,	1221,	862,	1517,	1138,	13895)
  data_inf_100 <- rbind.data.frame(filtered_by_ED, empty_by_ED)
  data
  colnames(data_inf_100) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8", "Summ_Sample")
  rownames(data_inf_100) <- c('filtered_by_ED', 'empty_by_ED')
  
  # все капли > 100 UMI
  all <- c(8486,	6792,	15813,	9413,	6197,	4114,	9069,	3980,	63864)
  filtered_by_ED <- c(6079,	4673,	12540,	7207,	4413,	2812,	7302,	2153,	47179)
  empty_by_ED <- all - filtered_by_ED
  data_100 <- rbind.data.frame(filtered_by_ED, empty_by_ED)
  data
  colnames(data_100) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8", "Summ_Sample")
  rownames(data_100) <- c('filtered_by_ED', 'empty_by_ED')
  
  
  fisher_100 <- fisher_test_results(data_100)
  fisher_inf <- fisher_test_results(data_inf)
  fisher_inf_100 <- fisher_test_results(data_inf_100)
  fisher_kn_inf <- fisher_test_results(data_kn_inf)
  
  fisher_100 
  fisher_inf 
  fisher_inf_100 
  fisher_kn_inf 
  
  