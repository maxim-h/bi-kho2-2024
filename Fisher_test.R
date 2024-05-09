filtered_by_ED <- c(2737,	2265,	1848,	1146,	2484,	2037,	1644,	610,	14771)
empty_by_ED <- c(431,	376,	1,	40,	563,	440,	250,	689,	2790)


ED_on_BR <- rbind.data.frame(filtered_by_ED, empty_by_ED)
data
colnames(ED_on_BR) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8", "Summ_Sample")
rownames(ED_on_BR) <- c('filtered_by_ED', 'empty_by_ED')
#test for one
#data1 <- rbind.data.frame(ED_on_BR$Sample4, ED_on_BR$Summ_Sample)
#result <- fisher.test(data1)
#print(result)

results = list()
for (i in 1:(ncol(ED_on_BR) - 1)) {  
  sample_data <- rbind.data.frame(ED_on_BR[,i], ED_on_BR$Summ_Sample)
  fisher_result <- fisher.test(sample_data)
  results[[i]] <- data.frame(
    sample = colnames(ED_on_BR)[i],
    odds_ratio = fisher_result$estimate,
    lower_ci = fisher_result$conf.int[1],
    upper_ci = fisher_result$conf.int[2]
  )
} 
results_df <- do.call(rbind, results)
  
fisher_plot <- ggplot(results_df, aes(x = sample, y = odds_ratio)) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs( y = "Odds Ratio", title = "Odds Ratio for each Sample vs. Sum") +
  scale_y_continuous(oob = scales::squish_infinite, 
                     limits = c(NA, 10))+
  scale_y_log10()
  theme_bw()

fisher_plot


###
#for Solo1 и Solo2

ED_yes <- c(768,	397)
ED_no <- c(120,	270)
ED_on_BR_Solo12 <- cbind.data.frame(ED_yes, ED_no)

rownames(ED_on_BR_Solo12) <- c("Solo1","Solo2")
colnames(ED_on_BR_Solo12) <- c('filtered_by_ED', 'empty_by_ED')
result <- fisher.test(ED_on_BR_Solo12)
# odds_ratio = 4.348







