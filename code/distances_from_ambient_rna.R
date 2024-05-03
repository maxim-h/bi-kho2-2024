library(DropletUtils)
library(glue)
library(MatrixExtra)
library(ggplot2)
###Загрузка данных
base_dir <- '/home/vgrigoriants/projects/ib_project/data/Species_mixing'
samples_dirs <- c(glue('{base_dir}/Solo.out_1/GeneFull/raw'),
                  glue('{base_dir}/Solo.out_2/GeneFull/raw'))
names(samples_dirs) <- c('1_raw', '2_raw')
sce_1 <- read10xCounts(samples_dirs[1], type = 'sparse')
sce_2 <- read10xCounts(samples_dirs[2], type = 'sparse')
###Запуск ED
sce_1_emptyDrops_100_df <- emptyDrops(sce_1, test.ambient = TRUE, lower = 100)
sce_2_emptyDrops_100_df <- emptyDrops(sce_2, test.ambient = TRUE, lower = 100)
###Отбор вероятностей и фильтрация по FDR
subset_1 <- sce_1_emptyDrops_100_df[!(is.na(sce_1_emptyDrops_100_df$LogProb)), ][c('LogProb', "FDR")]
subset_1 <- na.omit(subset_1)
subset_2 <- sce_2_emptyDrops_100_df[!(is.na(sce_2_emptyDrops_100_df$LogProb)), ][c('LogProb', "FDR")]
subset_2 <- na.omit(subset_2)
filtered_1 <- subset_1[subset_1$FDR < 0.05, ]
filtered_2 <- subset_2[subset_2$FDR < 0.05, ]
###Отрисовка распределений для сэмплов 1 и 2 
combined_df <- data.frame(value = c(filtered_1$LogProb, filtered_2$LogProb),
                          variable = c(rep("Sample 1", nrow(filtered_1)), rep("Sample 2", nrow(filtered_2))))
colnames(combined_df) <- c('value', 'Sample')
ggplot(combined_df, aes(x = value, fill = Sample)) +
  geom_density(alpha = 0.5) +
  xlim(-13000, 0) +
  labs(title = "Distribution of filtered cells LogProbs of two samples", 
       x = "log_Prob") +
  theme_minimal()
res <- wilcox.test(value ~ Sample, data = combined_df, exact = FALSE)

