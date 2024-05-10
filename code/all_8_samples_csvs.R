library(DropletUtils)
library(glue)
library(MatrixExtra)
library(ggplot2)
###Загрузка данных
base_dir <- '/home/vgrigoriants/projects/ib_project/data/Species_mixing'
samples_dirs <- c(glue('{base_dir}/Solo.out_1/GeneFull/raw'),
                  glue('{base_dir}/Solo.out_2/GeneFull/raw'),
                  glue('{base_dir}/Fixed_fresh_K562_PEG_R1/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_frozen_K562_noPEG_R2/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_fresh_K562_noPEG_R1/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_fresh_K562_PEG_R2/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_frozen_K562_PEG_R1/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_fresh_K562_noPEG_R2/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_frozen_K562_noPEG_R1/Solo.out/GeneFull/raw'),
                  glue('{base_dir}/Fixed_frozen_K562_PEG_R2/Solo.out/GeneFull/raw'))
names(samples_dirs) <- c('Solo.out_1', 'Solo.out_2', 'Fixed_fresh_K562_PEG_R1', 'Fixed_frozen_K562_noPEG_R2',
'Fixed_fresh_K562_noPEG_R1', 'Fixed_fresh_K562_PEG_R2', 'Fixed_frozen_K562_PEG_R1',
'Fixed_fresh_K562_noPEG_R2', 'Fixed_frozen_K562_noPEG_R1', 'Fixed_frozen_K562_PEG_R2')


######Запуск ED и сохранение всех данных в csv файлах
for (i in seq_along(samples_dirs)) {
  element_name <- names(samples_dirs)[i]
  element_value <- samples_dirs[i]
  sce_object <- read10xCounts(element_value, type = 'sparse')
  sce_object_emptyDrops_df <- emptyDrops(sce_object, test.ambient = TRUE, lower = 100)
  sce_object_emptyDrops_df$Barcodes <- sce_object@colData$Barcode
  sce_object_emptyDrops_df <- sce_object_emptyDrops_df[, c('Barcodes', 'Total', 'LogProb')]
  sce_object_emptyDrops_df <- na.omit(sce_object_emptyDrops_df)
  sce_object_emptyDrops_df$Sample <- element_name
  write.csv(sce_object_emptyDrops_df, glue('{base_dir}/csvs/{element_name}.csv'), row.names = FALSE)
}



#####Считывание csv файлов и отрисовка распрежедений LogProb для всех сэмплов
ED_Fixed_fresh_K562_PEG_R1 <- read.csv(glue('{base_dir}/csvs/Fixed_fresh_K562_PEG_R1.csv'))
ED_logs_dataset_8s <- data.frame(ED_Fixed_fresh_K562_PEG_R1)
for (i in seq_along(samples_dirs)[4:length(samples_dirs)]) {
  print(i)
  element_name <- names(samples_dirs)[i]
  df <- read.csv(glue('{base_dir}/csvs/{element_name}.csv'))
  ED_logs_dataset_8s <- rbind(ED_logs_dataset_8s, df)
}
ED_logs_dataset_8s$Grid <- 1
ED_logs_dataset_8s[ED_logs_dataset_8s$Sample %in% c('Fixed_fresh_K562_PEG_R1',
                                                    'Fixed_fresh_K562_PEG_R2'), ]$Grid <- 2
ED_logs_dataset_8s[ED_logs_dataset_8s$Sample %in% c('Fixed_frozen_K562_noPEG_R1',
                                                    'Fixed_frozen_K562_noPEG_R2'), ]$Grid <- 3
ED_logs_dataset_8s[ED_logs_dataset_8s$Sample %in% c('Fixed_frozen_K562_PEG_R1',
                                                    'Fixed_frozen_K562_PEG_R2'), ]$Grid <- 4
ED_logs_dataset_8s$Grid <- as.factor(ED_logs_dataset_8s$Grid)
ggplot(ED_logs_dataset_8s, aes(x = LogProb)) +
  geom_density(aes(fill = as.factor(Sample)), alpha = .5) +
  xlim(-200, 0) +
  # coord_cartesian(xlim=c(-100, 0)) +
  # scale_y_continuous(trans='log10') +
  # ylim(0, 0.05) +
  labs(title = "Distribution of LogProbs of all samples", 
       x = "log_Prob") +
  facet_wrap(~Grid) +
  theme_minimal()
