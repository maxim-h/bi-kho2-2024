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
subset_1 <- sce_1_emptyDrops_100_df[!(is.na(sce_1_emptyDrops_100_df$LogProb)), ]#[c('LogProb', "FDR")]
subset_1 <- na.omit(subset_1)
subset_2 <- sce_2_emptyDrops_100_df[!(is.na(sce_2_emptyDrops_100_df$LogProb)), ]#[c('LogProb', "FDR")]
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


#############
add_fdrs_and_lobProbs_to_barcodes <- function(sce_object, empty_drops_df) { 
  empty_drops_df$FDR[is.na(empty_drops_df$FDR)] <- 1 #выставляю всем NA FDR единичку
  empty_drops_df$LogProb[is.na(empty_drops_df$LogProb)] <- 0 #выставляю всем NA LogProbs ноль
  sce_object$EmptyDrops_fdr <- empty_drops_df$FDR
  sce_object$EmptyDrops_logProb <- empty_drops_df$LogProb
  return(sce_object)
}

sce_1 <- add_fdrs_and_lobProbs_to_barcodes(sce_1, sce_1_emptyDrops_100_df)
sce_2 <- add_fdrs_and_lobProbs_to_barcodes(sce_2, sce_2_emptyDrops_100_df)
############
plot_barcode_rank_plot <- function(sce_object, 
                                   title='title',
                                   logProb_threshold = FALSE,
                                   logProb_thresholds = c(0,0)) {
  if (logProb_threshold) { #можно красить по результатам empty_drops
    sample_colors <- ifelse((sce_object$EmptyDrops_logProb > logProb_thresholds[1]) & 
                            (sce_object$EmptyDrops_logProb < logProb_thresholds[2]), 
                            adjustcolor('red', alpha.f = 1), 
                            adjustcolor('blue', alpha.f = .2))
  }
  else {
    sce_object$EmptyDrops_logProb[sce_object$EmptyDrops_logProb > -5000] <- -5000
    # sce_object$EmptyDrops_logProb[sce_object$EmptyDrops_logProb < -13000] <- -13000
    color_value <- sce_object$EmptyDrops_logProb
    print(max(color_value))
    print(min(color_value))
    color_palette <- colorRampPalette(c("red", "blue"))
  }
  # title <- paste(title, sep = '_')
  br.out <- barcodeRanks(assays(sce_object)$counts)
  if (logProb_threshold != 0) {
  plot(br.out$rank, br.out$total, log="xy",
       xlab="Rank", ylab="UMI counts", main = title, pch = 19, col = sample_colors)
  }
  else {
  plot(br.out$rank, br.out$total, log="xy",
       xlab="Rank", ylab="UMI counts", main = title, pch = 19,
       col = color_palette(10)[findInterval(color_value, 
                                            seq(min(color_value), 
                                                max(color_value), length.out = 10))])
  }
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=br.out@metadata$knee, col="dodgerblue", lty=2)
  abline(h=br.out@metadata$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
}
par(mfrow=c(1,2))
plot_barcode_rank_plot(sce_object =  sce_1, 
                       title = 'Solo_1_raw',
                       logProb_threshold = TRUE,
                       logProb_thresholds = c(-13000, -5000))
plot_barcode_rank_plot(sce_object =  sce_2, 
                       title = 'Solo_2_raw', 
                       logProb_threshold = TRUE,
                       logProb_thresholds = c(-13000, -5000))
# 