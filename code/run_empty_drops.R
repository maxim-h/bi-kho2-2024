library(DropletUtils)
library(glue)

#Считывание данных

#Here I only work with GeneFull data, Gene data is not considered
base_dir <- './data/Species_mixing'
samples_dirs <- c(glue('{base_dir}/Solo.out_1/GeneFull/raw'),
             glue('{base_dir}/Solo.out_2/GeneFull/raw'))
names(samples_dirs) <- c('1_raw', '2_raw')
sce_1 <- read10xCounts(samples_dirs[1], type = 'sparse')
sce_2 <- read10xCounts(samples_dirs[2], type = 'sparse')


# Добавляю точку перегиба и inflection с knee-plot-а  в метадату SCE объекта - почему бы и нет
br.out <- barcodeRanks(assays(sce_1)$counts)
sce_1@metadata$knee <-  br.out@metadata$knee
sce_1@metadata$inflection <-  br.out@metadata$inflection
br.out <- barcodeRanks(assays(sce_2)$counts)
sce_2@metadata$knee <-  br.out@metadata$knee
sce_2@metadata$inflection <-  br.out@metadata$inflection
#П.С цикл фор не хочет работать фор сам ризон
# for (sce_object in c(sce_1, sce_2)) {
#   br.out <- barcodeRanks(assays(sce_object)$counts)
#   sce_object@metadata$knee <-  br.out@metadata$knee
#   sce_object@metadata$inflection <-  br.out@metadata$inflection
# }


# Запуск empty_Drops - считает не сказать чтоб супер быстро, пару минут мб. Думаю разумно запускать с lower 10-30. Дефолт - 100
sce_1_emptyDrops_10_df <- emptyDrops(sce_1, retain = sce_1@metadata$knee, 
                                  test.ambient = TRUE, lower = 10) # lower - порог по UMI counts, ниже которого точно пустые капли
                                                                   # retain - порог по UMI counts, выше которого точно хорошие капли
                                                                   # test.ambient - проверять даже те, что ниже lower
sce_2_emptyDrops_10_df <- emptyDrops(sce_2, retain = sce_1@metadata$knee, 
                                   test.ambient = TRUE, lower = 10)


add_fdrs_to_barcodes <- function(sce_object, empty_drops_df) { 
  empty_drops_df$FDR[is.na(empty_drops_df$FDR)] <- 1 #выставляю всем NA FDR единичку
  sce_object@colData$EmptyDrops_fdr <- empty_drops_df$FDR
  return(sce_object)
} #добавляет инфу по FDRs в SCE object, дальше можно фильтровать
sce_1 <- add_fdrs_to_barcodes(sce_1, sce_1_emptyDrops_10_df)
sce_2 <- add_fdrs_to_barcodes(sce_2, sce_2_emptyDrops_10_df)


plot_barcode_rank_plot <- function(sce_object, title='title', 
                                   color1='red', 
                                   color2='blue', 
                                   color_by = 'emptyDrops',
                                   fdr_threshold = 0.05) {
  if (color_by == 'origin') { #можно покрасить по происхождению (если сразу два образца лежит в SCE)
    sample_colors <- ifelse(as.factor(colData(sce_object)$Sample) == levels(as.factor(colData(sce_object)$Sample))[1], 
                            adjustcolor(color1, alpha.f = 1), 
                            adjustcolor(color2, alpha.f = .2))
  }
  if (color_by == 'emptyDrops') { #можно красить по результатам empty_drops
    sample_colors <- ifelse((sce_object$EmptyDrops_fdr < fdr_threshold), 
                            adjustcolor(color1, alpha.f = 1), 
                            adjustcolor(color2, alpha.f = .2))
  }
  title <- paste(title,fdr_threshold, sep = '_')
  br.out <- barcodeRanks(assays(sce_object)$counts)
  plot(br.out$rank, br.out$total, log="xy",
       xlab="Rank", ylab="UMI counts", main = title, pch = 19, col = sample_colors)
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=br.out@metadata$knee, col="dodgerblue", lty=2)
  abline(h=br.out@metadata$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
} #строит knee-plot, может красить капли по условию

#Можно делать Knee-плоты
par(mfrow=c(2,2))
for (fdr_threshold in c(0.05, 0.01, 0.001, 0.0001)) { 
plot_barcode_rank_plot(sce_object =  sce_1, 
                       title = 'Solo_1_raw',
                       color_by = 'emptyDrops',
                       fdr_threshold = fdr_threshold)
}

par(mfrow=c(2,2))
for (fdr_threshold in c(0.05, 0.01, 0.001, 0.0001)) { 
  plot_barcode_rank_plot(sce_object =  sce_2, 
                         title = 'Solo_2_raw',
                         color_by = 'emptyDrops',
                         fdr_threshold = fdr_threshold)
}
