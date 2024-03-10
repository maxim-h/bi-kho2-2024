library(ggplot2)
library(dplyr)

count_barcodes <- read.csv("../Data/Count_barcodes.csv", sep=';', col.names = c('Alignment', 'Filtering', 'Sample', 'Count'), header = FALSE)

count_barcodes$Groups <- paste(count_barcodes$Alignment, '_', count_barcodes$Filtering, '_', count_barcodes$Sample)

ggplot(count_barcodes, aes(y = Count, x = Groups, fill = Groups))+
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        legend.position = "none")
