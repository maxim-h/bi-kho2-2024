source('Basic_metrics.R')
library(scuttle)

counts <- c('../Data/Solo.out_1/Gene/raw/matrix.mtx.gz',
            '../Data/Solo.out_2/Gene/raw/matrix.mtx.gz')

genes <- c('../Data/Solo.out_1/Gene/raw/features.tsv.gz',
           '../Data/Solo.out_2/Gene/raw/features.tsv.gz')

barcodes <- c('../Data/Solo.out_1/Gene/raw/barcodes.tsv.gz',
              '../Data/Solo.out_2/Gene/raw/barcodes.tsv.gz')

filenames <- c('Solo1', 'Solo2')


count_tables <- Matrix::readMM(counts[1])
rownames(count_tables) <- readr::read_tsv(genes[1], col_names = FALSE)$X1
colnames(count_tables) <- readr::read_tsv(barcodes[1], col_names = FALSE)$X1

counts_data_empty <- count_tables[,Matrix::colSums(count_tables) < 100 & Matrix::colSums(count_tables) > 0]
counts_data_intermediate <- count_tables[,Matrix::colSums(count_tables) < 1000 & Matrix::colSums(count_tables) > 100]





