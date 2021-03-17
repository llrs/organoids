# Import the data from the files on servidor2-ciberehd.upc.es to an R object
library("tximport")
files <- list.files(path = "mapped/rsem", pattern = "*.genes.results",
                    full.names = TRUE)
samples <- gsub("\\.genes\\.results", "", gsub(".*/", "", files))
names(files) <- samples
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
saveRDS(txi.rsem, "data/tximport.RDS")
