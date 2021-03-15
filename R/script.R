# Script from the server to prepare a file for STAR+RSEM bash script

path <- "/home/lrevilla/data/organoids_amayorgas/HN00144024/RawFASTQ/"
lf <- list.files(path)
library("tidyverse")
df <- data.frame(files = lf) %>%
  separate(col = files, into = c("condition", "sample", "pair"), sep = "_",
           remove = FALSE) %>%
  mutate(sample = as.numeric(sample),
         condition = as.numeric(condition),
         pair = as.numeric(gsub("\\..+", "", pair))) %>%
  arrange(sample, condition, pair)

df2 <- df %>%
  pivot_wider(id_cols = c(sample, condition), values_from = files,
              names_prefix = "file", names_from = pair) %>%
  mutate(file1 = paste0(path, file1),
         file2 = paste0(path, file2),
         sample = paste(condition, sample, sep = "_"))

df2 %>%
  select(file1, file2, sample) %>%
  write_tsv(file = "paired_samples.tsv", col_names = FALSE)

