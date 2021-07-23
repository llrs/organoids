# Export data so Aida can make her own plots
# genes : "CHI3L1", "CEACAM6", "NOD2"
library("tidyr")
library("dplyr")

counts <- readRDS("output/counts.RDS")
meta <- readRDS("output/pheno.RDS")

# Select conditions
meta2 <- meta %>%
  mutate(PBS = endsWith(trimws(estimul), "PBS")) %>%
  filter((!is.na(PB) & !PBS) | (is.na(PB) & PBS) | (PBS & is.na(PB) & is.na(PBS)))

# Normalize

cn <- integration::norm_RNAseq(counts[, meta2$`Macrogen SAMPLE NAME`])

names_genes <- c("CHI3L1", "CEACAM6", "NOD2", "CARD9", "SAA2-SAA4")

# Genes matrix
genes <- rownames(cn)
g <- strsplit(genes, "\\.[0-9]*_")
g <- t(simplify2array(g))
colnames(g) <- c("ENSEMBL", "name")
g <- cbind.data.frame(genes = genes, g)
sel <- g[g$name %in% names_genes, ]

cn_genes <- cn[sel$genes, ]


as.data.frame(cn_genes) %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  left_join(meta2, by = c("name" = "Macrogen SAMPLE NAME")) %>%
  left_join(as.data.frame(sel), by = c("rowname" = "genes")) %>%
  select(name.y, value, cond, PB, estimul, name.x, SAMPLE) %>%
  rename(Gene = name.y, sample = name.x) %>%
  arrange(Gene, cond, PB, estimul) %>%
  writexl::write_xlsx("output/genes_CHI3L1_CEACAM6_NOD2_CARD9_SAA_21_conditions.xlsx")
