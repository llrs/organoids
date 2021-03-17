library("tidyr")
library("plyr")
library("dplyr")
library("stringr")
library("integration")
library("biomaRt")
source("R/02_functions_juanjo.R")

tx <- readRDS("data/amayorgas_organoids_counts.RDS")
counts <- tx$counts

# Prepare info from the samples
meta <- readxl::read_xlsx("data/DB_postbi2D.xlsx")
colnames(meta)[3:4] <- c("id", "cond")
meta2 <- fill(meta, c(1, 2)) %>%
  mutate(cond = gsub("PBSther", "PBSth", cond),
         PB = str_extract(cond, "PB.{2,}"),
         estimul = str_extract_all(cond, "INFg|TNFa|IL1b|FLA|(PBS$)")) %>%
  rowwise() %>%
  mutate(estimul = paste(estimul, collapse = "+")) %>%
  ungroup() %>%
  mutate(estimul = ifelse(!nzchar(estimul), NA, estimul))

saveRDS(meta2, "data_out/pheno.RDS")

# Reorder samples to match pheno data
counts <- counts[, match(meta2$`Macrogen SAMPLE NAME`, colnames(counts))]
# Filter genes without expression
var_counts <- apply(counts, 1, var)
counts2 <- counts[var_counts != 0, ]
dim(counts2)
counts3 <- filter_RNAseq(counts2)
dim(counts3)

r <- t(simplify2array(strsplit(rownames(counts3), "\\.[0-9]+_")))
colnames(r) <- c("Ensembl", "name")

# Filter by protein coding
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gb <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
            filters = c("ensembl_gene_id", "biotype"),
            values = list(r[, 1], "protein_coding"),
            mart = ensembl)
counts4 <- counts3[which(r[, 1] %in% gb$ensembl_gene_id), ]
saveRDS(counts4, "data_out/counts.RDS")

# Plots
png("Figures/library_size.png")
barplot(log10(sort(colSums(counts3))))
dev.off()

pdf("Figures/PCAs.pdf")
plotPCA(t(counts4), meta2$cond)
plotPCA(t(counts4), ifelse(!is.na(meta2$estimul), "Estímul", "No estímul"))
plotPCA(t(counts4), meta2$SAMPLE)
plotPCA(t(counts4), meta2$PB)
data.pca <- prcomp(t(counts4), scale. = TRUE)
ggbiplot(data.pca, groups = as.character(meta2$SAMPLE), var.axes = FALSE, ellipse = TRUE, circle = FALSE,
         labels = colnames(counts4)) +
  theme_minimal()
dev.off()
