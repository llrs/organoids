# Look at the quality of the data: library size and PCA
# No batch effects but sample 604 is different
# After talking with Aida is the only one that it is a biopsy, the other 3 are cirugy.
library("tidyr")
library("plyr")
library("dplyr")
library("stringr")
library("integration")
library("biomaRt")
library("sva")
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

batch <- ifelse(meta2$SAMPLE == "604", "surgery", "biopsy")
mod <- model.matrix(~cond + SAMPLE,  meta2)
sva_combat <- sva::ComBat(counts4, batch = batch, mod = mod)

pdf("Figures/PCAs_new.pdf")
plotPCA(t(counts4), meta2$cond)
plotPCA(t(counts4), ifelse(!is.na(meta2$estimul), "Estímul", "No estímul"))
plotPCA(t(counts4), meta2$SAMPLE)
plotPCA(t(counts4), meta2$PB)
data.pca <- prcomp(t(counts4), scale. = TRUE)
ggbiplot(data.pca, groups = as.character(meta2$SAMPLE), var.axes = FALSE, ellipse = TRUE, circle = FALSE,
         labels = colnames(counts4)) +
  theme_minimal()

plotPCA(t(sva_combat), meta2$cond) + labs(title = "Without sample effect")
plotPCA(t(sva_combat), ifelse(!is.na(meta2$estimul), "Estímul", "No estímul"))  + labs(title = "Without sample effect")

plotPCA(t(sva_combat), as.character(meta2$SAMPLE)) + labs(title = "Without sample effect")
plotPCA(t(sva_combat), meta2$PB) + labs(title = "Without sample effect")
data.pca <- prcomp(t(sva_combat), scale. = TRUE)
ggbiplot(data.pca, groups = as.character(meta2$SAMPLE), var.axes = FALSE, ellipse = TRUE, circle = FALSE,
         labels = colnames(sva_combat)) +
  theme_minimal() + labs(title = "Without sample effect")
dev.off()

# Adding PCAs 2021/05/03 based on TO_DO Potsti2D_bioinfo
tsc <- t(sva_combat)
meta2$PBS <- grepl("PBS$", meta2$cond)

pdf("Figures/PCAs_PBS_estímul_PBS.pdf")
keep <- !is.na(meta2$estimul) & meta2$PBS
plotPCA(tsc[keep, ], meta2$cond[keep]) + labs(title = "Without sample effect")
meta3 <- meta2[keep, ]
meta3$name <- meta3$`Macrogen SAMPLE NAME`

for (e in unique(meta3$estimul)) {
  keep2 <- meta3$estimul == e | meta3$estimul == "PBS"
  if (e == "PBS") {
    next
  }
  a <- plotPCA(tsc[meta3$name[keep2], ], meta3$cond[keep2]) +
    labs(title = paste("Without sample effect", e),
         col = "Condició")
  print(a)
}
dev.off()

pdf("Figures/PCAs_PB_estímul_PB.pdf")
keep <- !meta2$PBS & (!is.na(meta2$PB) | !is.na(meta2$estimul))
plotPCA(tsc[keep, ], meta2$cond[keep]) + labs(title = "Without sample effect")
meta3 <- meta2[keep, ]
meta3$name <- meta3$`Macrogen SAMPLE NAME`

for (e in unique(meta3$estimul)) {
  keep2 <- meta3$estimul %in% e
  if (e == "PBS") {
    next
  }
  a <- plotPCA(tsc[meta3$name[keep2], ], meta3$cond[keep2]) +
    labs(title = paste("Without sample effect, ", e), col = "Condició")
  print(a)
}
dev.off()
