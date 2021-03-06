# Compare the samples according to a file with Juanjo function
# Not many differences even accounting for paired data
library("dplyr")
source("R/02_functions_juanjo.R")
compar <- readxl::read_xlsx("output/comparatives.xlsx") # Manually made
meta2 <- readRDS("output/pheno.RDS")
counts <- readRDS("output/counts.RDS")
order_samples <- readxl::read_xlsx("output/ordre samples.xlsx") # Azu manually made

# On 21/04/2021 Decided to remove outliers
rm_outliers <- TRUE
outliers <- "605"

# On a mail the 2021/03/16, Azu decided to do this comparisons
# Adding missing comparison 49 on 2021/03/17 and placing it right after the PBS_vs_*
# On 2021/03/22 decided to add comparisons 59:61
# Adding Comparison 55 on 2021/03/25
# reordering comparisons on 2021/05/08
# On 16/04/2021 The PBS should be on the comparisons of the estimulated cells
# Manually added on comparisons.xlsx to be included here 61:74
# On 21/04/2021 Azu told comparison 49 shouldn't be used.
comp <- c(5:24, 49, 25:39, 59:61, 55, 51:54)
comp1 <- c(24, 23, 21, 22, 55, 61:59, # PBS
           19, 18, 16, 20, 17, # PBK
           14, 13, 11, 15, 12, #PBNissel
           9, 8, 6, 10, 7, # PBSth
           62:64, # aTNF
           65:67, # FLA
           68:70, # IL1B
           71:74, #INFg
           75:77,# INFg+TNFa
           54:51 # Controls
           )
stopifnot("43 comparatives" = length(comp1) == 43,
          "Each comparative only once" = length(unique(table(comp1))) == 1)
sub_compar <- compar[comp1, ]

all(sub_compar$`Referencia comparativa` %in% meta2$cond)
all(sub_compar$`Variable comparativa` %in% meta2$cond)

if (rm_outliers) {
  meta2 <- meta2[!meta2$SAMPLE %in% outliers, ]
  counts <- counts[, meta2$`Macrogen SAMPLE NAME`]
}

l <- vector("list", nrow(sub_compar))
for (comp in seq_along(l)) {
  x <- rep(NA, nrow(meta2))
  x[meta2$cond == sub_compar$`Referencia comparativa`[comp]] <- 1
  x[meta2$cond == sub_compar$`Variable comparativa`[comp]] <- 2
  l[[comp]] <- x
}

comp_01 <- simplify2array(l)
colnames(comp_01) <- paste0(sub_compar$`Variable comparativa`, "_vs_", sub_compar$`Referencia comparativa`)
rownames(comp_01) <- meta2$`Macrogen SAMPLE NAME`
if (rm_outliers)  {
  saveRDS(comp_01, "output/samples_comparisons_wo_outliers.RDS")
} else {
  saveRDS(comp_01, "output/samples_comparisons.RDS")
}

res <- multilimma(counts, comp_01, "cyclicloess", meta2$SAMPLE)
if (rm_outliers)  {
  saveRDS(res, "output/limma_juanjo_wo_outliers.RDS")
} else {
  saveRDS(res, "output/limma_juanjo.RDS")
}

tr <- make_TestResults(res) # Transform to limma TestResults object to see a summary

# GETS ####
# * Create matrix for significance ####
diff <- matrix(dimnames = dimnames(res$fc), ncol = ncol(res$fc), nrow = nrow(res$fc))
diff[res$p < 0.05 & res$fc > 0] <- "UP"
diff[res$p < 0.05 & res$fc < 0] <- "DW"
diff[res$fdr < 0.05 & res$fc > 0] <- "UUP"
diff[res$fdr < 0.05 & res$fc < 0] <- "DDW"

colnames(diff) <- gsub("fc_", "sign_", colnames(diff))

# * Number the comparisons ####
add_comparative <- function(comp, mat) {
  colnames(mat) <- paste0("c", comp, "_", colnames(mat))
  mat
}

diff <- add_comparative(sub_compar$`Número comparativa`, diff)
p <- add_comparative(sub_compar$`Número comparativa`, res$p)
fdr <- add_comparative(sub_compar$`Número comparativa`, res$fdr)
fc <- add_comparative(sub_compar$`Número comparativa`, res$fc)

# * Prepare gene information ####
r <- t(simplify2array(strsplit(rownames(diff), "\\.[0-9]+_")))
colnames(r) <- c("Ensembl", "name")
rownames(r) <- rownames(diff)

# Append gene information
out <- cbind(r, diff, fc, fdr, p, r)

m <- meta2[, c("cond", "Macrogen SAMPLE NAME", "SAMPLE", "PB", "estimul")] %>%
  arrange(cond, SAMPLE)

# Requests on order first control then PBS then the rest
order_samples2 <- cbind(order_samples, SAMPLE = c(602, 604, 605, 607))
m2 <- merge(order_samples2, m, all = TRUE, sort = FALSE)

if (rm_outliers) {
  m2 <- m2[!is.na(m2$`Macrogen SAMPLE NAME`), ]
}

stopifnot(all(table(m2$`Macrogen SAMPLE NAME`) == 1))
counts <- counts[, match(m2$`Macrogen SAMPLE NAME`, colnames(counts))]
stopifnot(all(table(colnames(counts)) == 1))

# Use normalized data to flatten those outliers and make more reduce the range
cn <- integration::norm_RNAseq(counts)
# Scale can be done with GETS or in R.
# Scale by median
# cn <- cn - apply(cn, 1, median)
# cm <- cbind(genes = rownames(counts), cn)
# Scale by scaling and centering
cn <- apply(cn, 1, scale)
rownames(cn) <- colnames(counts)
cm <- cbind(genes = rownames(counts), t(cn))

# * Export data in the right format for GETS ####
if (rm_outliers)  {
  write.table(cm, file = "output/GETS_matrix_wo_outliers.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
  system2("gzip", args = "-kf output/GETS_matrix_wo_outliers.tsv") # Compress it to upload to website
  write.table(out, file = "output/GETS_gene_wo_outliers.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
  system2("gzip", args = "-kf output/GETS_gene_wo_outliers.tsv") # Compress it to upload to website
  write.table(m2,
              file = "output/GETS_sample_wo_outliers.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
  df <- data.frame(type = c("SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO"),
                   column = c("PB", "PB", "PB", "SAMPLE"),
                   value = c("PBK12", "PBNissle", "PBSth", "604"),
                   color = c("MAGENTA", "BLUE", "GREEN", "RED"))
  write.table(df, file = "output/GETS_colors_wo_outliers.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
} else {
  write.table(cm, file = "output/GETS_matrix.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
  system2("gzip", args = "-kf output/GETS_matrix.tsv") # Compress it to upload to website
  write.table(out, file = "output/GETS_gene.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
  system2("gzip", args = "-kf output/GETS_gene.tsv") # Compress it to upload to website
  write.table(m2,
              file = "output/GETS_sample.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
  df <- data.frame(type = c("SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO"),
                   column = c("PB", "PB", "PB", "SAMPLE"),
                   value = c("PBK12", "PBNissle", "PBSth", "604"),
                   color = c("MAGENTA", "BLUE", "GREEN", "RED"))
  write.table(df, file = "output/GETS_colors.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}
# Process it on http://bioinfo.ciberehd.org/GETS/index.html
# Options: ....sample information file? & ... Rows and Centered?  True.
