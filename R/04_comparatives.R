# Compare the samples according to a file
library("dplyr")
source("R/02_functions_juanjo.R")
compar <- readxl::read_xlsx("data_out/comparatives.xlsx") # Manually made
meta2 <- readRDS("data_out/pheno.RDS")
counts <- readRDS("data_out/counts.RDS")
order_samples <- readxl::read_xlsx("data_out/ordre samples.xlsx") # Azu manually made

# On a mail the 2021/03/16, Azu decided to do this comparisons
# Adding missing comparison 49 on 2021/03/17 and placing it right after the PBS_vs_*
sub_compar <- compar[c(5:24, 49, 25:39), ]

all(sub_compar$`Referencia comparativa` %in% meta2$cond)
all(sub_compar$`Variable comparativa` %in% meta2$cond)

l <- vector("list", nrow(sub_compar))
for (comp in seq_along(l)) {
  x <- rep(NA, nrow(meta2))
  x[meta2$cond == sub_compar$`Referencia comparativa`[comp]] <- 1
  x[meta2$cond == sub_compar$`Variable comparativa`[comp]] <- 2
  l[[comp]] <- x
}

comp_01 <- simplify2array(l)
colnames(comp_01) <- paste0(sub_compar$`Referencia comparativa`, "_vs_", sub_compar$`Variable comparativa`)
rownames(comp_01) <- meta2$`Macrogen SAMPLE NAME`

res <- multilimma(counts, comp_01, "cyclicloess", meta2$SAMPLE)

saveRDS(res, "data_out/limma_juanjo.RDS")


diff <- matrix(dimnames = dimnames(res$fc), ncol = ncol(res$fc), nrow = nrow(res$fc))
diff[res$p < 0.05 & res$fc > 0] <- "UP"
diff[res$p < 0.05 & res$fc < 0] <- "DW"
diff[res$fdr < 0.05 & res$fc > 0] <- "UUP"
diff[res$fdr < 0.05 & res$fc < 0] <- "DDW"


colnames(diff) <- gsub("fc_", "sign_", colnames(diff))
add_comparative <- function(comp, mat) {
  colnames(mat) <- paste0("c", comp, "_", colnames(mat))
  mat
}

diff <- add_comparative(sub_compar$`Número comparativa`, diff)
p <- add_comparative(sub_compar$`Número comparativa`, res$p)
fdr <- add_comparative(sub_compar$`Número comparativa`, res$fdr)
fc <- add_comparative(sub_compar$`Número comparativa`, res$fc)
r <- t(simplify2array(strsplit(rownames(diff), "\\.[0-9]+_")))
colnames(r) <- c("Ensembl", "name")
rownames(r) <- rownames(diff)

out <- cbind(r, diff, fc, fdr, p, r)
m <- meta2[, c("cond", "Macrogen SAMPLE NAME", "SAMPLE", "PB", "estimul")] %>%
  arrange(cond, SAMPLE)

# Requests on order First control then PBS then the rest
new_order <- match(order_samples$cond, m$cond)
# TODO continue here as this is wrong as the same sample will be 4x.
m <- m[new_order, ]
counts <- counts[, match(m$`Macrogen SAMPLE NAME`, colnames(counts))]

# Use normalized data to flatten those outliers and make more reduce the range
cm <- cbind(genes = rownames(counts), integration::norm_RNAseq(counts))
write.table(cm, file = "data_out/GETS_matrix.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf data_out/GETS_matrix.tsv") # Compress it to upload to website
write.table(out, file = "data_out/GETS_gene.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf data_out/GETS_gene.tsv") # Compress it to upload to website
write.table(m,
            file = "data_out/GETS_sample.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
df <- data.frame(type = c("SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO"),
           column = c("PB", "PB", "PB", "SAMPLE"),
           value = c("PBK12", "PBNissle", "PBSth", "604"),
           color = c("MAGENTA", "BLUE", "GREEN", "RED"))
write.table(df, file = "data_out/GETS_colors.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
# Process it on http://bioinfo.ciberehd.org/GETS/index.html
# Options: ....sample information file? & ... Rows and Centered?  True.
