# Compare the samples according to a file with Juanjo function
# Not many differences even accounting for paired data
library("dplyr")
source("R/02_functions_juanjo.R")
compar <- readxl::read_xlsx("data_out/comparatives.xlsx") # Manually made
meta2 <- readRDS("data_out/pheno.RDS")
counts <- readRDS("data_out/counts.RDS")
order_samples <- readxl::read_xlsx("data_out/ordre samples.xlsx") # Azu manually made

# On a mail the 2021/03/16, Azu decided to do this comparisons
# Adding missing comparison 49 on 2021/03/17 and placing it right after the PBS_vs_*
# On 2021/03/22 decided to add comparisons 59:61
# Adding Comparison 55 on 2021/03/25
# reordering comparisons on 2021/05/08
comp <- c(5:24, 49, 25:39, 59:61, 55, 51:54)
comp1 <- c(24, 23, 21, 22, 49, 55, 61:59, # PBS
           19, 18, 16, 20, 17, # PBK
           14, 13, 11, 15, 12, #PBNissel
           9, 8, 6, 10, 7, # PBSth
           36, 31, 26, # aTNF
           38, 33, 28, # FLA
           37, 32, 27, # IL1B
           35, 30, 25, 5, #INFg
           39, 34, 29,# INFg+TNFa
           54:51 # Controls
           )
stopifnot("44 comparatives" = length(comp1) == 44,
          "Each comparative only once" = length(unique(table(comp1))) == 1)
sub_compar <- compar[comp1, ]

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
tr <- make_TestResults(res) # Transform to limma TestResults object to see a summary

mRP <- multiRankProd(counts, comp_01, "cyclicloess", meta2$SAMPLE)
saveRDS(mRP, "data_out/rankprod_paired.RDS")

genes <- strsplit(rownames(mRP[[1]][[1]]), "\\.[0-9]+_")
genes <- t(simplify2array(genes))
colnames(genes) <- c("ENSEMBL", "name")
res_RP <- lapply(mRP, function(x) {
  y <- x$pval
  colnames(y) <- paste0("p.val_", colnames(y))
  y2 <- apply(y, 2, p.adjust, method = "fdr")
  colnames(y2) <- paste0("fdr_", colnames(y2))
  o <- cbind(x$AveFC, y2, y)
  colnames(o)[1] <- "fc_"
  o
  })
cond <- paste0(sub_compar$`Referencia comparativa`, "_vs_",
               sub_compar$`Variable comparativa`)
for (i in seq_along(res_RP)) {
  colnames(res_RP[[i]]) <- paste0("c",sub_compar$`Número comparativa`[i], "_", colnames(res_RP[[i]]), cond[i])
}
sign_fc <- fc <- matrix(ncol = length(res_RP), nrow = nrow(res_RP[[1]]))
p.val <- p.adj <- matrix(ncol = length(res_RP)*2, nrow = nrow(res_RP[[1]]))

colnames(fc) <- seq_len(ncol(fc))
colnames(sign_fc) <- seq_len(ncol(sign_fc))
colnames(p.val) <- seq_len(ncol(p.val))
colnames(p.adj) <- seq_len(ncol(p.adj))
for (i in seq_along(res_RP)) {
  fc[, i] <- res_RP[[i]][, 1]
  colnames(fc)[i] <- colnames(res_RP[[i]])[1]
  sign_fc[, i] <- ifelse(res_RP[[i]][, 4] < 0.05 & res_RP[[i]][, 5] > 0.05, "UP",
                         ifelse(res_RP[[i]][, 5] < 0.05 & res_RP[[i]][, 4] > 0.05, "DW",
                                ifelse(res_RP[[i]][, 2] < 0.05, "UUP",
                                       ifelse(res_RP[[i]][, 3] < 0.05, "DW", NA))))
  colnames(sign_fc)[i] <- gsub("_fc_", "_sign_", colnames(res_RP[[i]])[1])
  # The strange i*2 is to account for 2 columns
  cols <- c(i*2 - 1, i*2)
  p.adj[, cols] <- res_RP[[i]][, 2:3]
  colnames(p.adj)[cols] <- colnames(res_RP[[i]])[2:3]
  p.val[, cols] <- res_RP[[i]][, 4:5]
  colnames(p.val)[cols] <- colnames(res_RP[[i]])[4:5]
}

o <- cbind.data.frame(genes, sign_fc, fc, genes)
write.table(o, file = "data_out/GETS_gene_RankProd.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf data_out/GETS_gene_RankProd.tsv") # Compress it to upload to website

# summary(tr)
# res_quantile <- multilimma(counts, comp_01, "quantile", meta2$SAMPLE)
# saveRDS(res_quantile, "data_out/limma_paired_quantile.RDS")
#
# res_unparied_cyclicloess <- multilimma(counts, comp_01, "cyclicloess")
# saveRDS(res_unparied_cyclicloess, "data_out/limma_unpaired_cyclicloess.RDS")
# res_unparied_quantile <- multilimma(counts, comp_01, "quantile")
# saveRDS(res_unparied_quantile, "data_out/limma_unpaired_quantile.RDS")

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
counts <- counts[, match(m2$`Macrogen SAMPLE NAME`, colnames(counts))]

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
write.table(cm, file = "data_out/GETS_matrix.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
system2("gzip", args = "-kf data_out/GETS_matrix.tsv") # Compress it to upload to website
write.table(out, file = "data_out/GETS_gene.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = "")
system2("gzip", args = "-kf data_out/GETS_gene.tsv") # Compress it to upload to website
write.table(m2,
            file = "data_out/GETS_sample.tsv", quote = FALSE, sep = "\t", row.names = FALSE, na = " ")
df <- data.frame(type = c("SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO", "SAMPLEINFO"),
           column = c("PB", "PB", "PB", "SAMPLE"),
           value = c("PBK12", "PBNissle", "PBSth", "604"),
           color = c("MAGENTA", "BLUE", "GREEN", "RED"))
write.table(df, file = "data_out/GETS_colors.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
# Process it on http://bioinfo.ciberehd.org/GETS/index.html
# Options: ....sample information file? & ... Rows and Centered?  True.
