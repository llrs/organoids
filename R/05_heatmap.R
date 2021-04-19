# Make heatmaps of some comparisons
library("dplyr")
library("ComplexHeatmap")
res <- readRDS("data_out/limma_juanjo.RDS")
counts <- readRDS("data_out/counts.RDS")
comp <- readRDS("data_out/samples_comparisons.RDS")
meta <- readRDS("data_out/pheno.RDS")

# * Create matrix for significance ####
diff <- matrix(dimnames = dimnames(res$fc), ncol = ncol(res$fc), nrow = nrow(res$fc))
diff[res$p < 0.05 & res$fc > 0] <- "UP"
diff[res$p < 0.05 & res$fc < 0] <- "DW"
diff[res$fdr < 0.05 & res$fc > 0] <- "UUP"
diff[res$fdr < 0.05 & res$fc < 0] <- "DDW"
colnames(diff) <- gsub("fc_", "", colnames(diff))

# Normalize
cn <- integration::norm_RNAseq(counts)

comp_name <- "INFg+PBS_vs_PBS"

sel_genes <- rownames(diff)[!is.na(diff[, comp_name])]
sel_genes_sorted <- names(sort(res$fc[sel_genes, paste0("fc_", comp_name)], decreasing = TRUE))
samples_diff <- colnames(counts)[!is.na(comp[, comp_name])]
df <- meta[meta$`Macrogen SAMPLE NAME` %in% samples_diff, c("cond", "SAMPLE"), drop = FALSE]
df$SAMPLE <- as.factor(df$SAMPLE)
ha <- HeatmapAnnotation(df = as.data.frame(df), show_annotation_name = FALSE,
                        name = "Condition", which = "column",
                        col = list("SAMPLE" = c("602" = "black",
                                   "604" = "yellow",
                                   "605" = "green",
                                   "607" = "orange")))

sel_n <- function(x, n = 50, two.sided = TRUE) {
  l <- length(x)
  if (two.sided){
    n <- ceiling(n/2)
    return(c(x[1:n], x[(l-n):l]))
  }
  x[1:n]
}

genes_heatmap <- sel_n(sel_genes_sorted, 50)
cn_scaled <- t(apply(cn[genes_heatmap, samples_diff], 1, scale))
Heatmap(cn_scaled,
        top_annotation = ha,
        row_labels = gsub("^ENSG[0-9]*\\.[0-9]+_", "", genes_heatmap),
        column_title = paste("Genes differential expressed on", comp_name),
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Z"))

