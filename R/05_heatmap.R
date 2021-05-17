# Make heatmaps of some comparisons
library("dplyr")
library("ComplexHeatmap")
library("circlize")
res <- readRDS("output/limma_juanjo.RDS")
counts <- readRDS("output/counts.RDS")
comp <- readRDS("output/samples_comparisons.RDS")
meta <- readRDS("output/pheno.RDS")

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
comp_name2 <- "TNFa+PBS_vs_PBS"
comp_name3 <- "INFg+TNFa+PBS_vs_PBS"

sel_genes <- function(x, comp) {
  comps <- grepl(comp, x = colnames(x), fixed = TRUE)
  if (sum(comps) > 1) {
    comps <- colnames(x) %in% comp
  }
  if (sum(comps) > 1) {
    warning("More than 1 comparison selected")
  }
  rownames(x)[!is.na(x[, comps])]
}
c1 <- sel_genes(diff, comp_name)
c3 <- sel_genes(diff, comp_name3)
genes_heatmap <- intersect(c1, c3)

coi <- c(comp_name, comp_name2, comp_name3)
sel <- apply(comp[, coi], 1, function(x){any(!is.na(x))})
samples_diff <- names(sel)[sel]
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

cn_scaled <- t(apply(cn[genes_heatmap, samples_diff], 1, scale))
colnames(cn_scaled) <- colnames(cn[, samples_diff])
hms <- Heatmap(cn_scaled,
        top_annotation = ha,
        row_labels = gsub("^ENSG[0-9]*\\.[0-9]+_", "", genes_heatmap),
        column_title = paste("Genes differential expressed on", comp_name, "and", comp_name3),
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Z"))

order_heatmap <- colnames(cn_scaled)[column_order(hms)]
meta[match(order_heatmap, meta$`Macrogen SAMPLE NAME`), ] %>% View("dendo")

s <- sapply(coi, grepl, x = colnames(res$fc), fixed = TRUE)
s <- apply(s, 1, any)
df2 <- res$fc[genes_heatmap, s]
col_fun <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(df2,
        row_labels = gsub("^ENSG[0-9]*\\.[0-9]+_", "", genes_heatmap),
        column_title = paste("Genes differential expressed on", comp_name, "and", comp_name3),
        column_labels = gsub("fc_", "", colnames(df2)),
        column_names_rot = 0,
        column_names_centered = TRUE,
        row_names_gp = gpar(fontsize = 8),
        col = col_fun,
        heatmap_legend_param = list(title = "fc"))

# 2021/05/03 adding a heatmap with diff between three postbiotics vs PBS
# Sther vs PBS, Nissle vs PBS and K12 vs PBS
comp_name <- "PBSth_vs_PBS"
comp_name2 <- "PBNissle_vs_PBS"
comp_name3 <- "PBK12_vs_PBS"

c1 <- sel_genes(diff, comp_name)
c2 <- sel_genes(diff, comp_name2)
c3 <- sel_genes(diff, comp_name3)

g <- unique(c(c1, c2, c3)) # Genes DEG in any
g1 <- unique(c(setdiff(c1, c2), setdiff(c2, c3), setdiff(c1, c3))) # Genes DEG in two
g3 <- intersect(intersect(c1, c2), c3) # Genes DEG in the three


make_heatmaps <- function(g, name) {
  comps <- c(comp_name, comp_name2, comp_name3)
  s <- colnames(res$fc) %in% paste0("fc_", comps)
  df2 <- res$fc[g, s]
  col_fun <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))

  h_fc <- Heatmap(df2,
                  show_row_names = FALSE,
                  column_title = paste("Genes differential expressed on postbiotics vs PBS"),
                  column_labels = gsub("fc_", "", colnames(df2)),
                  column_names_rot = 0,
                  column_names_centered = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  col = col_fun,
                  heatmap_legend_param = list(title = "fc"))

  png(paste0("Figures/", name,"_comparatives.png"))
  print(h_fc)
  dev.off()

  sel <- apply(comp[, comps], 1, function(x){any(!is.na(x))})
  samples_diff <- names(sel)[sel]
  cn_scaled <- t(apply(cn[g, samples_diff], 1, scale))
  colnames(cn_scaled) <- colnames(cn[, samples_diff])
  df <- meta[meta$`Macrogen SAMPLE NAME` %in% samples_diff, c("cond", "SAMPLE", "Macrogen SAMPLE NAME"), drop = FALSE]
  df <- df %>%
    mutate(cond = factor(cond, levels = c("PBS", "PBK12", "PBNissle", "PBSth"))) %>%
    rename(sample = `Macrogen SAMPLE NAME`) %>%
    arrange(cond, SAMPLE, sample)

  cn_scaled <- cn_scaled[, match(df$sample, colnames(cn_scaled))]
  ha <- HeatmapAnnotation(df = as.data.frame(df[, -3]), show_annotation_name = FALSE,
                          name = "Condition", which = "column",
                          col = list("SAMPLE" = c("602" = "black",
                                                  "604" = "yellow",
                                                  "605" = "green",
                                                  "607" = "orange")))
  hms <- Heatmap(cn_scaled,
                 top_annotation = ha,
                 show_row_names = FALSE,
                 column_title = paste("Genes differential expressed on postbiotics vs PBS"),
                 row_names_gp = gpar(fontsize = 8),
                 show_column_names = FALSE,
                 cluster_columns = FALSE,
                 heatmap_legend_param = list(title = "Z"))
  png(paste0("Figures/", name, "_samples.png"))
  print(hms)
  dev.off()

  # Filtering for those genes not DEG on PBS vs controls
  comp_name4 <- "PBS_vs_Control"
  c4 <- sel_genes(diff, comp_name4)
  g2 <- setdiff(g, c4)

  df2 <- res$fc[g2, s]
  h_fc <- Heatmap(df2,
                  show_row_names = FALSE,
                  column_title = paste("Genes differential expressed on postbiotics vs PBS"),
                  column_labels = gsub("fc_", "", colnames(df2)),
                  column_names_rot = 0,
                  column_names_centered = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  col = col_fun,
                  heatmap_legend_param = list(title = "fc"))

  png(paste0("Figures/", name, "_comparatives_filtered.png"))
  print(h_fc)
  dev.off()

  sel <- apply(comp[, comps], 1, function(x){any(!is.na(x))})
  samples_diff <- names(sel)[sel]
  cn_scaled <- t(apply(cn[g2, samples_diff], 1, scale))
  colnames(cn_scaled) <- colnames(cn[, samples_diff])
  df <- meta[meta$`Macrogen SAMPLE NAME` %in% samples_diff, c("cond", "SAMPLE", "Macrogen SAMPLE NAME"), drop = FALSE]
  df <- df %>%
    mutate(cond = factor(cond, levels = c("PBS", "PBK12", "PBNissle", "PBSth"))) %>%
    rename(sample = `Macrogen SAMPLE NAME`) %>%
    arrange(cond, SAMPLE, sample)

  cn_scaled <- cn_scaled[, match(df$sample, colnames(cn_scaled))]

  ha <- HeatmapAnnotation(df = as.data.frame(df[, -3]), show_annotation_name = FALSE,
                          name = "Condition", which = "column",
                          col = list("SAMPLE" = c("602" = "black",
                                                  "604" = "yellow",
                                                  "605" = "green",
                                                  "607" = "orange")))
  hms <- Heatmap(cn_scaled,
                 top_annotation = ha,
                 show_row_names = FALSE,
                 column_title = paste("Genes differential expressed on postbiotics vs PBS"),
                 row_names_gp = gpar(fontsize = 8),
                 show_column_names = FALSE,
                 cluster_columns = FALSE,
                 heatmap_legend_param = list(title = "Z"))
  png(paste0("Figures/", name, "_samples_filtered.png"))
  print(hms)
  dev.off()
}

make_heatmaps(g, "heatmap_genes_diff_any")
make_heatmaps(g1, "heatmap_genes_diff_2")
make_heatmaps(g3, "heatmap_genes_diff_3")
