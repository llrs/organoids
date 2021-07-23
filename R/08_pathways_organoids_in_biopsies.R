# Explores the pathways from IPA on the biopsyes RNAseq
# Got a reminder on 03/06/2021 to look at the pathways from 2D organoids on
# biopsies
library("limma")
library("readxl")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("EnsDb.Hsapiens.v86")
library("dplyr")
library("tidyr")
library("fgsea")
library("ComplexHeatmap")
library("fgsea")
library("GSVA")
trimVer <- integration::trimVer

# RNAseq ####
db <- readxl::read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = "n.a.")
# Using this data as there shouldn't be no problem with w0
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)
rna <- rna[ , !grepl(" reseq$", colnames(rna))] # Remove as said
db2 <- db[db$Sample_id %in% colnames(rna), ]
# Clean unnecessary columns
keep_cols <- c(colnames(db2)[1:13], colnames(db2)[20:58])
db2 <- unique(db2[, keep_cols])

#Same metadata than samples ####
stopifnot(nrow(db2) == ncol(rna))
db2 <- db2[match(db2$Sample_id, colnames(rna)), ]

db2 <- mutate(db2, Week = case_when(week == 0 & IBD != "ctrl" ~ "0",
                                    week != 0 & IBD != "ctrl" ~ "14/46",
                                    IBD == "ctrl" ~ "ctrl",
                                    TRUE ~ NA_character_))

# Show samples by disease, location and ulcers
db2 %>%
  count(Ulcers, sample_location, week, IBD) %>%
  arrange(sample_location, IBD, week, Ulcers) %>%
  dplyr::select(sample_location, IBD, week, Ulcers, n)

# Limit to protein coding genes
genes <- trimVer(rownames(rna))
meta_genes <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = genes,
                       keytype = "GENEID",
                       columns = c("SYMBOL", "GENEBIOTYPE"))
raw_genes <- data.frame(rownames = rownames(rna), genes = genes)
meta_genes2 <- merge(raw_genes, meta_genes, by.x = "genes", by.y = "GENEID", all = TRUE)
rna <- rna[meta_genes2$GENEBIOTYPE %in% "protein_coding", ]

# Split the data by sample location
db_ileum <- db2[db2$sample_location == "ileum" & db2$week %in% c("0", "ctrl"), ]
db_colon <- db2[db2$sample_location == "colon" & db2$week %in% c("0", "ctrl"), ]
db_colon_cd <- db_colon[db_colon$IBD == "CD", ]
db_colon_uc <- db_colon[db_colon$IBD == "UC", ]
rna_ileum <- rna[ , db_ileum$Sample_id]
rna_colon <- rna[ , db_colon$Sample_id]
rna_colon_cd <- rna[ , db_colon_cd$Sample_id]
rna_colon_uc <- rna[ , db_colon_uc$Sample_id]

remove_unvariant <- function(x) {
  x[apply(x, 1, var) != 0, ]
}

rna_ileum <- remove_unvariant(rna_ileum)
rna_colon <- remove_unvariant(rna_colon)
rna_colon_cd <- remove_unvariant(rna_colon_cd)
rna_colon_uc <- remove_unvariant(rna_colon_uc)

# Comparisons ####

categ_fc <- function(p.value, adj.p.val, FC, threshold = 1.5) {
  binary <- p.value
  binary[p.value > 0.05 | abs(FC) < threshold] <- ""
  binary[p.value <= 0.05 & FC >= threshold] <- "UP"
  binary[p.value <= 0.05 & FC <= -threshold] <- "DW"
  binary[adj.p.val <= 0.05 & FC >= threshold] <- "UUP"
  binary[adj.p.val <= 0.05 & FC <= -threshold] <- "DDW"
  binary
}

clean_limma <- function(rna, db, meta_genes, plot = FALSE) {
  browse()
  model <- model.matrix(~ Week + Ulcers, db)
  fit <- eBayes(lmFit(rna, model))
  if (plot) {
    volcanoplot(fit)
  }
  tt <- topTable(fit, coef = 2, number = Inf, resort.by = "logFC")
  tt <- merge(tt, meta_genes, by.x = "row.names", by.y = "rownames", sort = FALSE)
  FC <- gtools::logratio2foldchange(tt$logFC)
  tt$FC <- FC
  tt$sign <- categ_fc(tt$P.Value, tt$adj.P.val, FC)
  tt[c("SYMBOL", "sign", "FC", "logFC", "adj.P.Val", "P.Value")]
}

# tt_ileum <- clean_limma(rna_ileum, db_ileum, meta_genes2)
# tt_colon_cd <- clean_limma(rna_colon_cd, db_colon_cd, meta_genes2)
# tt_colon_uc <- clean_limma(rna_colon_uc, db_colon_uc, meta_genes2)


# DEG organoids ####
# From comparatives_13.xls
r <- readRDS("output/limma_juanjo.RDS")
s <- readRDS("output/samples_comparisons.RDS")
a <- read_xlsx("output/comparatives_v13.xlsx", skip = 5)
compar <- readxl::read_xlsx("output/comparatives.xlsx") # Manually made
comp <- c(24, 23, 21, 22, 49, 55)
compar <- compar[comp, ]

com <- paste0(compar$`Variable comparativa`, "_vs_", compar$`Referencia comparativa`)
com2 <- paste0("c", compar$`Número comparativa`, "_sign_", com)

g1 <- categ_fc(r$p[, 1], r$fdr[, 1], r$fc[, 1], 0)
g2 <- categ_fc(r$p[, 2], r$fdr[, 2], r$fc[, 2], 0)
g3 <- categ_fc(r$p[, 3], r$fdr[, 3], r$fc[, 3], 0)
g4 <- categ_fc(r$p[, 4], r$fdr[, 4], r$fc[, 4], 0)
g5 <- categ_fc(r$p[, 5], r$fdr[, 5], r$fc[, 5], 0)
g6 <- categ_fc(r$p[, 6], r$fdr[, 6], r$fc[, 6], 0)

sign_comp <- cbind(g1, g2, g3, g4, g5, g6)
colnames(sign_comp) <- com
k <- apply(sign_comp, 1, function(x){any(x != "")})
sign_comp_org <- sign_comp[k, ]

# 2021/06/04 Azu: heatmap antiTNF week 0 & week 14/46 responders with regulated genes for each stimuli
# heatmap antiTNF ####
meta_genes_excel_org <- a[, c("genes", "Ensembl...122", "name...123")]
colnames(meta_genes_excel_org) <- c("genes", "Ensembl", "name")
rownames(sign_comp_org) <- meta_genes_excel_org$Ensembl[match(rownames(sign_comp_org), meta_genes_excel_org$genes)]
genes_excel_org <- as.matrix(a[, com2])

plot_heat <- function(mat, db, title, file_head) {
  for (i in 1:6) {
    genes_of_interest <- meta_genes_excel_org$Ensembl...122[!is.na(genes_excel_org[, i])]
    mat <- as.matrix(mat)
    mat_median <- apply(mat, 1, function(x){median(x)})
    mat_norm <- sweep(mat, 1, mat_median)
    # heatmap(rna_colon_cd_norm[trimVer(rownames(rna_colon_cd)) %in% genes_of_interest, ])

    ha <- HeatmapAnnotation(df = as.data.frame(db[, c("Ulcers", "week", "Responder")]),
                            col = list(Ulcers = c("yes" = "red", "no" = "green"),
                                       Responder = c("R" = "pink", "NR" = "gray"),
                                       week = c("0" = "brown", "14" = "blue", "46" = "darkblue")))
    png(paste0("Figures/", file_head, "heatmap_antiTNF_genes_organoids_", com[i], ".png"))
    h <- Heatmap(mat_norm[trimVer(rownames(mat)) %in% genes_of_interest, ],
                 show_row_names = FALSE, show_column_names = FALSE,
                 top_annotation = ha, column_title = paste(com[i], title),
                 name = "norm. expr.")
    draw(h)
    dev.off()
  }
}

plot_heat(rna_ileum, db_ileum, "on ileum CD", "20210609_ileum_CD_")
plot_heat(rna_colon_cd, db_colon_cd, "on colon CD", "20210609_colon_CD_")
plot_heat(rna_colon_uc, db_colon_uc, "on colon UC", "20210609_colon_UC_")


# 2021/06/08 Summarize the content of the plots, but only from colon.
# Sort samples according to week and ulcers.
meta_genes_excel_org$genes_ver <- gsub("(ENSG[0-9]*.[0-9]+).+", "\\1", meta_genes_excel_org$genes)

diff_organoids_complete <- list(g1, g2, g3, g4, g5, g6)
names(diff_organoids_complete) <- com
diff_organoids <- lapply(diff_organoids_complete, function(x, meta){
  diff <- names(x)[x != ""]
  meta[, 2, drop = TRUE][meta$genes %in% diff]
}, meta = meta_genes_excel_org)
g_organoids <- unique(unlist(diff_organoids, FALSE, FALSE))


g <- strcapture("(^(ENSG[0-9]*).[0-9]+)", rownames(rna_colon),
           proto = data.frame(ensembl_ver = character(), ensembl = character()))

g_colon <- g[g$ensembl %in% g_organoids, ]
rna_colon_org <- rna_colon[g_colon$ensembl_ver, ]
rownames(rna_colon_org) <- g_colon$ensembl

col_ann <- list(Ulcers = c("yes" = "red", "no" = "green"),
                Responder = c("R" = "pink", "NR" = "black"),
                IBD = c("UC" = "darkgreen", "CD" = "blue", "ctrl" = "lightblue"),
                week = c("0" = "brown", "14" = "blue", "46" = "darkblue"),
                Week = c("0" = "brown", "14/46" = "blue"))

# * GSVA ####
gsva_colon <- gsva(as.matrix(rna_colon_org), diff_organoids)

db_colon <- arrange(db_colon, Responder, Week, Ulcers, IBD, Gender)

top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[, c("Ulcers", "IBD", "Responder")]),
                             col = col_ann)
png("Figures/heatmap_GSVA_antiTNF_colon_all_genes_organoids.png")
Heatmap(gsva_colon[ , db_colon$Sample_id],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = top_ann, column_title = "GSVA colon",
        name = "GSVA")
dev.off()

png("Figures/heatmap_GSVA_antiTNF_colon_uc_genes_organoids.png")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[db_colon$IBD == "UC",
                                                         c("Ulcers", "week", "Responder")]),
                             col = col_ann)
Heatmap(gsva_colon[ , db_colon$Sample_id[db_colon$IBD == "UC"]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = top_ann, column_title = "Colon UC; GSVA",
        name = "GSV")
dev.off()

png("Figures/heatmap_GSVA_antiTNF_colon_cd_genes_organoids.png")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[db_colon$IBD == "CD",
                                                         c("Ulcers", "Responder")]),
                             col = col_ann)
Heatmap(gsva_colon[ , db_colon$Sample_id[db_colon$IBD == "CD"]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE,
        top_annotation = top_ann, column_title = "Colon CD; GSVA",
        name = "GSV")
dev.off()


# * All genes ####
# On 2021/06/09 Summarize but not using GSVA plot every gene that is present
db_colon <- arrange(db_colon, week, Ulcers, IBD, Responder, Gender)

png("Figures/heatmap_antiTNF_colon_uc_genes_organoids.png")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[db_colon$IBD == "UC", c("Ulcers")]),
                             col = col_ann)
Heatmap(rna_colon_org[ , db_colon$Sample_id[db_colon$IBD == "UC"]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        top_annotation = top_ann, column_title = "Colon UC w0",
        name = "norm. expr.")
dev.off()

png("Figures/heatmap_antiTNF_colon_cd_genes_organoids.png")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[db_colon$IBD == "CD", c("Ulcers")]),
                             col = col_ann)
Heatmap(rna_colon_org[ , db_colon$Sample_id[db_colon$IBD == "CD"]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        top_annotation = top_ann, column_title = "Colon CD w0",
        name = "norm. expr.")
dev.off()

# On 2021/06/10 reminder to make a heatmap of w0 and controls with the organoids genes
png("Figures/heatmap_antiTNF_colon_cd_w0_ctrls_genes_organoids.png")
db_colon <- arrange(db_colon, IBD, Ulcers)
keep_samples_cd <- db_colon$IBD %in% c("CD", "ctrl")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[keep_samples_cd, c("Ulcers", "IBD")]),
                             col = col_ann)
Heatmap(rna_colon_org[ , db_colon$Sample_id[keep_samples_cd]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        top_annotation = top_ann, column_title = "Colon CD w0 or controls",
        name = "norm. expr.")
dev.off()

norm_colon_cd_ctrl <- rna_colon[, match(db_colon$Sample_id[keep_samples_cd], colnames(rna_colon))]
design <- model.matrix(~Week, data = db_colon[keep_samples_cd, ])
fit1 <- lmFit(norm_colon_cd_ctrl, design)
eb1 <- eBayes(fit1)
plotSA(fit1, main = "Final model: Mean-variance trend")
tt <- topTable(eb1, number = Inf, coef = 2)
tt_g <- merge(tt, g, by.x = "row.names", by.y = "ensembl_ver")
tt_g_organoids_cd <- tt_g[tt_g$ensembl %in% g_organoids, ]
diff_cd <- categ_fc(tt_g_organoids_cd$P.Value, adj.p.val = tt_g_organoids_cd$adj.P.Val,
         FC = gtools::logratio2foldchange(-tt_g_organoids_cd$logFC))
names(diff_cd) <- tt_g_organoids_cd$ensembl


diff_organoids_sign <- lapply(diff_organoids_complete, function(x, meta){
  diff <- x[x != ""]
  names(diff) <- meta$Ensembl[match(names(diff), meta$genes)]
  diff
}, meta = meta_genes_excel_org)

s <- unique(unlist(sapply(diff_organoids_sign, names)))




png("Figures/heatmap_antiTNF_colon_UC_w0_ctrls_genes_organoids.png")
db_colon <- arrange(db_colon, IBD, Ulcers)
keep_samples_uc <- db_colon$IBD %in% c("UC", "ctrl")
top_ann <- HeatmapAnnotation(df = as.data.frame(db_colon[keep_samples_uc, c("Ulcers", "IBD")]),
                             col = col_ann)
Heatmap(rna_colon_org[ , db_colon$Sample_id[keep_samples_uc]],
        cluster_columns = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
        show_column_names = FALSE, show_row_names = FALSE,
        top_annotation = top_ann, column_title = "Colon UC w0 or controls",
        name = "norm. expr.")
dev.off()

norm_colon_uc_ctrl <- rna_colon[, match(db_colon$Sample_id[keep_samples_uc], colnames(rna_colon))]
design <- model.matrix(~Week, data = db_colon[keep_samples_uc, ])
contrasts <- makeContrasts(Weekctrl-Intercept, levels = design)
fit1 <- lmFit(norm_colon_uc_ctrl, design)
eb1 <- eBayes(fit1)
plotSA(fit1, main = "Final model: Mean-variance trend")
tt2 <- topTable(eb1, number = Inf, coef = 2)
tt_g <- merge(tt2, g, by.x = "row.names", by.y = "ensembl_ver")
tt_g_organoids_uc <- tt_g[tt_g$ensembl %in% g_organoids, ]
diff_uc <- categ_fc(tt_g_organoids_uc$P.Value, adj.p.val = tt_g_organoids_uc$adj.P.Val,
         FC = gtools::logratio2foldchange(-tt_g_organoids_uc$logFC))
names(diff_uc) <- tt_g_organoids_uc$ensembl

m_genes <- meta_genes_excel_org[match(rownames(sign_comp_org), meta_genes_excel_org$Ensembl), ]
diff_antiTNF <- cbind(diff_uc, diff_cd)
diff_antiTNF <- diff_antiTNF[match(rownames(sign_comp_org), rownames(diff_antiTNF)), ]

comp_all <- cbind(m_genes[, 2:3], diff_antiTNF, sign_comp_org)
saveRDS(comp_all, "output/w0_vs_ctrl_biopsies_genes_orgnoids.RDS")
writexl::write_xlsx(comp_all, path = "output/w0_vs_ctrl_biopsies_genes_orgnoids.xlsx")


# Genes shared bx organoids ####
s <- sapply(com, function(x, comp_all) {
  comp_co <- as.character(comp_all[, x])
  diff_co <- nzchar(comp_co) & !is.na(comp_co)
  sub_m <- comp_all[diff_co, ]

  comp_uc <- as.character(sub_m[, "diff_uc"])
  comp_cd <- as.character(sub_m[, "diff_cd"])

  diff_uc <- nzchar(comp_uc) & !is.na(comp_uc)
  diff_cd <- nzchar(comp_cd) & !is.na(comp_cd)

  up <- sum(startsWith(comp_uc, "U") & startsWith(as.character(sub_m[, x]), "U"), na.rm = TRUE)
  dw <- sum(startsWith(comp_uc, "D") & startsWith(as.character(sub_m[, x]), "D"), na.rm = TRUE)
  diff_co_uc <- sum(up, dw)
  up <- sum(startsWith(comp_cd, "U") & startsWith(as.character(sub_m[, x]), "U"), na.rm = TRUE)
  dw <- sum(startsWith(comp_cd, "D") & startsWith(as.character(sub_m[, x]), "D"), na.rm = TRUE)
  diff_co_cd <- sum(up, dw)
  c(genes = sum(diff_co),
    UC = sum(diff_uc),
    CD = sum(diff_cd),
    UC_shared_perc =  sum(diff_uc)/sum(diff_co)*100,
    CD_shared_perc =  sum(diff_cd)/sum(diff_co)*100,
    UC_similar_perc = diff_co_uc/sum(diff_uc)*100,
    CD_similar_perc = diff_co_cd/sum(diff_cd)*100
  )
}, comp_all = comp_all)
s2 <- t(s)
s2 <- as.data.frame(s2)
s2$comparative <- rownames(s2)
writexl::write_xlsx(s2, "output/resum_gens_org_biopsies.xlsx")
