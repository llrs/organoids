# Explores the pathways from IPA on the biopsyes RNAseq
# Got a reminder on 03/06/2021 to look at the pathways from 2D organoids on
# biopsies
library("limma")
library("readxl")
library("dplyr")
library("stringr")
library("purrr")
library("ggplot2")
library("lubridate")
library("EnsDb.Hsapiens.v86")
library("fgsea")
library("ComplexHeatmap")

# RNAseq ####
db <- readxl::read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = "n.a.")
conn <- gzfile("data/voom.RNAseq.data.all.cal.noduplications.tsv.gz")
rna <- read.table(conn, sep = "\t", check.names = FALSE)
rna <- rna[ , !grepl(" reseq$", colnames(rna))] # Remove as said
db2 <- db[db$Sample_id %in% colnames(rna), ]
# Clean unnecessary rows
keep_cols <- c(colnames(db2)[1:13], colnames(db2)[20:58])
db2 <- unique(db2[, keep_cols])
#Same metadata than samples
stopifnot(nrow(db2) == ncol(rna))
db2 <- db2[match(db2$Sample_id, colnames(rna)), ]

# Show samples by disease, location and ulcers
db2 %>%
  count(Ulcers, sample_location, IBD) %>%
  arrange(sample_location, IBD, Ulcers) %>%
  dplyr::select(sample_location, IBD, Ulcers, n)

# Limit to protein coding genes
genes <- integration::trimVer(rownames(rna))
meta_genes <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = genes,
                       keytype = "GENEID",
                       columns = c("SYMBOL", "GENEBIOTYPE"))
raw_genes <- data.frame(rownames = rownames(rna), genes = genes)
meta_genes2 <- merge(raw_genes, meta_genes, by.x = "genes", by.y = "GENEID", all = TRUE)
rna <- rna[meta_genes2$GENEBIOTYPE %in% "protein_coding", ]

# Split the data by sample location
db_ileum <- db2[db2$sample_location == "ileum" & db2$IBD != "ctrl", ]
db_colon_cd <- db2[db2$sample_location != "ileum" & db2$IBD == "CD", ]
db_colon_uc <- db2[db2$sample_location != "ileum" & db2$IBD == "UC", ]
rna_ileum <- rna[ , db_ileum$Sample_id]
rna_colon_cd <- rna[ , db_colon_cd$Sample_id]
rna_colon_uc <- rna[ , db_colon_uc$Sample_id]

remove_unvariant <- function(x) {
  x[apply(x, 1, var) != 0, ]
}

rna_ileum <- remove_unvariant(rna_ileum)
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
paste0("c", comp, "_")
# cocktail
  c24_sign_FLA+PBS_vs_PBS	c23_sign_IL1b+PBS_vs_PBS	c21_sign_INFg+PBS_vs_PBS	c22_sign_TNFa+PBS_vs_PBS	c49_sign_INFg+TNFa_vs_PBS	c55_sign_INFg+TNFa+PBS_vs_PBS

  )

com <- paste0(compar$`Variable comparativa`, "_vs_", compar$`Referencia comparativa`)
com2 <- paste0("c", compar$`Número comparativa`, "_sign_", com)
which(gsub("fc_", "", colnames(r$fc)) %in% com)
# TODO select genes UP, etc and then  make heatmaps on the next section
g1 <- categ_fc(r$p[, 1], r$fdr[, 1], r$fc[, 1], 0)
g2 <- categ_fc(r$p[, 2], r$fdr[, 2], r$fc[, 2], 0)
g3 <- categ_fc(r$p[, 3], r$fdr[, 3], r$fc[, 3], 0)
g4 <- categ_fc(r$p[, 4], r$fdr[, 4], r$fc[, 4], 0)
g5 <- categ_fc(r$p[, 5], r$fdr[, 5], r$fc[, 5], 0)
g6 <- categ_fc(r$p[, 6], r$fdr[, 6], r$fc[, 6], 0)
g1e <- trimVer(a[, com2[1], drop = TRUE])

# Check that it is the same as from excel files
names(g1) <- NULL
g1[g1 == ""] <- NA
ts <- table(g1e, g1, useNA = "ifany" )
diag(ts) <- 0
stopifnot(all(ts == 0))

# 2021/06/04 Azu: heatmap antiTNF week 0 & week 4/46 responders with regulated genes for each stimuli
# heatmap antiTNF ####
meta_genes_excel_org <- a[, c("genes", "Ensembl...122", "name...123")]
genes_excel_org <- as.matrix(a[, com2])

db2 %>% dplyr::select(IBD, Sample_id, Ulcers, sample_location, week)
db_colon_cd$Week <- ifelse(db_colon_cd$week == 0, "0", "14/46")
db_colon_uc$Week <- ifelse(db_colon_uc$week == 0, "0", "14/46")
db_ileum$Week <- ifelse(db_ileum$week == 0, "0", "14/46")

plot_heat <- function(mat, db, title, file_head) {
  for (i in 1:6) {
    genes_of_interest <- meta_genes_excel_org$Ensembl...122[!is.na(genes_excel_org[, i])]
    mat <- as.matrix(mat)
    mat_median <- apply(mat, 1, function(x){median(x)})
    mat_norm <- sweep(mat, 1, mat_median)
    # heatmap(rna_colon_cd_norm[trimVer(rownames(rna_colon_cd)) %in% genes_of_interest, ])

    ha <- HeatmapAnnotation(df = as.data.frame(db[, c("Ulcers", "Week")]),
                            col = list(Ulcers = c("yes" = "red", "no" = "green"),
                                       Week = c("0" = "brown", "14/46" = "blue")))
    png(paste0("Figures/", file_head, "heatmap_antiTNF_genes_organoids_", com[i], ".png"))
    h <- Heatmap(mat_norm[trimVer(rownames(mat)) %in% genes_of_interest, ],
                 show_row_names = FALSE, show_column_names = FALSE,
                 top_annotation = ha, column_title = paste(com[i], title),
                 name = "norm. expr.")
    draw(h)
    dev.off()
  }
}

plot_heat(rna_ileum, db_ileum, "on ileum CD", "ileum_CD_")
plot_heat(rna_colon_cd, db_colon_cd, "on colon CD", "colon_CD_")
plot_heat(rna_colon_uc, db_colon_uc, "on colon UC", "colon_UC_")
