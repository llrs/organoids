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
clean_limma <- function(rna, db, meta_genes, plot = FALSE) {
  model <- model.matrix(~ Ulcers, db)
  fit <- eBayes(lmFit(rna, model))
  if (plot) {
    volcanoplot(fit)
  }
  tt <- topTable(fit, coef = 2, number = Inf, resort.by = "logFC")
  tt <- merge(tt, meta_genes, by.x = "row.names", by.y = "rownames", sort = FALSE)
  FC <- gtools::logratio2foldchange(tt$logFC)
  tt$FC <- FC
  binary <- tt$P.Value
  binary[tt$P.Value > 0.05 | abs(FC) < 1.5] <- ""
  binary[tt$P.Value <= 0.05 & FC >= 1.5] <- "UP"
  binary[tt$P.Value <= 0.05 & FC <= -1.5] <- "DW"
  binary[tt$adj.P.Val <= 0.05 & FC >= 1.5] <- "UUP"
  binary[tt$adj.P.Val <= 0.05 & FC <= -1.5] <- "DDW"
  tt$sign <- binary
  tt[c("SYMBOL", "sign", "FC", "logFC", "adj.P.Val", "P.Value")]
}

tt_ileum <- clean_limma(rna_ileum, db_ileum, meta_genes2)
tt_colon_cd <- clean_limma(rna_colon_cd, db_colon_cd, meta_genes2)
tt_colon_uc <- clean_limma(rna_colon_uc, db_colon_uc, meta_genes2)


# raw IPA data ####
excels <- list.files("data/IPA/", full.names = TRUE)

raw <- lapply(excels, read_xls, skip = 1, .name_repair = "universal")
names(raw) <- sub(".xls", "", basename(excels))
upstream_genes <- lapply(raw[c("cock", "fla", "il1b", "inf", "tnf")], function(x) {
  unique(toupper(x$Upstream.Regulator))
})

ipa_sel <- function(tt, genes) {
  u_genes <- unique(toupper(unlist(upstream_genes, FALSE, FALSE)))
  tt[tt$SYMBOL %in% u_genes, ]
}

tt_fgsea <- function(tt, genes, plot = FALSE, ...) {
  x <- tt$FC
  names(x) <- tt$SYMBOL
  r <- fgsea(genes, x, ...)
  if (plot) {
   plotGseaTable(genes, x, r)
  }
  r
}

tt_fgsea(tt_ileum, upstream_genes, nperm = 1000)
tt_fgsea(tt_colon_cd, upstream_genes, nperm = 1000)
tt_fgsea(tt_colon_uc, upstream_genes, nperm = 1000)
