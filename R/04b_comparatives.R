# Make comparisons as I would do them (defining contrast and design matrices)
library("dplyr")
library("limma")
compar <- readxl::read_xlsx("data_out/comparatives.xlsx") # Manually made
meta <- readRDS("data_out/pheno.RDS")
counts <- readRDS("data_out/counts.RDS")
order_samples <- readxl::read_xlsx("data_out/ordre samples.xlsx") # Azu manually made

# On a mail the 2021/03/16, Azu decided to do this comparisons
# Adding missing comparison 49 on 2021/03/17 and placing it right after the PBS_vs_*
# On 2021/03/22 decided to add comparisons 59:61
# Adding Comparison 55 on 2021/03/25
comp <- c(5:24, 49, 25:39, 59:61, 55, 51:54)
sub_compar <- compar[comp, ]

stopifnot(all(sub_compar$`Referencia comparativa` %in% meta$cond))
stopifnot(all(sub_compar$`Variable comparativa` %in% meta$cond))

expr_edge <- edgeR::DGEList(counts)
expr_edge <- edgeR::calcNormFactors(expr_edge)
cpm <- edgeR::cpm(expr_edge, normalized.lib.sizes = TRUE, log = TRUE)
v <- voom(expr_edge, normalize.method = "cyclicloess", plot = TRUE)
cpm <- v$E # v$E <-normalised matrix


mm <- model.matrix(~cond, meta)


# Apply function to get contrasts
cm <- apply(sub_compar, 1, function(x, m) {
  y <- rep(0, length(m))
  y[x[["Referencia comparativa"]] == m] <- 1
  y[x[["Variable comparativa"]] == m] <- -1
  y
}, m = gsub("^cond", "", colnames(mm)))


colnames(cm) <- paste0(sub_compar$`Referencia comparativa`, "_vs_", sub_compar$`Variable comparativa`)
rownames(cm) <- colnames(mm)
# Set contrasts manually ####
# * Unpaired ####
fit <- lmFit(cpm, design = mm)
fit2 <- contrasts.fit(fit, cm)
fit3 <- eBayes(fit2)
sdt <- summary(decideTests(fit3))

# * Paired ####
mm_p <- cbind(mm, paired  = as.numeric(meta$SAMPLE))
fit_p <- lmFit(cpm, design = mm)
fit2_p <- contrasts.fit(fit_p, cm)
fit3_p <- eBayes(fit2_p)
sdt_p <- summary(decideTests(fit3_p))

stopifnot("Summary is not the same" = all(sdt_p == sdt))
stopifnot("Paired and unpaired return different results" = all(decideTests(fit3_p) == decideTests(fit3)))

# Complex systems ####
# * each factor as separate ####
mm_complex <- inteRmodel::model_RGCCA(meta, c("PB", "estimul"), intercept = TRUE)
mm_complex <- cbind( mm_complex, PBSth = as.numeric(meta$PB == "PBSth"),
                    INFg = as.numeric(meta$estimul == "INFg"),
                    paired  = as.numeric(meta$SAMPLE))
colnames(mm_complex)[1] <- "Intercept"
stopifnot("Check removing interactions" = all(grepl("+", colnames(mm_complex)[7:12])))
mm_complex <- mm_complex[, -c(7:12)]
mm_complex[grep("TNFa", meta$cond), "TNFa"] <- 1
mm_complex[grep("INFg", meta$cond), "INFg"] <- 1
mm_complex[grep("PBS\b", meta$cond), "PBS"] <- 1
mm_complex[grep("FLA", meta$cond), "FLA"] <- 1

mm_complex[is.na(mm_complex)] <- 0

fit_c1 <- lmFit(cpm, design = mm_complex)

mc2 <- makeContrasts("INFg_vs_INFg+TNFa" = (INFg + TNFa) - INFg,
                     "PBSth_vs_INFg+PBSth" = (INFg + PBSth) - PBSth,
                     "PBSth_vs_TNFa+PBSth" = (TNFa + PBSth) - PBSth,
                     "PBSth_vs_IL1b+PBSth" = (IL1b + PBSth) - PBSth,
                     "PBSth_vs_FLA+PBSth" = (FLA + PBSth) - PBSth,
                     "PBSth_vs_INFg+TNFa+PBSth" = (INFg + TNFa + PBSth) - PBSth,
                     "PBNissle_vs_INFg+PBNissle" = (INFg + PBNissle) - PBNissle,
                     "PBNissle_vs_TNFa+PBNissle" = (TNFa + PBNissle) - PBNissle,
                     "PBNissle_vs_IL1b+PBNissle" = (IL1b + PBNissle) - PBNissle,
                     "PBNissle_vs_FLA+PBNissle" = (FLA + PBNissle) - PBNissle,
                     "PBNissle_vs_INFg+TNFa+PBNissle" = (INFg + TNFa + PBNissle) - PBNissle,
                     "PBK12_vs_INFg+PBK12" = (INFg + PBK12) - PBK12,
                     "PBK12_vs_TNFa+PBK12" = (TNFa + PBK12) - PBK12,
                     "PBK12_vs_IL1b+PBK12" = (IL1b + PBK12) - PBK12,
                     "PBK12_vs_FLA+PBK12" = (FLA + PBK12) - PBK12,
                     "PBK12_vs_INFg+TNFa+PBK12" = (INFg + TNFa + PBK12) - PBK12,
                     "PBS_vs_INFg+PBS" = (INFg + PBS) - PBS,
                     "PBS_vs_TNFa+PBS" = (TNFa + PBS) - PBS,
                     "PBS_vs_IL1b+PBS" = (IL1b + PBS) - PBS,
                     "PBS_vs_FLA+PBS" = (FLA + PBS) - PBS,
                     "PBS_vs_INFg+TNFa" = (INFg + TNFa) - PBS,
                     "INFg_vs_INFg+PBSth" = (INFg + PBSth) - INFg,
                     "TNFa_vs_TNFa+PBSth" = (TNFa + PBSth) - TNFa,
                     "IL1b_vs_IL1b+PBSth" = (IL1b + PBSth) - IL1b,
                     "FLA_vs_FLA+PBSth" = (FLA + PBSth) - FLA,
                     "INFg+TNFa_vs_INFg+TNFa+PBSth" = (INFg + TNFa + PBSth) - (INFg + TNFa),
                     "INFg_vs_INFg+PBNissle" = (INFg + PBNissle) - INFg,
                     "TNFa_vs_TNFa+PBNissle" = (TNFa + PBNissle) - TNFa,
                     "IL1b_vs_IL1b+PBNissle" = (IL1b + PBNissle) - IL1b,
                     "FLA_vs_FLA+PBNissle" = (FLA + PBNissle) - FLA,
                     "INFg+TNFa_vs_INFg+TNFa+PBNissle" = (INFg + TNFa + PBNissle) - (INFg + TNFa),
                     "INFg_vs_INFg+PBK12" = (PBK12 + INFg) - INFg,
                     "TNFa_vs_TNFa+PBK12" = (PBK12 + TNFa) - TNFa,
                     "IL1b_vs_IL1b+PBK12" = (PBK12 + IL1b) - IL1b,
                     "FLA_vs_FLA+PBK12" = (PBK12 + FLA) - FLA,
                     "INFg+TNFa_vs_INFg+TNFa+PBK12" = (INFg + TNFa + PBK12) - (INFg + TNFa),
                     "PBS_vs_PBSth" = PBSth - PBS,
                     "PBS_vs_PBNissle" = PBNissle - PBS,
                     "PBS_vs_PBK12" = PBK12 - PBS,
                     "Control_vs_PBSth" = PBSth,
                     "Control_vs_PBNissle" = PBNissle,
                     "Control_vs_PBK12" = PBK12,
                     "Control_vs_PBS" = PBS,
                     levels = mm_complex)

stopifnot("Not all comparisons included" = all(colnames(mc2) == colnames(cm)))
fit2_c1 <- contrasts.fit(fit_c1, mc2)
fit3_c1 <- eBayes(fit2_c1)
sdt_c1 <- summary(decideTests(fit3_c1))
# We obtain different results! so

# * Factor pasted ####
mm_complex2 <- model.matrix(~0+cond, meta)
# mm_complex2 <- cbind(mm_complex2, "paired" = meta$SAMPLE)

v <- voom(expr_edge, design = mm_complex2, normalize.method = "cyclicloess", plot = TRUE)
cpm <- v$E # v$E <-normalised matrix

colnames(mm_complex2) <- gsub("^cond", "", colnames(mm_complex2))
colnames(mm_complex2) <- make.names(colnames(mm_complex2))

mc3 <- makeContrasts("INFg_vs_INFg+TNFa" = INFg - INFg.TNFa,
                     "PBSth_vs_INFg+PBSth" = PBSth - INFg.PBSth,
                     "PBSth_vs_TNFa+PBSth" = PBSth - TNFa.PBSth,
                     "PBSth_vs_IL1b+PBSth" = PBSth - IL1b.PBSth,
                     "PBSth_vs_FLA+PBSth" = PBSth - FLA.PBSth,
                     "PBSth_vs_INFg+TNFa+PBSth" = PBSth - INFg.TNFa.PBSth,
                     "PBNissle_vs_INFg+PBNissle" = PBNissle - INFg.PBNissle,
                     "PBNissle_vs_TNFa+PBNissle" = PBNissle - TNFa.PBNissle,
                     "PBNissle_vs_IL1b+PBNissle" = PBNissle - IL1b.PBNissle,
                     "PBNissle_vs_FLA+PBNissle" = PBNissle - FLA.PBNissle,
                     "PBNissle_vs_INFg+TNFa+PBNissle" = PBNissle - INFg.TNFa.PBNissle,
                     "PBK12_vs_INFg+PBK12" = PBK12 - INFg.PBK12,
                     "PBK12_vs_TNFa+PBK12" = PBK12 - TNFa.PBK12,
                     "PBK12_vs_IL1b+PBK12" = PBK12 - IL1b.PBK12,
                     "PBK12_vs_FLA+PBK12" = PBK12 - FLA.PBK12,
                     "PBK12_vs_INFg+TNFa+PBK12" = PBK12 - INFg.TNFa.PBK12,
                     "PBS_vs_INFg+PBS" = PBS - INFg.PBS,
                     "PBS_vs_TNFa+PBS" = PBS - TNFa.PBS,
                     "PBS_vs_IL1b+PBS" = PBS - IL1b.PBS,
                     "PBS_vs_FLA+PBS" = PBS - FLA.PBS,
                     "PBS_vs_INFg+TNFa" = PBS - INFg.TNFa,
                     "INFg_vs_INFg+PBSth" = INFg - INFg.PBSth,
                     "TNFa_vs_TNFa+PBSth" = TNFa - TNFa.PBSth,
                     "IL1b_vs_IL1b+PBSth" = IL1b - IL1b.PBSth,
                     "FLA_vs_FLA+PBSth" = FLA - FLA.PBSth,
                     "INFg+TNFa_vs_INFg+TNFa+PBSth" = INFg.TNFa - INFg.TNFa.PBSth,
                     "INFg_vs_INFg+PBNissle" = INFg - INFg.PBNissle,
                     "TNFa_vs_TNFa+PBNissle" = TNFa - TNFa.PBNissle,
                     "IL1b_vs_IL1b+PBNissle" = IL1b - IL1b.PBNissle,
                     "FLA_vs_FLA+PBNissle" = FLA - FLA.PBNissle,
                     "INFg+TNFa_vs_INFg+TNFa+PBNissle" = INFg.TNFa - INFg.TNFa.PBNissle,
                     "INFg_vs_INFg+PBK12" = INFg - INFg.PBK12,
                     "TNFa_vs_TNFa+PBK12" = TNFa - TNFa.PBK12,
                     "IL1b_vs_IL1b+PBK12" = IL1b - IL1b.PBK12,
                     "FLA_vs_FLA+PBK12" = FLA - FLA.PBK12,
                     "INFg+TNFa_vs_INFg+TNFa+PBK12" =  INFg.TNFa - INFg.TNFa.PBK12,
                     "PBS_vs_PBSth" = PBS - PBSth,
                     "PBS_vs_PBNissle" = PBS - PBNissle,
                     "PBS_vs_PBK12" = PBS - PBK12,
                     "Control_vs_PBSth" = Control - PBSth,
                     "Control_vs_PBNissle" = Control - PBNissle,
                     "Control_vs_PBK12" = Control - PBK12,
                     "Control_vs_PBS" = Control - PBS,
                     levels = mm_complex2)
stopifnot("Not all comparisons included" = all(colnames(mc3) == colnames(cm)))
corfit <- duplicateCorrelation(cpm, design = mm_complex2, block = meta$SAMPLE)
fit_c2 <- lmFit(cpm, design = mm_complex2, block = meta$SAMPLE, correlation = corfit$consensus.correlation)
fit_c2b <- lmFit(cpm, design = mm_complex2)

# Section 9.7 of limma user guide
fit2_c2 <- contrasts.fit(fit_c2, mc3)
fit3_c2 <- eBayes(fit2_c2)
sdt_c2 <- summary(decideTests(fit3_c2))

fit2_c2b <- contrasts.fit(fit_c2b, mc3)
fit3_c2b <- eBayes(fit2_c2b)
sdt_c2b <- summary(decideTests(fit3_c2b))

# Compare multiple methods ####
# Check samples picked  up
c0 <- comp_01 # Samples on each comparison manually made for multilimma
# Changing to compare with the method limma really uses
c0[is.na(c0)] <- 0
c0[c0 == 2] <- -1
colnames(c0) <- colnames(mc2)
c1 <- samples_contrasts(mm_complex, mc2)
c2 <- samples_contrasts(mm_complex2, mc3)
if (!all(c0 == c2)) {
  c02_col <- unique(which(c0 != c2, arr.ind = TRUE)[, "col"])
  colnames(mc3)[c02_col] # complex 2 and original are the same
} else {
  message("Selecting same samples as manually")
}
c01_col <- unique(which(c0 != c1, arr.ind = TRUE)[, "col"])
colnames(mc3)[c01_col] # Check these
# Model complex 1 is different from model complex 2!

# Compare results visually
resum_multilimma <- t(summary(tr)) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq)
View(resum_multilimma, "juanjo")
write.csv(resum_multilimma, "data_out/resum_multilimma_aparellat.csv",
          row.names = FALSE)
resum_limma_aparellat_9.7 <- t(sdt_c2) %>%
  as.data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq)
View(resum_limma_aparellat_9.7, "complex2")
write.csv(resum_limma_aparellat_9.7, "data_out/resum_limma_aparellat_9.7.csv",
          row.names = FALSE)
l <- lapply(colnames(mc3), topTable, fit = fit3_c2, number = Inf)
logFC <- aveExpr <- ts <- p <- fdr <- matrix(nrow = nrow(cpm), ncol = length(l),
                                             dimnames = list(rownames(cpm), colnames(mc3)))
for (i in seq_along(l)) {
  logFC[, i] <- l[[i]][, "logFC"]
  aveExpr[, i] <- l[[i]][, "AveExpr"]
  ts[, i] <- l[[i]][, "t"]
  p[, i] <- l[[i]][, "P.Value"]
  fdr[, i] <- l[[i]][, "adj.P.Val"]
}
colnames(logFC) <- paste0("c", sub_compar$`Número comparativa`, "_logFC_", colnames(logFC))
colnames(aveExpr) <- paste0("c", sub_compar$`Número comparativa`, "_AveExpr_", colnames(aveExpr))
colnames(ts) <- paste0("c", sub_compar$`Número comparativa`, "_t_", colnames(ts))
colnames(p) <- paste0("c", sub_compar$`Número comparativa`, "_p_", colnames(p))
colnames(fdr) <- paste0("c", sub_compar$`Número comparativa`, "_fdr_", colnames(fdr))
fc <- gtools::logratio2foldchange(logFC)
colnames(fc) <- gsub("_logFC_", "_fc_", colnames(fc))


# * Create matrix for significance ####
diff <- matrix(dimnames = dimnames(fc), ncol = ncol(fc), nrow = nrow(fc))
diff[p < 0.05 & fc > 1.5] <- "UP"
diff[p < 0.05 & fc < 1.5] <- "DW"
diff[fdr < 0.05 & fc > 1.5] <- "UUP"
diff[fdr < 0.05 & fc < 1.5] <- "DDW"
colnames(diff) <- gsub("fc_", "sign_", colnames(diff))


# * Prepare gene information ####
r <- t(simplify2array(strsplit(rownames(diff), "\\.[0-9]+_")))
colnames(r) <- c("Ensembl", "name")
rownames(r) <- rownames(diff)

# Append gene information
out <- cbind(r, diff, fc, fdr, p, r)
