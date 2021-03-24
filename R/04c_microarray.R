# Download and use a known dataset for testing what might be happening
library("GEOquery")
library("limma")
library("umap")
library("maptools") # point labels without overlaps

# load series and platform data from GEO

gset <- getGEO("GSE123553", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) {
  idx <- grep("GPL20650", attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
    ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
}

# box-and-whisker plot
png(filename = "Figures/microarray_elena.png")
title <- paste("GSE123553", "/", annotation(gset), sep = "")
boxplot(ex, boxwex = 0.7, notch = T, main = title, outline = FALSE, las = 2)
dev.off()

# expression value distribution plot
par(mar = c(4, 4, 2, 1))
title <- paste("GSE123553", "/", annotation(gset), " value distribution", sep = "")
plotDensities(ex, main = title, legend = F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main = "Mean variance trend, GSE123553")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ] # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main = "UMAP plot, nbrs=15", xlab = "", ylab = "", pch = 20, cex = 1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method = "SANN", cex = 0.6)


meta <- pData(gset)

mm <- model.matrix(~`stimulus_time:ch1`+`stimulus_conc:ch1`+ `characteristics_ch1.1`, meta)
mm2 <- cbind(`(Intercept)` = mm[, 1], meta[, c("patient id:ch1", "stimulus_time:ch1", "stimulus_conc:ch1", "characteristics_ch1.1")])

mm3 <- cbind(mm, patient = as.numeric(meta[, c("patient id:ch1")]))
colnames(mm3)[2:8] <- c("24h", "3h", "6h", "5mM", "c", "EpOCs", "patient")

# Paired
fit <- lmFit(ex, design = mm3)
fit3 <- eBayes(fit)
sdt_m <- summary(decideTests(fit3))

# Unpaired
fit2 <- lmFit(ex, design = mm3[, 1:7])
fit4 <- eBayes(fit2)
sdt2 <- summary(decideTests(fit4))

# Difference of ~ 200 genes

# Multilimma doesn't work because is tailred to RNAseq
res_cyclicloess_unpaired <- multilimma(ex, mm3[ , 2:7], nmethod = "cyclicloess")
res_cyclicloess_paired <- multilimma(ex, mm3[ , 2:7], nmethod = "cyclicloess", paired = mm3[, "patient"])
res_quantile_unpaired <- multilimma(ex, mm3[ , 2:7], nmethod = "quantile")
res_quantile_paired <- multilimma(ex, mm3[ , 2:7], nmethod = "quantile", paired = mm3[, "patient"])
