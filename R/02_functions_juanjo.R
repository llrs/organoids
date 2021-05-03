# Functions adapted from Juanjo Lozano  to make plots etc
library("ggbiplot")
library("edgeR")

l2fc <- function(logratio, base = 2) {
  retval <- base^(logratio)
  retval <- ifelse(retval < 1, -1 / retval, retval)
  retval
}

assigncol <- function(pheno) {
  colx <- c("blue", "red", "green", "orange", "yellow", "gray", "pink", "white", "brown")
  colxm <- matrix("gray", nrow = length(pheno))
  g <- unique(pheno)

  for (i in 1:length(g)) {
    ss <- which(pheno == g[i])
    colxm[ss] <- colx[i]
  }
  if (length(colx) > 9) {
    stop("N>9: Don't Run")
  }

  return(colxm)
}



plotPCA <- function(xmat, gg) {

  data.class <- gg
  data.pca <- prcomp(xmat, scale. = TRUE)
  colx <- assigncol(data.class)
  # g <- ggbiplot(data.pca,
  #   obs.scale = 0.25, var.scale = 0.5, labels.size = 3,
  #   groups = data.class, var.axes = FALSE, ellipse = TRUE, circle = FALSE
  # )
  #
  # g <- g + theme(
  #   axis.text = element_text(size = 12),
  #   axis.title = element_text(size = 13)
  # )
  #
  #
  # print(g)

  g <- ggbiplot(data.pca,
    obs.scale = 0.25, var.scale = 0.5, labels.size = 3,
    groups = data.class, var.axes = FALSE, labels = rownames(xmat),
    ellipse = TRUE, circle = FALSE
  )

  g <- g + theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )


  g
}

# Classes must be a matrix with each column being a comparative with:
# 1 as reference and 2 as the other samples
# paired can be a vector with the name/factor of the paired samples
multilimma <- function(xdat, classes, nmethod, paired = NULL) {
  stopifnot(ncol(xdat) == nrow(classes))
  # nmethod=
  # cyclicloess
  # none
  # quantile
  fc <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  p <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  fdr <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  t <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])

  colnames(fdr) <- paste("fdr_", colnames(classes), sep = "")
  colnames(fc) <- paste("fc_", colnames(classes), sep = "")
  colnames(p) <- paste("p_", colnames(classes), sep = "")
  colnames(t) <- paste("t_", colnames(classes), sep = "")

  rownames(fc) <- rownames(fdr) <- rownames(p) <- rownames(t) <- rownames(xdat)
  for (i in seq_len(ncol(classes))) {
    nona <- which(!is.na(classes[, i]))
    myclassx <- na.omit(classes[, i])
    myxdat <- xdat[, nona]
    ### norm step
    y <- DGEList(myxdat)
    y <- calcNormFactors(y)
    pdf(paste("data_out/QC_", colnames(classes)[i], "_", nmethod, ".pdf", sep = ""))
    on.exit({if (length(dev.list()) >= 2) dev.off()})
    v <- voom(y, normalize.method = nmethod, plot = TRUE) # v$E <-normalised matrix

    message(i, " ", colnames(classes)[i])
    print(table(classes[, i]))

    norm <- v$E

    BaseMean <- apply(norm, 1, mean)
    Log2R <- apply(norm[, which(myclassx == 2)], 1, mean) - apply(norm[, which(myclassx == 1)], 1, mean)

    boxplot(norm, las = 3)
    plot(BaseMean, Log2R, main = colnames(classes)[i], cex = 0.3)
    abline(h = 0, lwd = 2, col = "blue")
    # to speed only most 1000 genes
    vars <- apply(norm, 1, var)
    sel1000 <- order(vars, decreasing = TRUE)[1:1000]
    plotPCA(t(norm[sel1000, ]), factor(myclassx))
    fc[, i] <- l2fc(Log2R)
    if (!is.null(paired)) {
      stopifnot(length(paired[nona]) == length(myclassx))
      paired_contrast <- paired[nona]
      dd <- data.frame(class = myclassx, paired = paired_contrast)
      mod <- model.matrix(~ class + paired, dd)
      fit1 <- lmFit(norm, mod)

    } else {
      mod <- model.matrix(~ factor(myclassx, levels = c("1", "2")))
      fit1 <- lmFit(norm, mod)
    }
    eb1 <- eBayes(fit1)
    plotSA(fit1, main = "Final model: Mean-variance trend")
    p[, i] <- eb1$p.value[, 2]
    fdr[, i] <- p.adjust(p[, i], method = "BH")
    tt <- topTable(eb1, number = Inf, coef = 2)
    plot(tt$AveExpr, tt$logFC)
    dev.off()
    t[, i] <- eb1$t[, 2]

    lup <- length(which(p[, i] < 0.05 & fc[, i] > 1.5))
    ldw <- length(which(p[, i] < 0.05 & fc[, i] < -1.5))

    message("UP=", lup, ": DW=", ldw, "\n")
  }
  return(list(
    fc = fc, p = p, fdr = fdr, t = t, BaseMean = as.matrix(BaseMean),
    Amean = as.matrix(eb1$Amean)
  ))
}


myplclust <- function(hclust, lab = hclust$labels,
                      lab.col = rep(1, length(hclust$labels)),
                      hang = 0.1, ...) {
  ## modifiction of plclust for plotting hclust objects *in colour*!
  ## Copyright Eva KF Chan 2009
  ## Arguments:
  ##    hclust:    hclust object
  ##    lab:        a character vector of labels of the leaves of the tree
  ##    lab.col:    colour for the labels; NA=default device foreground colour
  ##    hang:     as in hclust & plclust
  ## Side effect:
  ##    A display of hierarchical cluster with coloured leaf labels.


  y <- rep(hclust$height, 2)
  x <- as.numeric(hclust$merge)
  y <- y[which(x < 0)]
  x <- x[which(x < 0)]
  x <- abs(x)
  y <- y[order(x)]
  x <- x[order(x)]
  plot(hclust, labels = FALSE, hang = hang, ...)
  text(x = x, y = y[hclust$order] - (max(hclust$height) * hang),
       labels = lab[hclust$order], col = lab.col[hclust$order], srt = 90,
       adj = c(1, 0.5), xpd = NA, ...)
}


colorhclust <- function(dat, var) {
  dd <- dist(t(dat))
  hh <- hclust(dd, method = "average")
  myplclust(hh, lab = colnames(dat), lab.col = as.numeric(as.factor(var)),
            main = "", lwd = 2, cex = 1, ylab = "", xlab = "",
            ann = FALSE, las = 1)
}


extractPC123 <- function(mydata) {
  # extract from ggbiplot
  pcobj <- prcomp(mydata, scale. = TRUE)
  choices <- 1:3
  obs.scale <- 1
  var.scale <- 1
  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  d <- pcobj$sdev
  u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = "*")
  v <- pcobj$rotation
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
  df.u <- df.u * nobs.factor
  colnames(df.u) <- c("PC1", "PC2", "PC3")
  rownames(df.u) <- rownames(mydata)
  return(df.u)
}

make_TestResults <- function(x) {
  s_fc <- sign(x$fc)
  s_fc[x$fdr > 0.05] <- 0
  colnames(s_fc) <- gsub("^fc_", "", colnames(s_fc))
  class(s_fc) <-  "TestResults"
  s_fc
}


samples_contrasts <- function(mod, contrasts) {
  mod %*% contrasts
}


# Classes must be a matrix with each column being a comparative with:
# 1 as reference and 2 as the other samples
# paired can be a vector with the name/factor of the paired samples
multiRankProd <- function(xdat, classes, nmethod, paired = NULL) {
  stopifnot(ncol(xdat) == nrow(classes))
  # nmethod=
  # cyclicloess
  # none
  # quantile

  l <- vector("list", ncol(classes))
  for (i in seq_len(ncol(classes))) {
    nona <- which(!is.na(classes[, i]))
    myclassx <- na.omit(classes[, i])
    myxdat <- xdat[, nona]
    ### norm step
    y <- DGEList(myxdat)
    y <- calcNormFactors(y)
    pdf(paste("data_out/QC_rankProd", colnames(classes)[i], "_", nmethod, ".pdf", sep = ""))
    on.exit({if (length(dev.list()) >= 2) dev.off()})
    v <- voom(y, normalize.method = nmethod, plot = TRUE) # v$E <-normalised matrix

    message(i, " ", colnames(classes)[i])
    print(table(classes[, i]))

    norm <- v$E

    # Assume paired samples are on the same order
    norm2 <- norm[, which(myclassx == 2)] - norm[, which(myclassx == 1)]
    RP <- RankProducts(norm2, rep(1, ncol(norm2)), logged = TRUE, plot = TRUE, na.rm = FALSE,
                 rand = 246, gene.names = rownames(norm))
    l[[i]] <- RP
    dev.off()
  }
  names(l) <- names(classes)
  l
}
