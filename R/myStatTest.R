myStatTest <- function(exprs,
                       meta,
                       test,
                       request.type,
                       comparison,
                       group.ctrl.test,
                       group.colors,
                       outdir,
                       diff.space,
                       compute.log2fc = T,
                       samplesAreRows = F,
                       test.group = NULL,
                       ctrl.group = NULL,
                       pair_id = NULL) {

  if (is.character(exprs)) {
    exprs <- read.xlsx(exprs, rowNames = 1, check.names = F)
  }

  if (is.character(meta)) {
    meta <- read.xlsx(meta, rowNames = 1, check.names = F)
  }

  # Ensuring data has rows of samples and columns of features
  if (!samplesAreRows) {
    exprs <- t(exprs)
  }

  meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  if (!is.null(pair_id)) {
    pair_id <- pair_id[rownames(exprs)]
  }

  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}
  }

  createDir(outdir)

  if (test == "t-test") {
    # Use explicit groups if provided, otherwise parse from comparison name
    if (is.null(test.group) || is.null(ctrl.group)) {
      test.group <- strsplit(comparison,"_over_")[[1]][1]
      ctrl.group <- strsplit(comparison,"_over_")[[1]][2]
    }

    if ((length(test.group) == 0)|(length(ctrl.group) == 0)) {
      stop(cat("\n######## ERROR WITH COMPARISON NAME:", comparison,"########",
               "\n######## T-test comparison name format is 'test_over_control' ########\n"))
    }
    if (!is.null(pair_id)) {
      pair_test <- pair_id[meta[, 1] == test.group]
      pair_ctrl <- pair_id[meta[, 1] == ctrl.group]
      if (!setequal(pair_test, pair_ctrl)) {
        stop(paste0("Paired t-test for '", comparison, "': pair IDs do not match between '",
                    test.group, "' and '", ctrl.group, "'. Check the Pairing column."))
      }
      ord_test <- order(pair_test)
      ord_ctrl <- order(pair_ctrl)
      all.pvals <- apply(exprs, 2, function(x) {
        xt <- x[meta[, 1] == test.group][ord_test]
        xc <- x[meta[, 1] == ctrl.group][ord_ctrl]
        ifelse(var(xt - xc, na.rm = TRUE) > 0,
               t.test(xt, xc, paired = TRUE)$p.value,
               1)
      })
    } else {
      all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
                                                       t.test(x~meta[,1])$p.value,
                                                       1)})
    }

    if (diff.space == "log2") {
      # here, data is already transformed to log2 space in the pre processing step
      log2FC <- apply(exprs[meta[,1] == test.group, ,drop=F], 2, function(x){mean(x, na.rm=T)}) - apply(exprs[meta[,1] == ctrl.group, ,drop=F], 2, function(x){mean(x, na.rm=T)})
      linearFC <- 2^log2FC
    } else if (diff.space == "linear") {
      if (request.type == "biocrates") {
        diffMeans <- apply(exprs[meta[,1] == test.group, ,drop=F], 2, function(x){mean(x, na.rm=T)}) - apply(exprs[meta[,1] == ctrl.group, ,drop=F], 2, function(x){mean(x, na.rm=T)})
      } else {
        # here, data is already transformed to linear space in the pre processing step
        linearFC <- apply(exprs[meta[,1] == test.group, ,drop=F], 2, function(x){mean(x, na.rm=T)})/apply(exprs[meta[,1] == ctrl.group, ,drop=F], 2, function(x){mean(x, na.rm=T)})
        if (compute.log2fc) {
          log2FC <- log2(linearFC)
        }
      }
    }

    all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction

    if (compute.log2fc) {
      reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, linear_fc = linearFC, log2_fc = log2FC, row.names = NULL)
    } else {
      if (request.type == "biocrates") {
        reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, diff_means = diffMeans, row.names = NULL)
      } else {
        reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, linear_fc = linearFC, row.names = NULL)
      }
    }

    reportdf <- reportdf[order(reportdf$fdr), ,drop=F]
    write.xlsx(reportdf,
               paste0(outdir,"Report_",test,"_",comparison,".xlsx"),
               rowNames = F, overwrite = T)
    if (compute.log2fc) {
      fullreportdf <- rbind(c(NA,NA,NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID,drop=F])))
    } else {
      fullreportdf <- rbind(c(NA,NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID,drop=F])))
    }
    write.xlsx(fullreportdf,
               paste0(outdir,"FullReport_",test,"_",comparison,".xlsx"),
               rowNames = F, overwrite = T)
  } else if (test == "anova") {
    all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
                                                     summary(aov(x~meta[,1]))[[1]][["Pr(>F)"]][1],
                                                     1)}) # ANOVA test across study sites for each chemical
    all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction

    reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, row.names = NULL)
    reportdf <- reportdf[order(reportdf$fdr),,drop=F]
    write.xlsx(reportdf,
               paste0(outdir,"Report_",test,"_",comparison,".xlsx"),
               rowNames = F, overwrite = T)

    fullreportdf <- rbind(c(NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID,drop=F])))
    write.xlsx(fullreportdf,
               paste0(outdir,"FullReport_",test,"_",comparison,".xlsx"),
               rowNames = F, overwrite = T)

  } else if (test == "limma") {
    # Moderated t-test via limma lmFit/eBayes.
    # Data must be in log2 space — runModac ensures this before calling.
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
    library(limma)

    # Sanitize group names for use as design-matrix column names (handles hyphens, spaces, etc.)
    clean_ctrl <- make.names(ctrl.group)
    clean_test  <- make.names(test.group)
    group_factor <- factor(meta[, 1], levels = c(ctrl.group, test.group))
    design <- model.matrix(~0 + group_factor)
    colnames(design) <- c(clean_ctrl, clean_test)

    if (!is.null(pair_id)) {
      corfit <- duplicateCorrelation(t(exprs), design, block = pair_id)
      fit <- lmFit(t(exprs), design, block = pair_id, correlation = corfit$consensus)
    } else {
      fit <- lmFit(t(exprs), design)
    }
    cont <- makeContrasts(
      contrasts = paste0(clean_test, "-", clean_ctrl),
      levels    = design
    )
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- eBayes(fit2)

    # Align results to the original feature order
    all.pvals <- fit2$p.value[, 1][colnames(exprs)]
    log2FC    <- fit2$coefficients[, 1][colnames(exprs)]
    linearFC  <- 2^log2FC
    all.fdr   <- p.adjust(all.pvals, method = "BH")

    reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr,
                           linear_fc = linearFC, log2_fc = log2FC, row.names = NULL)
    reportdf <- reportdf[order(reportdf$fdr), , drop = FALSE]
    write.xlsx(reportdf,
               paste0(outdir, "Report_", test, "_", comparison, ".xlsx"),
               rowNames = FALSE, overwrite = TRUE)

    fullreportdf <- rbind(c(NA, NA, NA, NA, NA, meta[, 1]),
                          cbind(reportdf, t(exprs[, reportdf$ID, drop = FALSE])))
    write.xlsx(fullreportdf,
               paste0(outdir, "FullReport_", test, "_", comparison, ".xlsx"),
               rowNames = FALSE, overwrite = TRUE)

  } else if (test == "limma-full") {
    # Moderated t-test with the model fit on ALL groups supplied in meta.
    # eBayes variance shrinkage borrows information across all groups, not just the two
    # being contrasted. Data must already be in log2 space (caller's responsibility).
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
    library(limma)

    # Build design matrix across every group present in meta
    all_levels   <- unique(as.character(meta[, 1]))
    clean_names  <- make.names(all_levels)
    name_map     <- setNames(clean_names, all_levels)
    clean_ctrl   <- name_map[ctrl.group]
    clean_test   <- name_map[test.group]

    group_factor <- factor(meta[, 1], levels = all_levels)
    design <- model.matrix(~0 + group_factor)
    colnames(design) <- clean_names

    if (!is.null(pair_id)) {
      corfit <- duplicateCorrelation(t(exprs), design, block = pair_id)
      fit <- lmFit(t(exprs), design, block = pair_id, correlation = corfit$consensus)
    } else {
      fit <- lmFit(t(exprs), design)
    }
    cont <- makeContrasts(
      contrasts = paste0(clean_test, "-", clean_ctrl),
      levels    = design
    )
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- eBayes(fit2)

    # Align results to original feature order
    all.pvals <- fit2$p.value[, 1][colnames(exprs)]
    log2FC    <- fit2$coefficients[, 1][colnames(exprs)]
    linearFC  <- 2^log2FC
    all.fdr   <- p.adjust(all.pvals, method = "BH")

    reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr,
                           linear_fc = linearFC, log2_fc = log2FC, row.names = NULL)
    reportdf <- reportdf[order(reportdf$fdr), , drop = FALSE]
    write.xlsx(reportdf,
               paste0(outdir, "Report_", test, "_", comparison, ".xlsx"),
               rowNames = FALSE, overwrite = TRUE)

    fullreportdf <- rbind(c(NA, NA, NA, NA, NA, meta[, 1]),
                          cbind(reportdf, t(exprs[, reportdf$ID, drop = FALSE])))
    write.xlsx(fullreportdf,
               paste0(outdir, "FullReport_", test, "_", comparison, ".xlsx"),
               rowNames = FALSE, overwrite = TRUE)
  }

  pdf(paste0(outdir,"/boxplots_features_",comparison,".pdf"))
  for (i in 1:ncol(exprs)) {
    mytitle <- paste0(colnames(exprs)[i]," (FDR = ",round(all.fdr[i],2),")")
    par(cex.main=0.7)
    par(cex.axis=0.5)
    if (!(all(is.na(unlist(exprs[,i]))))) {
      boxplot(unlist(exprs[,i]) ~ factor(meta[,1], levels = group.ctrl.test),
              col = group.colors,
              xlab = NULL,
              ylab = NULL,
              main=mytitle,
              las=2,
              outline = F)
    }
  }
  dev.off()
}
