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
                       samplesAreRows = F) {
  
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
  
  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
  }
  
  createDir(outdir)
  
  if (test == "t-test") {
    test.group <- strsplit(comparison,"_over_")[[1]][1]
    ctrl.group <- strsplit(comparison,"_over_")[[1]][2]
    
    if ((length(test.group) == 0)|(length(ctrl.group) == 0)) {
      stop(cat("\n######## ERROR WITH COMPARISON NAME:", comparison,"########",
               "\n######## T-test comparison name format is 'test_over_control' ########\n"))
    }
    all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
                                                     t.test(x~meta[,1])$p.value,
                                                     1)})
    
    # code for paired t-test. make sure same number of samples in both groups
    # all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
    #                                                  t.test(x[meta[,1] == test.group],
    #                                                         x[meta[,1] == ctrl.group],
    #                                                         paired = T)$p.value,
    #                                                  1)})
    
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
  }
  
  # all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction
  # 
  # reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, row.names = NULL)
  # reportdf <- reportdf[order(reportdf$fdr),]
  # write.csv(reportdf, paste0(outdir,"/Report_",test,"_",comparison,".csv"), row.names = F)
  # 
  # fullreportdf <- rbind(c(NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID])))
  # write.csv(fullreportdf, paste0(outdir,"/FullReport_",test,"_",comparison,".csv"), row.names = F)
  
  pdf(paste0(outdir,"/boxplots_features_",comparison,".pdf"))
  for (i in 1:ncol(exprs)) {
    mytitle <- paste0(colnames(exprs)[i]," (FDR = ",round(all.fdr[i],2),")")
    par(cex.main=0.7)
    par(cex.axis=0.5)
    boxplot(unlist(exprs[,i]) ~ factor(meta[,1], levels = group.ctrl.test),
            col = group.colors,
            xlab = NULL,
            ylab = NULL,
            main=mytitle,
            las=2,
            outline = F)
  }
  dev.off()
}
