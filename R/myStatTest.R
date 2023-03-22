myStatTest <- function(exprs, meta, test, comparison, group.ctrl.test, group.colors, outdir, fc.type, samplesAreRows = F) {
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
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
    all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
                                                     t.test(x~meta[,1])$p.value,
                                                     1)})
    
    if (fc.type == "log2") {
      all.fc <- apply(exprs[meta[,1] == group.ctrl.test[2], ], 2, mean) - apply(exprs[meta[,1] == group.ctrl.test[1], ], 2, mean)
    } else if (fc.type == "linear") {
      all.fc <- apply(exprs[meta[,1] == group.ctrl.test[2], ], 2, mean)/apply(exprs[meta[,1] == group.ctrl.test[1], ], 2, mean)
    }
    
    all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction
    
    reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, fc = all.fc, row.names = NULL)
    names(reportdf)[names(reportdf) == "fc"] <- paste0(fc.type,"_fc")
    reportdf <- reportdf[order(reportdf$fdr),]
    write.csv(reportdf, paste0(outdir,"/Report_",test,"_",comparison,".csv"), row.names = F)
    
    fullreportdf <- rbind(c(NA,NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID])))
    write.csv(fullreportdf, paste0(outdir,"/FullReport_",test,"_",comparison,".csv"), row.names = F)
  } else if (test == "anova") {
    all.pvals <- apply(exprs, 2, function(x) {ifelse(var(x) > 0,
                                                     summary(aov(x~meta[,1]))[[1]][["Pr(>F)"]][1],
                                                     1)}) # ANOVA test across study sites for each chemical
    all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction
    
    reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, row.names = NULL)
    reportdf <- reportdf[order(reportdf$fdr),]
    write.csv(reportdf, paste0(outdir,"/Report_",test,"_",comparison,".csv"), row.names = F)
    
    fullreportdf <- rbind(c(NA,NA,NA,meta[,1]),cbind(reportdf, t(exprs[,reportdf$ID])))
    write.csv(fullreportdf, paste0(outdir,"/FullReport_",test,"_",comparison,".csv"), row.names = F)
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