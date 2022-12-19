myttest <- function(exprs, meta, comparison, outdir) {
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
    meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  }
  
  # mytest <- strsplit(comparison, "_over_")[[1]][1]
  # myctrl <- strsplit(comparison, "_over_")[[1]][2]
  # grouping <- meta[meta[,1] == mytest | meta[,1] == myctrl, , drop = F]
  mygroups <- unlist(strsplit(comparison, "_over_"))
  
  mygroups <- unlist(strsplit(comparison, "_"))
  if ("over" %in% mygroups) {
    mygroups <- mygroups[!(mygroups %in% "over")]
    mytest <- "t-test"
  } else {
    mytest <- "ANOVA"
  }
  
  grouping <- meta[meta[,1] %in% mygroups, , drop = F]
  myexprs <- exprs[rownames(grouping),]
  
  if (mytest == "t-test") {
    all.pvals <- apply(myexprs, 2, function(x) {t.test(x~grouping[,1])$p.value}) # ANOVA test across study sites for each chemical 
  } else if (mytest == "ANOVA") {
    all.pvals <- apply(myexprs, 2, function(x) {summary(aov(x~grouping[,1]))[[1]][["Pr(>F)"]][1]}) # ANOVA test across study sites for each chemical 
  }
  
  all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction
  
  reportdf <- data.frame(ID = colnames(myexprs), pval = all.pvals, fdr = all.fdr, row.names = NULL)
  reportdf <- reportdf[order(reportdf$fdr),]
  
  write.csv(reportdf, paste0(outdir,"/Report_",mytest,"_",comparison,".csv"), row.names = F)
  
  mycolors <- c("red","blue","green","yellow","black","grey","orange","turqoise")
  
  pdf(paste0(outdir,"/boxplots_features_",comparison,".pdf"))
  for (i in 1:ncol(myexprs)) {
    mytitle <- paste0(colnames(myexprs)[i]," (FDR = ",round(all.fdr[i],2),")")
    par(cex.main=0.7)
    par(cex.axis=0.5)
    boxplot(unlist(myexprs[,i]) ~ grouping[,1],
            col = mycolors,
            xlab = NULL,
            ylab = NULL,
            main=mytitle,
            las=2,
            outline = F)
  }
  dev.off()
}