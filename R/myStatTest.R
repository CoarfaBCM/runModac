myStatTest <- function(exprs, meta, test, comparison, outdir) {
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
    meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  }
  
  if (test == "t-test") {
    all.pvals <- apply(exprs, 2, function(x) {t.test(x~meta[,1])$p.value}) # ANOVA test across study sites for each chemical 
  } else if (test == "anova") {
    all.pvals <- apply(exprs, 2, function(x) {summary(aov(x~meta[,1]))[[1]][["Pr(>F)"]][1]}) # ANOVA test across study sites for each chemical 
  }
  
  all.fdr <- p.adjust(all.pvals, method = "BH") # FDR correction
  
  reportdf <- data.frame(ID = colnames(exprs), pval = all.pvals, fdr = all.fdr, row.names = NULL)
  reportdf <- reportdf[order(reportdf$fdr),]
  
  write.csv(reportdf, paste0(outdir,"/Report_",test,"_",comparison,".csv"), row.names = F)
  
  mycolors <- c("red","blue","green","yellow","black","grey","orange","turqoise")
  
  pdf(paste0(outdir,"/boxplots_features_",comparison,".pdf"))
  for (i in 1:ncol(exprs)) {
    mytitle <- paste0(colnames(exprs)[i]," (FDR = ",round(all.fdr[i],2),")")
    par(cex.main=0.7)
    par(cex.axis=0.5)
    boxplot(unlist(exprs[,i]) ~ meta[,1],
            col = mycolors,
            xlab = NULL,
            ylab = NULL,
            main=mytitle,
            las=2,
            outline = F)
  }
  dev.off()
}