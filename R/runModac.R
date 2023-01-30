
runModac <- function(inputFile, comparisonsFile, type = "metabolomics", outdir = "results", scriptPath) {
  library(readxl)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  
  source(paste0(scriptPath,"/createDir.R"))
  source(paste0(scriptPath,"/preprocessMetab.R"))
  
  # Creating output dir for QA plots
  myoutdir <- paste(outdir,"QA_plots",sep = "/")
  createDir(myoutdir)
  
  # Creating output dir for report files
  myoutdir <- paste(outdir,"report",sep = "/")
  createDir(myoutdir)
  
  exprsdf <- preprocessMetab(inputFile = inputFile,
                             comparisonsFile = comparisonsFile,
                             outdir = outdir,
                             scriptPath = scriptPath)
  
  # Boxplot of all samples
  pdf(paste0(outdir,"/QA_plots/boxplot_all_samples.pdf"))
  mar.def <- par()$mar
  par(mar = mar.def + c(5,0,-3,0))
  par(cex.axis=0.6)
  boxplot(t(exprsdf[["norm"]]),
          ylab="Relative abundance",
          main="Comparison of samples (all methods)",
          las=2,
          outline = T)
  dev.off()
  par(mar = mar.def)
  
  # reading in all comparisons
  all.comparisons <- excel_sheets(comparisonsFile)[-1]
  
  source(paste0(scriptPath,"/myStatTest.R"))
  source(paste0(scriptPath,"/plotPCA.R"))
  for (i in seq_along(all.comparisons)) {
    # readying inputs
    all.info <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[i], n_max = 2, col_names = F), row.names = 1, check.rows = F))
    mycomparison <- all.info[1,]
    mytest <- all.info[2,]
    
    mygroups <- suppressMessages(as.character(read_excel(comparisonsFile, trim_ws = T, sheet = i+1, skip = 2, n_max = 1, col_names = F))[-1])
    mymeta <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[i], skip = 3), row.names = 2, check.rows = F))
    mymeta <- mymeta[mymeta[,1] %in% mygroups, , drop = F]
    myexprs.raw <- exprsdf[["raw"]][rownames(mymeta),]
    myexprs.norm <- exprsdf[["norm"]][rownames(mymeta),]
    
    # PCA plot
    plotPCA(exprs = myexprs.raw, meta = mymeta, outdir = paste0(outdir, "/pca"), suffix = paste0("_raw_",mycomparison), samplesAreRows = F)
    plotPCA(exprs = myexprs.norm, meta = mymeta, outdir = paste0(outdir, "/pca"), suffix = paste0("_norm_",mycomparison), samplesAreRows = F)
    
    # stat test
    myStatTest(exprs = myexprs.norm,
            meta = mymeta,
            test = mytest,
            comparison = mycomparison,
            outdir = paste0(outdir, "/report/"))
    
    # heatmaps
    settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 7, col_names = F), row.names = 1, check.rows = F))
    source(paste0(scriptPath,"/plotHeatmap.R"))
    plotHeatmap(exprs = myexprs.norm,
                meta = mymeta,
                test = mytest,
                comparison = mycomparison,
                outdir = paste0(outdir, "/heatmaps/"),
                groupOrder = mygroups,
                reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                cutoffStat = settings["statistic_selector",1],
                cutoff = settings["statistic_cutoff",1])
    plotHeatmap(exprs = myexprs.norm,
                meta = mymeta,
                test = mytest,
                comparison = mycomparison,
                outdir = paste0(outdir, "/heatmaps/"),
                groupOrder = mygroups,
                reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                cutoffStat = settings["statistic_selector",1],
                cutoff = 1)
  }
}
