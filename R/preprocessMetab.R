preprocessMetab <- function(inputFile, comparisonsFile, outdir, scriptPath) {
  
  all.methods <- excel_sheets(inputFile)
  settings <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 8, col_names = F), row.names = 1, check.rows = F)
  normalization <- settings["normalization",1]
  log2tr <- as.logical(settings["log2transform",1])
  
  all.exprs.norm <- list()
  all.exprs.raw <- list()
  for (i in seq_along(all.methods)) {
    
    mymethod <- all.methods[i]
    
    # reading in expression data per method
    inputdf <- data.frame(read_excel(inputFile, trim_ws = T, sheet = mymethod), row.names = 1, check.rows = F, check.names = F)
    
    # saving qc row start idx
    qc.idx <- as.numeric(settings["QC_row",1])-1
    
    all.exprs.raw[[i]] <- inputdf[-c(qc.idx:nrow(inputdf)),]
    
    # Creating output dir for QA plots
    myoutdir <- paste(outdir,"QA_plots",sep = "/")
    createDir(myoutdir)
    
    if (normalization == "istd") {
      source(paste0(scriptPath,"/normISTD.R"))
      inputdf <- normISTD(inputdf = inputdf,
                          comparisonsFile = comparisonsFile,
                          myoutdir = myoutdir,
                          mymethod = mymethod)
    }
    
    if (log2tr) {
      # Log2 transformation
      inputdf <- log2(inputdf) 
    }
    
    # Plotting liver QC samples
    temp <- pivot_longer(data.frame(feature=colnames(inputdf),t(inputdf[qc.idx:nrow(inputdf),])),cols = !feature, names_to = "sample")
    pdf(paste(myoutdir,paste0("liver_dist_",mymethod,".pdf"),sep = "/"))
    print(ggplot(temp, aes(x = feature, y = value, group = feature, colour = sample)) +
            geom_point(size = 2, alpha = 0.75) +
            scale_y_continuous(n.breaks = 25, labels = function(x) {
              scales::label_number_si(accuracy = 0.1)(x)
            }) +
            theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
            labs(title = paste0(mymethod, " liver control distribution"), x = 'Compounds', y = 'Liver control'))
    dev.off()
    
    all.exprs.norm[[i]] <- inputdf[-c(qc.idx:nrow(inputdf)),] #dropping QC samples
    
    pdf(paste0(outdir,"/QA_plots/boxplot_",mymethod,".pdf"))
    mar.def <- par()$mar
    par(mar = mar.def + c(5,0,-3,0))
    par(cex.axis=0.6)
    boxplot(t(all.exprs.norm[[1]]),
            ylab="Relative abundance",
            main=paste0("Comparison of samples (",mymethod,")"),
            las=2,
            outline = T)
    dev.off()
    par(mar = mar.def)
  }
  
  # merging all methods
  exprsdf.raw <- Reduce(function(dtf1, dtf2) {
    ans <- merge(dtf1, dtf2, by = "row.names", all.x = TRUE)
    row.names(ans) <- ans[,"Row.names"]
    ans[,!names(ans) %in% "Row.names"]
  },
  all.exprs.raw)
  
  exprsdf.norm <- Reduce(function(dtf1, dtf2) {
    ans <- merge(dtf1, dtf2, by = "row.names", all.x = TRUE)
    row.names(ans) <- ans[,"Row.names"]
    ans[,!names(ans) %in% "Row.names"]
  },
  all.exprs.norm)
  
  pdf(paste0(outdir,"/QA_plots/boxplot_all_samples.pdf"))
  mar.def <- par()$mar
  par(mar = mar.def + c(5,0,-3,0))
  par(cex.axis=0.6)
  boxplot(t(exprsdf.norm),
          ylab="Relative abundance",
          main="Comparison of samples (all methods)",
          las=2,
          outline = T)
  dev.off()
  par(mar = mar.def)
  
  myoutdir <- paste(outdir,"report",sep = "/")
  createDir(myoutdir)
  
  write.xlsx(exprsdf.raw, paste0(myoutdir,"/raw_data.xlsx"), overwrite = T, rowNames = T)
  write.xlsx(exprsdf.norm, paste0(myoutdir,"/normalized_data.xlsx"), overwrite = T, rowNames = T)
  
  return(list(raw=exprsdf.raw, norm=exprsdf.norm))
}