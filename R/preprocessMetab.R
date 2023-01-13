preprocessMetab <- function(inputFile, comparisonsFile, outdir, scriptPath) {
  
  all.methods <- excel_sheets(inputFile)
  settings <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 7, col_names = F), row.names = 1, check.rows = F)
  normalization <- settings["normalization",1]
  log2tr <- as.logical(settings["log2transform",1])
  
  all.exprs.norm <- list()
  all.exprs.raw <- list()
  for (i in seq_along(all.methods)) {
    
    mymethod <- all.methods[i]
    myoutdir <- paste(outdir,"QA_plots",sep = "/")
    
    # reading in expression data per method
    inputdf <- data.frame(read_excel(inputFile, trim_ws = T, sheet = mymethod), row.names = 1, check.rows = F, check.names = F)
    
    if (normalization == "none") {
      # log2 transform
      if(log2tr){
        inputdf <- log2(inputdf)
      }
      all.exprs.raw[[i]] <- inputdf
      all.exprs.norm[[i]] <- inputdf
    } else {
      # saving qc row start idx
      qc.idx <- as.numeric(settings["QC_row",1])-1
      
      all.exprs.raw[[i]] <- inputdf[-c(qc.idx:nrow(inputdf)),]
      
      if (normalization == "istd") {
        source(paste0(scriptPath,"/normISTD.R"))
        inputdf <- normISTD(inputdf = inputdf,
                            comparisonsFile = comparisonsFile,
                            myoutdir = myoutdir,
                            mymethod = mymethod)
      }
      
      if (normalization == "iqr") {
        source(paste0(scriptPath,"/normIQR.R"))
        inputdf <- normIQR(inputdf = inputdf,
                           comparisonsFile = comparisonsFile)
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
    }
    
    pdf(paste0(myoutdir,"/boxplot_",mymethod,".pdf"))
    mar.def <- par()$mar
    par(mar = mar.def + c(5,0,-3,0))
    par(cex.axis=0.6)
    boxplot(t(all.exprs.norm[[i]]),
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
  
  myoutdir <- paste(outdir,"report",sep = "/")
  
  write.xlsx(exprsdf.raw, paste0(myoutdir,"/raw_data.xlsx"), overwrite = T, rowNames = T)
  write.xlsx(exprsdf.norm, paste0(myoutdir,"/normalized_data.xlsx"), overwrite = T, rowNames = T)
  
  return(list(raw=exprsdf.raw, norm=exprsdf.norm))
}