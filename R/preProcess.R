preProcess <- function(inputFile, comparisonsFile, outdir, scriptPath, samplesAreRows = T) {
  
  all.methods <- excel_sheets(inputFile)
  settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 9, col_names = F), row.names = 1, check.rows = F))
  normalization <- settings["normalization",1]
  input_space <- settings["input_data_space",1]
  diff_space <- settings["differential_analysis_space",1]
  
  log2TF <- input_space == "linear" & diff_space == "log2"
  linearTF <- input_space == "log2" & diff_space == "linear"
  
  all.exprs.norm <- list()
  all.exprs.raw <- list()
  all.istd <- c()
  all.cv <- c()
  all.qc.num <- c()
  
  print(cat("##### Normalization type:", normalization, "#####\n"))
  print(cat("##### Input data space:", input_space, "#####\n"))
  print(cat("##### Differential analysis space:", diff_space, "#####\n"))
  
  for (i in seq_along(all.methods)) {
    mymethod <- all.methods[i]
    print(cat("##### Processing method:", mymethod,"#####\n"))
    
    myoutdir <- paste(outdir,"QA_plots",sep = "/")
    
    # reading in expression data per method
    inputdf <- suppressMessages(data.frame(read_excel(inputFile, trim_ws = T, sheet = mymethod, na = c("","N/A","NA")), row.names = 1, check.rows = F, check.names = F))
    
    if (!samplesAreRows) {
      inputdf <- data.frame(t(inputdf), check.rows = F, check.names = F)
    }
    
    # convert all columns to numeric
    inputdf[] <- lapply(inputdf, as.numeric)
    
    if(any(is.na(inputdf))) {
      print(cat("##### Replacing", sum(is.na(inputdf)),"NAs with 1 #####\n"))
      inputdf[is.na(inputdf)] <- 1 
    }
    
    if(log2TF & any(inputdf == 0)) {
      print(cat("##### Since data has to be log2 transformed, replacing", sum(inputdf == 0),"zeros with 1 #####\n"))
      inputdf[inputdf == 0] <- 1 
    }
    
    if (normalization == "none") {
      
      all.exprs.raw[[i]] <- inputdf
      
      # log2 transform
      if(log2TF){
        inputdf <- log2(inputdf)
      }
      
      #linear transform
      if(linearTF){
        inputdf <- 2^inputdf
      }
      
      all.exprs.norm[[i]] <- inputdf
      
      all.istd <- 0
      all.cv <- 0
      all.qc.num <- 0
    } else {
      # saving qc row start idx
      qc.idx <- as.numeric(settings["QC_row",1])-1
      if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
        all.exprs.raw[[i]] <- inputdf
      } else {
        all.exprs.raw[[i]] <- inputdf[-c(qc.idx:nrow(inputdf)),] 
      }
      
      if (normalization == "istd") {
        source(paste0(scriptPath,"/normISTD.R"))
        tempdf <- normISTD(inputdf = inputdf,
                           comparisonsFile = comparisonsFile,
                           myoutdir = myoutdir,
                           mymethod = mymethod)
        
        inputdf <- tempdf$data
        all.istd <- c(all.istd, tempdf$istd)
        all.cv <- c(all.cv, tempdf$cv)
        all.qc.num <- tempdf$qc_num
      }
      
      if (normalization == "iqr") {
        source(paste0(scriptPath,"/normIQR.R"))
        tempdf <- normIQR(inputdf = inputdf,
                          comparisonsFile = comparisonsFile)
        
        inputdf <- tempdf$data
        all.istd <- 0
        all.cv <- 0
        all.qc.num <- tempdf$qc_num
      }
      
      # Plotting QC samples
      
      # # OLD LIVER QC PLOTS CODE
      # temp <- pivot_longer(data.frame(feature=colnames(inputdf),t(inputdf[qc.idx:nrow(inputdf),])),cols = !feature, names_to = "sample")
      # myplot <- ggplot(temp, aes(x = feature, y = value, group = feature, colour = sample)) +
      #   geom_point(size = 2, alpha = 0.75) +
      #   scale_y_continuous(n.breaks = 25, labels = function(x) {
      #     scales::label_number_si(accuracy = 0.1)(x)
      #   }) +
      #   theme(axis.text.x = element_blank(),
      #         legend.position = "bottom",
      #         legend.box = "vertical",
      #         legend.margin = margin(),
      #         panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(),
      #         panel.background = element_blank(),
      #         axis.line = element_line(colour = "black")) +
      #   labs(title = paste0(mymethod, " quality control distribution"), x = 'Compounds', y = 'Quality control')
      # 
      # pdf(paste(myoutdir,paste0("qc_dist_",mymethod,".pdf"),sep = "/"))
      # print(myplot)
      # dev.off()
      # 
      # jpeg(paste(myoutdir,paste0("qc_dist_",mymethod,".jpg"),sep = "/"))
      # print(myplot)
      # dev.off()
      # 
      # write.csv(myplot$data,
      #           paste(myoutdir,paste0("qc_dist_",mymethod,".csv"),sep = "/"),
      #           quote = F,
      #           row.names = F)
      
      if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
        all.exprs.norm[[i]] <- inputdf
      } else {
        temp <- data.frame(feature=colnames(inputdf),t(inputdf[qc.idx:nrow(inputdf),]))
        
        mar.def <- par()$mar
        
        pdf(paste(myoutdir,paste0("qc_dist_",mymethod,".pdf"),sep = "/"))
        par(mar = mar.def + c(5,0,-3,0))
        boxplot(temp[,-1], ylab="Relative abundance", las = 2, cex.axis = 0.8)
        dev.off()
        
        jpeg(paste(myoutdir,paste0("qc_dist_",mymethod,".jpg"),sep = "/"))
        par(mar = mar.def + c(5,0,-3,0))
        boxplot(temp[,-1], ylab="Relative abundance", las = 2, cex.axis = 0.8)
        dev.off()
        
        par(mar = mar.def)
        dev.off()
        
        # write.table(temp,
        #             paste(myoutdir,paste0("qc_dist_",mymethod,".xls"),sep = "/"),
        #             sep = "\t",
        #             quote = F,
        #             row.names = F)
        write.xlsx(temp,
                   file.path(myoutdir,paste0("qc_dist_",mymethod,".xlsx")),
                   rowNames = F,
                   overwrite = T)
        
        all.exprs.norm[[i]] <- inputdf[-c(qc.idx:nrow(inputdf)),,drop=F] #dropping QC samples
      }
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
    
    jpeg(paste0(myoutdir,"/boxplot_",mymethod,".jpg"))
    par(mar = mar.def + c(5,0,-3,0))
    boxplot(t(all.exprs.norm[[i]]),
            ylab="Relative abundance",
            main=paste0("Comparison of samples (",mymethod,")"),
            las=2,
            outline = T)
    dev.off()
    
    par(mar = mar.def)
    dev.off()
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
  
  exprsdf.raw <- data.frame(exprsdf.raw, check.rows = F, check.names = F)
  exprsdf.norm <- data.frame(exprsdf.norm, check.rows = F, check.names = F)
  
  print(cat("##### Total number of samples:", nrow(exprsdf.norm), "#####\n"))
  print(cat("##### Total number of features:", ncol(exprsdf.norm), "#####\n"))
  
  myoutdir <- paste(outdir,"report",sep = "/")
  
  # write.table(data.frame(ID=rownames(exprsdf.raw), exprsdf.raw, check.names = F), paste0(myoutdir,"/raw_data.xls"), quote = F, row.names = F, sep = "\t", col.names = T)
  # write.table(data.frame(ID=rownames(exprsdf.norm), exprsdf.norm, check.names = F), paste0(myoutdir,"/normalized_data.xls"), quote = F, row.names = F, sep = "\t", col.names = T)
  write.xlsx(data.frame(ID=rownames(exprsdf.raw), exprsdf.raw, check.names = F), paste0(myoutdir,"/raw_data.xlsx"), rowNames = F, colNames = T, overwrite = T)
  write.xlsx(data.frame(ID=rownames(exprsdf.norm), exprsdf.norm, check.names = F), paste0(myoutdir,"/normalized_data.xlsx"), rowNames = F, colNames = T, overwrite = T)
  
  return(list(raw = exprsdf.raw,
              norm=exprsdf.norm,
              istd = all.istd,
              cv = all.cv,
              qc_num = all.qc.num,
              methods = all.methods))
}