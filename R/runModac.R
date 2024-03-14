
runModac <- function(inputFile,
                     comparisonsFile,
                     type,
                     outdir = "results",
                     project_title = "Report",
                     project_subtitle = "by",
                     heatmap_color_scale = c("blue", "black", "yellow"),
                     scriptPath,
                     samplesAreRows = T,
                     sampleIDRow = 2,
                     min_signal = NULL) {
  
  # Loading required packages and installing ones not present
  list.of.packages <- c("foreach","doParallel")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  # Set the number of cores to use
  num_cores <- 4

  # Register the parallel backend
  cl <- makeCluster(num_cores, outfile="")
  registerDoParallel(cl)

  # Convert the loop to parallel using foreach
  foreach(i = 1) %dopar% {
  # for (i in 1) {
    # Loading required packages and installing ones not present
    list.of.packages <- c("ggplot2", "readxl", "openxlsx","tidyr","foreach","doParallel")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
    
    source(paste0(scriptPath,"/createDir.R"))
    
    # Creating output dir for QA plots
    myoutdir <- paste(outdir,"QA_plots",sep = "/")
    createDir(myoutdir)
    
    # Creating output dir for report files
    myoutdir <- paste(outdir,"report",sep = "/")
    createDir(myoutdir)
    
    print(cat("##### Preprocessing started #####\n\n"))
    source(paste0(scriptPath,"/preProcess.R"))
    if (type == "rppa") {
      source(paste0(scriptPath,"/rppaAggr.R"))
      rppaAggr(inputFile = inputFile,
               outputFile = paste0(myoutdir, "/full_aggregate_data.xlsx"),
               sampleIDRow = sampleIDRow,
               cv_cutoff = 0.25,
               replace.cvcutoff = NULL,
               replace.na = NULL)
      
      tempdf <- read.xlsx(paste0(myoutdir, "/full_aggregate_data.xlsx"), sheet="Norm_Median")
      write.xlsx(list(rppa = data.frame(tempdf[,-c(1,3,4)],check.rows = F,check.names = F)),
                 paste0(myoutdir, "/aggregate_data.xlsx"),
                 rowNames = F)
      exprsdf <- preProcess(inputFile = paste0(myoutdir, "/aggregate_data.xlsx"),
                            comparisonsFile = comparisonsFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            samplesAreRows = F)
    } else {
      exprsdf <- preProcess(inputFile = inputFile,
                            comparisonsFile = comparisonsFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            samplesAreRows = samplesAreRows)
    }
    print(cat("##### Preprocessing complete #####\n\n"))
    
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
    
    jpeg(paste0(outdir,"/QA_plots/boxplot_all_samples.jpg"))
    par(mar = mar.def + c(5,0,-3,0))
    boxplot(t(exprsdf[["norm"]]),
            ylab="Relative abundance",
            main="Comparison of samples (all methods)",
            las=2,
            outline = T)
    dev.off()
    
    par(mar = mar.def)
    dev.off()
    
    # reading in all comparisons
    all.comparisons <- excel_sheets(comparisonsFile)[-1]
    all.comparison.labels <- c()
    
    # assigning group colors for heatmaps
    all.groups <- c()
    for (i in seq_along(all.comparisons)) {
      mymeta <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[i], skip = 3), row.names = 2, check.rows = F))
      all.groups <- c(all.groups, mymeta[,1])
    }
    all.groups <- unique(all.groups)
    
    # This color set restricted colors to 12 groups
    # colors_set <- c("red","steelblue","green", "orange", "pink","aquamarine","purple","grey","black","khaki","maroon","yellow")
    
    # Expanding color set to color 36 groups
    # library(Polychrome)
    # unname(createPalette(36,  c("#ff0000", "#00ff00", "#0000ff")))
    colors_set <- c("#FD001C", "#16FF00", "#0000FF", "#63422A", "#FE16D4", "#00F7FF",
                    "#F1E500", "#D7B0FD", "#FD8E3D", "#008E3B", "#CD1660", "#000099",
                    "#CB0DF7", "#225C75", "#FC9DBD", "#22AAFE", "#920D83", "#FCCB95",
                    "#512EA3", "#ACF25F", "#00FFC4", "#94352A", "#858B22", "#A2D1F5",
                    "#FC83FD", "#946A8E", "#22A593", "#F2532A", "#FDBB0D", "#D63B9A",
                    "#ABEFAE", "#B37AFF", "#FCD5F2", "#005332", "#1CE26D", "#FC7C7D")
    col_list <- list(group=colors_set[seq_along(all.groups)])
    names(col_list$group) <- all.groups
    
    source(paste0(scriptPath,"/myStatTest.R"))
    source(paste0(scriptPath,"/plotPCA.R"))
    for (i in seq_along(all.comparisons)) {
      # readying inputs
      all.info <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[i], n_max = 2, col_names = F), row.names = 1, check.rows = F))
      mycomparison <- all.info[1,]
      all.comparison.labels <- c(all.comparison.labels, mycomparison)
      mytest <- all.info[2,]
      
      print(cat("##### Performing", mytest, "for comparison:", mycomparison, "#####\n"))
      
      mygroups <- suppressMessages(as.character(read_excel(comparisonsFile, trim_ws = T, sheet = i+1, skip = 2, n_max = 1, col_names = F))[-1])
      mymeta <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[i], skip = 3), row.names = 2, check.rows = F))
      mymeta <- mymeta[mymeta[,1] %in% mygroups, , drop = F]
      myexprs.raw <- exprsdf[["raw"]][rownames(mymeta),]
      myexprs.norm <- exprsdf[["norm"]][rownames(mymeta),]
      
      # filtering out features for which none of the samples in this comparison have signal > min_signal
      if (!(is.null(min_signal))) {
        print(cat("##### Filtering out features by min_signal cutoff =", min_signal,"#####\n"))
        if (min_signal == 0) {
          keep <- apply(myexprs.raw,2,function(x){any(x>min_signal)})
        } else {
          keep <- apply(myexprs.raw,2,function(x){any(x>=min_signal)}) 
        }
        print(cat("##### Dropping", length(keep)-sum(keep),"out of the", length(keep)," total features #####\n"))
        myexprs.raw <- myexprs.raw[,keep,drop=F]
        myexprs.norm <- myexprs.norm[,keep,drop=F] 
      }
      
      # PCA plot
      plotPCA(exprs = myexprs.raw,
              meta = mymeta,
              outdir = paste0(outdir, "/pca"),
              suffix = paste0("_raw_",mycomparison),
              groupColors = col_list$group[mygroups],
              samplesAreRows = T)
      plotPCA(exprs = myexprs.norm,
              meta = mymeta,
              outdir = paste0(outdir, "/pca"),
              suffix = paste0("_norm_",mycomparison),
              groupColors = col_list$group[mygroups],
              samplesAreRows = T)
      
      settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 9, col_names = F), row.names = 1, check.rows = F))
      input_space <- settings["input_data_space",1]
      diff_space <- settings["differential_analysis_space",1]
      compute_log2fc <- as.logical(settings["compute_log2fc",1])
      
      myStatTest(exprs = myexprs.norm,
                 meta = mymeta,
                 test = mytest,
                 comparison = mycomparison,
                 group.ctrl.test = mygroups,
                 group.colors = col_list$group[mygroups],
                 diff.space = diff_space,
                 compute.log2fc = compute_log2fc,
                 samplesAreRows = T,
                 outdir = paste0(outdir, "/report/"))
      
      if (type == "rppa") {
        fullreport <- read.csv(paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison,".csv"), header = T, row.names = NULL)
        tempdf <- read.xlsx(paste0(outdir,"/report/full_aggregate_data.xlsx"), sheet = "Norm_Median", rowNames = F)
        geneSymbols <- unlist(unname(sapply(fullreport$ID[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]})))
        finalreport <- data.frame(GeneSymbol = c(NA,geneSymbols), fullreport)
        names(finalreport)[2] <- "AB_name"
        write.csv(finalreport, paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison,".csv"), row.names = F) 
      }
      
      if (mytest == "t-test") {
        createDir(paste0(outdir,"/Rnk"))
        createDir(paste0(outdir,"/Signature"))
        
        source(paste0(scriptPath,"/createRNK.R"))
        createRNK(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                  outFile = paste0(outdir,"/Rnk/Rnk_",mycomparison,".rnk"))
        
        if (compute_log2fc) {
          temp_fc_cutoff <- log2(as.numeric(settings["linear_fc_cutoff",1]))
          outFileName <- paste0(mycomparison,"_log2FC",round(temp_fc_cutoff,2),"_",settings["padj_method",1],settings["padj_cutoff",1])
        } else {
          temp_fc_cutoff <- as.numeric(settings["linear_fc_cutoff",1])
          outFileName <- paste0(mycomparison,"_linearFC",round(temp_fc_cutoff,2),"_",settings["padj_method",1],settings["padj_cutoff",1])
        }
        
        source(paste0(scriptPath,"/createSignature.R"))
        createSignature(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                        outFileName = outFileName,
                        outDir = paste0(outdir,"/Signature/"),
                        fcCutoff = temp_fc_cutoff,
                        statType = settings["padj_method",1],
                        statCutoff = as.numeric(settings["padj_cutoff",1]))
        
        if (type == "rppa") {
          # For RPPA analysis, saving 1 signature file with antibody names and another with gene symbols
          sigdf <- read_tsv(paste0(outdir,"/Signature/Sig_",outFileName,".txt"),
                            col_names = F,
                            show_col_types = FALSE)
          tempName <- paste0(outFileName,"-AB")
          tempFileName <- paste0(outdir,"/Signature/Sig_",tempName,".txt")
          sigdf1 <- sigdf
          sigdf1[1,1] <- tempName
          write_tsv(x = sigdf1,
                    file = tempFileName,
                    col_names = FALSE,
                    na = "")
          
          # sigdf[1,1] <- tempName
          
          sigdf$X1[-1] <- unname(sapply(sigdf$X1[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]}))
          
          tempFileName <- paste0(outdir,"/Signature/Sig_",outFileName,".txt")
          write_tsv(x = sigdf,
                    file = tempFileName,
                    col_names = FALSE,
                    na = "")
          
        }
      }
      
      # volcano plots
      # if compute_log2fc is F, volcano plots should not be generated and a message should be displayed
      source(paste0(scriptPath,"/plotVolcano.R"))
      if (compute_log2fc & mytest == "t-test") {
        plotVolcano(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                    myComparison = mycomparison,
                    outDir = paste0(outdir, "/volcano_plots/"),
                    fcCutoff = as.numeric(settings["linear_fc_cutoff",1]),
                    padjCutoff = as.numeric(settings["padj_cutoff",1]),
                    padjMethod = settings["padj_method",1])
      }
      # heatmaps
      source(paste0(scriptPath,"/plotHeatmap.R"))
      plotHeatmap(exprs = myexprs.norm,
                  meta = mymeta,
                  test = mytest,
                  comparison = mycomparison,
                  outdir = paste0(outdir, "/heatmaps/"),
                  groupOrder = mygroups,
                  groupColors = col_list$group[mygroups],
                  heatmapColorScale = heatmap_color_scale,
                  reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                  cutoffStat = settings["padj_method",1],
                  cutoff = as.numeric(settings["padj_cutoff",1]),
                  samplesAreRows = T)
      plotHeatmap(exprs = myexprs.norm,
                  meta = mymeta,
                  test = mytest,
                  comparison = mycomparison,
                  outdir = paste0(outdir, "/heatmaps/"),
                  groupOrder = mygroups,
                  groupColors = col_list$group[mygroups],
                  heatmapColorScale = heatmap_color_scale,
                  reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                  cutoffStat = settings["padj_method",1],
                  cutoff = 1,
                  samplesAreRows = T)
    }
    
    source(paste0(scriptPath,"/createPptx.R"))
    template_pptx_path <- paste0(gsub("/R$|/R/$","/data",scriptPath),"/Project_Report_Template.pptx")
    createPptx(project_title = project_title,
               project_subtitle = project_subtitle,
               exprsdf = exprsdf,
               settings = settings,
               all.comparison.labels = all.comparison.labels,
               outdir = outdir,
               template_pptx_path = template_pptx_path)
    
    writeLines(capture.output(sessionInfo()), paste0(outdir,"/sessionInfo.txt"))
  }
  
  # Check if "Rplots.pdf" exists in the current directory
  file_path <- file.path(getwd(), "Rplots.pdf")
  
  # If the file exists, delete it
  if (file.exists(file_path)) {
    file.remove(file_path)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
}
