
runModac <- function(inputFile,
                     comparisonsFile,
                     type = "metabolomics",
                     outdir = "results",
                     project_title = "Project Report",
                     project_subtitle = "Metabolomics Analysis",
                     heatmap_color_scale = c("blue", "black", "yellow"),
                     scriptPath,
                     samplesAreRows = T,
                     sampleIDRow = 2) {
  
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
               outdir = myoutdir,
               geneNames = F,
               sampleIDRow = sampleIDRow)
      exprsdf <- preProcess(inputFile = paste0(myoutdir, "/aggregate_data.xlsx"),
                            comparisonsFile = comparisonsFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            samplesAreRows = T)
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
    colors_set <- c("red","steelblue","green", "orange", "pink","aquamarine","purple","grey","black","khaki","maroon","yellow")
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
      
      settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 8, col_names = F), row.names = 1, check.rows = F))
      fcType <- settings["fc_type",1]
      if (settings["log2transform",1]) {
        # stat test
        myStatTest(exprs = myexprs.norm,
                   meta = mymeta,
                   test = mytest,
                   comparison = mycomparison,
                   group.ctrl.test = mygroups,
                   group.colors = col_list$group[mygroups],
                   fc.type = "log2",
                   samplesAreRows = T,
                   outdir = paste0(outdir, "/report/"))
      } else {
        # stat test
        myStatTest(exprs = myexprs.norm,
                   meta = mymeta,
                   test = mytest,
                   comparison = mycomparison,
                   group.ctrl.test = mygroups,
                   group.colors = col_list$group[mygroups],
                   fc.type = fcType,
                   samplesAreRows = T,
                   outdir = paste0(outdir, "/report/")) 
      }
      
      if (type == "rppa") {
        fullreport <- read.csv(paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison,".csv"), header = T, row.names = NULL)
        tempdf <- read.xlsx(paste0(outdir,"/report/full_aggregate_data.xlsx"), rowNames = F)
        geneSymbols <- unname(sapply(fullreport$ID[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]}))
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
        
        if (fcType == "log2") {
          temp_fc_cutoff <- log2(as.numeric(settings["linear_fc_cutoff",1]))
        } else if (fcType == "linear") {
          temp_fc_cutoff <- settings["linear_fc_cutoff",1]
        }
        outFileName <- paste0(outdir,"/Signature/Sig_",mycomparison,"_",fcType,"FC",temp_fc_cutoff,"_",settings["padj_method",1],settings["padj_cutoff",1],".txt")
        
        source(paste0(scriptPath,"/createSignature.R"))
        createSignature(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".csv"),
                        outFile = outFileName,
                        fcCutoff = temp_fc_cutoff,
                        statType = settings["padj_method",1],
                        statCutoff = settings["padj_cutoff",1])
      }
      
      # volcano plots
      source(paste0(scriptPath,"/plotVolcano.R"))
      if (fcType == "log2" & mytest == "t-test") {
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
                  cutoff = settings["padj_cutoff",1],
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
  
  # Stop the parallel backend
  stopCluster(cl)
}
