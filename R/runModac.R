
runModac <- function(inputFile, # path to input data Excel file
                     comparisonsFile, # path to Excel file defining comparisons, groups, and settings
                     type, # analysis type: "metabolomics", "rppa", "biocrates", or "none"
                     outdir = "results", # directory where all outputs are written
                     project_title = "Report", # title shown on the PPTX report
                     project_subtitle = "by", # subtitle shown on the PPTX report
                     heatmap_color_scale = c("blue", "black", "yellow"), # low/mid/high colors for heatmaps
                     samplesAreRows = T, # TRUE if input sheets have samples as rows (FALSE if transposed)
                     sampleIDRow = 2, # which metadata row/column holds sample IDs (used for RPPA aggregation)
                     max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
                     min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
                     replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
                     replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
                     onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
                     scriptPath # directory containing the pipeline's R scripts
                     ) {
  
  # Loading required packages and installing ones not present
  # list.of.packages <- c("foreach","doParallel")
  # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  # Creating output root directory
  source(paste0(scriptPath,"/createDir.R"))
  createDir(outdir)
  
  # Start log
  logfilename <- paste0("runMODAC.log.",format(Sys.time(), "%m%d%y_%H%M%S"),".txt")
  sink(file = paste(outdir, logfilename, sep = "/"), split = TRUE)
  
  # Creating function to install and load packages
  install_pkg <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Setting CRAN mirror
  options(repos = c(CRAN = "https://cran.r-project.org"))
  
  # Installing required packages
  list.of.packages <- c("foreach","doParallel")
  lapply(list.of.packages, install_pkg)
  
  # # Set the number of cores to use
  # num_cores <- 4
  # 
  # # Register the parallel backend
  # cl <- makeCluster(num_cores, outfile="")
  # registerDoParallel(cl)
  # 
  # # Convert the loop to parallel using foreach
  # foreach(i = 1) %dopar% {
  for (i in 1) {
    # Loading required packages and installing ones not present
    # list.of.packages <- c("ggplot2", "readxl", "openxlsx","tidyr","foreach","doParallel")
    # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    # if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
    
    # setting CRAN mirror
    options(repos = c(CRAN = "https://cran.r-project.org"))
    
    # installing required packages
    list.of.packages <- c("ggplot2", "readxl", "openxlsx","tidyr","dplyr","foreach","doParallel")
    lapply(list.of.packages, install_pkg)
    
    # installing complex heatmap
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install(c("ComplexHeatmap"))
    
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
                 rowNames = F, overwrite = T)
      exprsdf <- preProcess(inputFile = paste0(myoutdir, "/aggregate_data.xlsx"),
                            comparisonsFile = comparisonsFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            samplesAreRows = F,
                            replaceNA = replaceNA,
                            replaceZeros = replaceZeros)
    } else {
      exprsdf <- preProcess(inputFile = inputFile,
                            comparisonsFile = comparisonsFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            samplesAreRows = samplesAreRows,
                            replaceNA = replaceNA,
                            replaceZeros = replaceZeros)
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
    
    # stop here if only normalization requested
    if (onlyNorm) {
      print("##### Normalization complete. No further downstream analysis requested. #####")
    } else {
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

    # collect per-comparison dropped-feature reports (written to one workbook after the loop)
    dropped_reports <- list()

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
      
      commonsamples <- intersect(rownames(exprsdf[["raw"]]), rownames(mymeta))
      print(cat("##### Number of samples in current comparison:", nrow(mymeta), "#####\n"))
      print(cat("##### Number of samples in input data:", nrow(exprsdf[["raw"]]), "#####\n"))
      print(cat("##### Number of samples common between current comparison and input data:", length(commonsamples), "#####\n"))
      
      if (nrow(mymeta) > length(commonsamples)){
        print(cat("##### Number of samples present in current comparison but not in input:", nrow(mymeta)-length(commonsamples), "#####\n"))
        print(cat("##### Samples present in current comparison but not in input:", rownames(mymeta)[!(rownames(mymeta) %in% rownames(exprsdf[["raw"]]))], "#####\n"))
        print(cat("##### The above mentioned non matching samples are dropped from analysis of this comparison #####\n"))
      }
      
      if (nrow(exprsdf[["raw"]]) > length(commonsamples)){
        print(cat("##### Number of samples present in input but not in current comparison:", nrow(exprsdf[["raw"]])-length(commonsamples), "#####\n"))
        print(cat("##### Samples present in input but not in current comparison:", rownames(exprsdf[["raw"]])[!(rownames(exprsdf[["raw"]]) %in% rownames(mymeta))], "#####\n"))
        print(cat("##### The above mentioned non matching samples are dropped from analysis of this comparison #####\n"))
      }
      
      myexprs.raw <- exprsdf[["raw"]][commonsamples,, drop=F]
      myexprs.norm <- exprsdf[["norm"]][commonsamples,, drop=F]
      
      # filtering out features for which none of the samples in this comparison have signal > min_signal
      min_signal_report <- NULL
      if (!(is.null(min_signal))) {
        print(cat("##### Filtering out features by min_signal cutoff =", min_signal,"#####\n"))
        if (min_signal == 0) {
          keep <- apply(myexprs.raw,2,function(x){any(x>min_signal, na.rm = TRUE)})
        } else {
          keep <- apply(myexprs.raw,2,function(x){any(x>=min_signal, na.rm = TRUE)})
        }
        print(cat("##### Dropping", length(keep)-sum(keep),"out of the", length(keep)," total features #####\n"))
        dropped_ms <- colnames(myexprs.raw)[!keep]
        if (length(dropped_ms) > 0) {
          min_signal_report <- data.frame(feature = dropped_ms,
                                          reason  = paste0("max signal < min_signal (", min_signal, ")"),
                                          row.names = NULL, check.names = F)
        }
        myexprs.raw <- myexprs.raw[,keep,drop=F]
        myexprs.norm <- myexprs.norm[,keep,drop=F]
      }
      
      dropZeroVarArrays <- function(inputdf, checkRows = F, checkCols = T){
        if (checkCols) {
          flag <- colSums(inputdf, na.rm = T) == 0 | colSums(!is.na(inputdf)) == 0
          if (any(flag, na.rm = T)) {
            print(cat("##### Dropping features with no expression or missing values across all samples for this comparison #####\n"))
            print(cat("##### Number of features dropped:", sum(flag), "#####\n"))
            print(cat("##### Exact features dropped:", colnames(inputdf)[flag], "#####\n"))
            inputdf <- inputdf[,!flag, drop=F]
          }
        }
        
        if (checkRows) {
          flag <- rowSums(inputdf, na.rm = T) == 0 | rowSums(!is.na(inputdf)) == 0
          if (any(flag, na.rm = T)) {
            print(cat("##### Dropping samples with no expression or missing values across all features for this comparison #####\n"))
            print(cat("##### Number of samples dropped:", sum(flag), "#####\n"))
            print(cat("##### Exact samples dropped:", rownames(inputdf)[flag], "#####\n"))
            inputdf <- inputdf[!flag,, drop=F]
          }
        }
        
        return(inputdf)
      }
      
      myexprs.raw <- dropZeroVarArrays(myexprs.raw, checkRows = T, checkCols = F)
      myexprs.norm <- dropZeroVarArrays(myexprs.norm, checkRows = T, checkCols = F)

      # Drop features too sparse within any comparison group:
      #  - missing fraction strictly greater than 'threshold' (the `max_missing_frac` arg) in any group, OR
      #  - fewer than 'minNonMissing' non-missing samples in any group (always enforced, even when
      #    max_missing_frac = 1) -- this is what keeps t.test/aov from failing on a near-empty group.
      dropByMissingness <- function(exprs, meta, groups, threshold, minNonMissing = 2){
        drop   <- rep(FALSE, ncol(exprs))
        reason <- rep(NA_character_, ncol(exprs))
        miss.by.group <- list()
        nok.by.group  <- list()
        for (g in groups) {
          samps <- intersect(rownames(meta)[meta[,1] == g], rownames(exprs))
          sub   <- exprs[samps, , drop = F]
          miss  <- colMeans(is.na(sub))    # fraction missing per feature in this group
          nok   <- colSums(!is.na(sub))     # non-missing count per feature in this group
          miss.by.group[[g]] <- miss
          nok.by.group[[g]]  <- nok
          fail  <- (miss > threshold) | (nok < minNonMissing)
          new   <- fail & is.na(reason)     # keep the first (most upstream) failing reason
          reason[new] <- ifelse(nok[new] < minNonMissing,
                                paste0("<", minNonMissing, " non-missing in group '", g, "'"),
                                paste0(">", round(threshold*100), "% missing in group '", g, "'"))
          drop  <- drop | fail
        }
        report <- NULL
        if (any(drop)) {
          report <- data.frame(feature = colnames(exprs)[drop],
                               reason  = reason[drop],
                               row.names = NULL, check.names = F)
          for (g in groups) {
            report[[paste0("missing_frac_", g)]] <- round(miss.by.group[[g]][drop], 3)
            report[[paste0("n_nonmissing_", g)]] <- nok.by.group[[g]][drop]
          }
          print(cat("##### Missingness filter: dropping", sum(drop), "features for this comparison (threshold =", threshold, ", min non-missing per group =", minNonMissing, ") #####\n"))
          print(cat("##### Exact features dropped:", colnames(exprs)[drop], "#####\n"))
        }
        return(list(keep = !drop, report = report))
      }

      # Apply the per-group missingness filter. Compute the mask on the normalized
      # (analysis) matrix, then drop the same features from raw so PCA/stats/heatmaps
      # stay on a consistent feature set.
      missfilt <- dropByMissingness(myexprs.norm, mymeta, mygroups, threshold = max_missing_frac)
      myexprs.norm <- myexprs.norm[, missfilt$keep, drop = F]
      if (!is.null(min_signal_report) || !is.null(missfilt$report)) {
        dropped_reports[[mycomparison]] <- bind_rows(min_signal_report, missfilt$report)
      }
      
      # PCA plot
      if (ncol(myexprs.raw) > 1) {
        plotPCA(exprs = myexprs.raw,
                meta = mymeta,
                outdir = paste0(outdir, "/pca"),
                suffix = paste0("_raw_",mycomparison),
                groupColors = col_list$group[mygroups],
                samplesAreRows = T)
      }
      
      if (ncol(myexprs.norm) > 1) {
        plotPCA(exprs = myexprs.norm,
                meta = mymeta,
                outdir = paste0(outdir, "/pca"),
                suffix = paste0("_norm_",mycomparison),
                groupColors = col_list$group[mygroups],
                samplesAreRows = T)
      }
      
      settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 9, col_names = F), row.names = 1, check.rows = F))
      input_space <- settings["input_data_space",1]
      diff_space <- settings["differential_analysis_space",1]
      compute_log2fc <- as.logical(settings["compute_log2fc",1])
      
      myStatTest(exprs = myexprs.norm,
                 meta = mymeta,
                 test = mytest,
                 request.type = type,
                 comparison = mycomparison,
                 group.ctrl.test = mygroups,
                 group.colors = col_list$group[mygroups],
                 diff.space = diff_space,
                 compute.log2fc = compute_log2fc,
                 samplesAreRows = T,
                 outdir = paste0(outdir, "/report/"))
      
      if (type == "rppa") {
        fullreport <- read.xlsx(paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison,".xlsx"), colNames = T, rowNames = F)
        tempdf <- read.xlsx(paste0(outdir,"/report/full_aggregate_data.xlsx"), sheet = "Norm_Median", rowNames = F)
        geneSymbols <- unlist(unname(sapply(fullreport$ID[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]})))
        finalreport <- data.frame(GeneSymbol = c(NA,geneSymbols), fullreport)
        names(finalreport)[2] <- "AB_name"
        write.xlsx(finalreport,
                   paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison,".xlsx"),
                   rowNames = F, overwrite = T)
      }
      
      if (mytest == "t-test") {
        createDir(paste0(outdir,"/Rnk"))
        createDir(paste0(outdir,"/Signature"))
        
        source(paste0(scriptPath,"/createRNK.R"))
        createRNK(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".xlsx"),
                  outFile = paste0(outdir,"/Rnk/Rnk_",mycomparison,".rnk"))
        
        if (type == "biocrates") {
          cat("#### No Signatures for Biocrates analysis since we are using difference of means and not linear or log2 fold change. #### \n")
        } else {
          if (compute_log2fc) {
            temp_fc_cutoff <- log2(as.numeric(settings["linear_fc_cutoff",1]))
            outFileName <- paste0(mycomparison,"_log2FC",round(temp_fc_cutoff,2),"_",settings["padj_method",1],settings["padj_cutoff",1])
          } else {
            temp_fc_cutoff <- as.numeric(settings["linear_fc_cutoff",1])
            outFileName <- paste0(mycomparison,"_linearFC",round(temp_fc_cutoff,2),"_",settings["padj_method",1],settings["padj_cutoff",1])
          }
          
          source(paste0(scriptPath,"/createSignature.R"))
          createSignature(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".xlsx"),
                          outFileName = outFileName,
                          outDir = paste0(outdir,"/Signature/"),
                          fcCutoff = temp_fc_cutoff,
                          statType = settings["padj_method",1],
                          statCutoff = as.numeric(settings["padj_cutoff",1]))
        }
        
        if (type == "rppa") {
          # For RPPA analysis, saving 1 signature file with antibody names and another with gene symbols
          sigdf <- read_tsv(paste0(outdir,"/Signature/sig.",outFileName,".txt"),
                            col_names = F,
                            show_col_types = FALSE)
          tempName <- paste0(outFileName,"-AB")
          tempFileName <- paste0(outdir,"/Signature/sig.",tempName,".txt")
          sigdf1 <- sigdf
          sigdf1[1,1] <- tempName
          write_tsv(x = sigdf1,
                    file = tempFileName,
                    col_names = FALSE)
          
          # sigdf[1,1] <- tempName
          
          sigdf$X1[-1] <- unname(sapply(sigdf$X1[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]}))
          sigdf <- sigdf[!is.na(sigdf[,1]), ,drop=F] # drop rows where antibody gene symbol is NA
          tempFileName <- paste0(outdir,"/Signature/sig.",outFileName,".txt")
          write_tsv(x = sigdf,
                    file = tempFileName,
                    col_names = FALSE)
          
        }
      }
      if (type == "biocrates") {
        cat("#### No Volcano Plots for Biocrates analysis since we are using difference of means and not linear or log2 fold change. #### \n")
      } else {
        # volcano plots
        # if compute_log2fc is F, volcano plots should not be generated and a message should be displayed
        source(paste0(scriptPath,"/plotVolcano.R"))
        if (compute_log2fc & mytest == "t-test") {
          plotVolcano(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".xlsx"),
                      myComparison = mycomparison,
                      outDir = paste0(outdir, "/volcano_plots/"),
                      fcCutoff = as.numeric(settings["linear_fc_cutoff",1]),
                      padjCutoff = as.numeric(settings["padj_cutoff",1]),
                      padjMethod = settings["padj_method",1])
        }
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
                  reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".xlsx"),
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
                  reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison,".xlsx"),
                  cutoffStat = settings["padj_method",1],
                  cutoff = 1,
                  samplesAreRows = T)
    }
    
    # Write one workbook with a sheet per comparison listing features dropped by the missingness filter
    if (length(dropped_reports) > 0) {
      names(dropped_reports) <- make.unique(substr(gsub("[^A-Za-z0-9_. -]", "_", names(dropped_reports)), 1, 31))
      write.xlsx(dropped_reports,
                 paste0(outdir, "/report/DroppedFeatures.xlsx"),
                 rowNames = F, overwrite = T)
      print(cat("##### Missingness filter: dropped-feature report written to report/DroppedFeatures.xlsx #####\n"))
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
  }
  
  # Check if "Rplots.pdf" exists in the current directory
  file_path <- file.path(getwd(), "Rplots.pdf")
  
  # If the file exists, delete it
  if (file.exists(file_path)) {
    file.remove(file_path)
  }
  
  # # Stop the parallel backend
  # stopCluster(cl)

  # End log
  sink(NULL)
}
