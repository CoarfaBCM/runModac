# runModac — Multiomics differential analysis pipeline
#
# PARAMETER REFERENCE
# ===================
#
# inputFile               Path to the input data Excel file (.xlsx).
#                         Each sheet = one analytical method (e.g. "AminoAcids", "Lipids").
#                         Rows = samples (if samplesAreRows = TRUE), columns = features.
#
# comparisonsFile         Path to the comparisons Excel file (.xlsx) in v2 format:
#                           Sheet "samples"     — SampleID | GroupingVar1 | GroupingVar2 | ...
#                           Sheet "comparisons" — ComparisonName | TestType | Test | Control | GroupingVariable | Pairing (optional)
#                         Pairing column (optional):
#                           Empty / NA  — unpaired analysis (default)
#                           Column name — name of a column in the "samples" sheet whose
#                                         values identify matched pairs (e.g. "SubjectID").
#                                         Each pair ID must appear exactly once in the Test
#                                         group and once in the Control group.
#                                         Supported by: t-test, limma, limma-full.
#                                         limma / limma-full use duplicateCorrelation + block.
#                         TestType values:
#                           "t-test" — Welch two-sample t-test; Test and Control must be filled
#                           "anova"  — one-way ANOVA across all levels of GroupingVariable;
#                                      Test and Control columns are ignored
#                           "limma"  — limma moderated t-test (lmFit / makeContrasts / eBayes)
#                                      fitted on Test + Control samples only;
#                                      Test and Control must be filled;     
#                                      differential_analysis_space must be "log2" — the
#                                      comparison is skipped with a warning otherwise
#                           "limma-full" — same as "limma" but lmFit is run on ALL samples
#                                      in the GroupingVariable (all group levels); the
#                                      Test vs Control contrast is then extracted via
#                                      makeContrasts; eBayes borrows variance information
#                                      across all groups rather than just the two being
#                                      compared; Test and Control must be filled;
#                                      differential_analysis_space must be "log2"
#
# type                    Data type. Controls RPPA aggregation and biocrates-specific FC logic.
#                         Values: "metabolomics" | "rppa" | "biocrates" | "none"
#
# outdir                  Path to the output directory. Created if it does not exist.
#                         Default: "results"
#
# project_title           Title string for the PowerPoint report cover slide.
#                         Default: "Report"
#
# project_subtitle        Subtitle string for the PowerPoint report cover slide.
#                         Default: "by"
#
# heatmap_color_scale     Character vector of 3 colors defining the heatmap gradient
#                         (low → mid → high z-score).
#                         Default: c("blue", "black", "yellow")
#                         Example: c("blue", "white", "red")
#
# scriptPath              Path to the directory containing all runModac R scripts.
#                         Must end without a trailing slash, e.g. "~/Box/runModac/R"
#
# samplesAreRows          Logical. TRUE if rows = samples in the input file; FALSE if
#                         columns = samples (transposed layout).
#                         Default: TRUE
#
# sampleIDRow             Integer. Row number containing sample IDs. Used only when
#                         type = "rppa" (passed to rppaAggr).
#                         Default: 2
#
# min_signal              Numeric or NULL. Features with max signal below this threshold
#                         are dropped before analysis. NULL = no filtering.
#                         Default: NULL
#
# replaceNA               Numeric or NULL. Value to substitute for NA entries in the
#                         input matrix before normalization.
#                         NULL = keep NAs as-is.
#                         Default: NULL
#
# replaceZeros            Numeric or NULL. Value to substitute for zero entries in the
#                         input matrix before normalization.
#                         NULL = keep zeros as-is.
#                         Common choices: 1 (linear), 0.001 (near-zero imputation)
#                         Default: NULL
#
# onlyNorm                Logical. If TRUE, run preprocessing/normalization only and
#                         skip all downstream statistical analysis and plotting.
#                         Default: FALSE
#
# --- Normalization & data space settings ---
#
# normalization           Normalization method to apply.
#                         Values: "none" | "istd" | "iqr"
#                           "none" — no normalization; apply space transform only
#                           "istd" — divide each feature by an internal standard column
#                                    (requires ISTD_columns)
#                           "iqr"  — IQR normalization: (value - median) / (Q3 - Q1)
#                                    per sample
#                         Default: "none"
#
# input_data_space        Space the raw input data is stored in.
#                         Values: "linear" | "log2"
#                         Default: "log2"
#
# differential_analysis_space
#                         Space in which differential analysis is performed.
#                         Values: "linear" | "log2"
#                         If input_data_space = "linear" and differential_analysis_space = "log2",
#                         data is log2-transformed before analysis.
#                         If input_data_space = "log2"  and differential_analysis_space = "linear",
#                         data is back-transformed (2^x) before analysis.
#                         Default: "log2"
#
# QC_row                  Integer or NA. 1-based row number (in the input sheet) where
#                         quality-control samples begin. Rows from QC_row onward are
#                         excluded from differential analysis but included in normalization.
#                         NA = no QC rows.
#                         Default: NA
#
# cv_cutoff_internal_standard
#                         Numeric. Maximum allowable coefficient of variation (CV) for the
#                         internal standard across samples. Analysis stops with an error if
#                         exceeded. Only used when normalization = "istd".
#                         Default: 0.25  (i.e. 25%)
#
# ISTD_columns            Named integer vector mapping each method (sheet name) to the
#                         1-based column number of its internal standard in that sheet.
#                         Required when normalization = "istd"; ignored otherwise.
#                         Example: c(AminoAcids = 25, Lipids = 18)
#                         Default: NULL
#
# --- Statistical analysis settings ---
#
# padj_method             Method for p-value adjustment (multiple testing correction).
#                         Values: any method accepted by p.adjust():
#                           "fdr" / "BH" — Benjamini-Hochberg FDR (recommended)
#                           "bonferroni"  — conservative Bonferroni
#                           "holm", "hochberg", "hommel", "BY", "none"
#                         Default: "fdr"
#
# padj_cutoff             Numeric threshold applied to adjusted p-values for significance
#                         calls in heatmaps, volcano plots, and signature files.
#                         Default: 0.25
#
# linear_fc_cutoff        Numeric. Fold-change threshold in linear scale used to define
#                         up/down-regulated features in volcano plots and signature files.
#                         Internally converted to log2 when compute_log2fc = TRUE.
#                         Example: 1.5 means |log2FC| > log2(1.5) ≈ 0.58
#                         Default: 1.5
#
# compute_log2fc          Logical. If TRUE, compute and report log2 fold-change in addition
#                         to linear fold-change for t-test comparisons. Set to FALSE for
#                         biocrates data (difference of means is used instead).
#                         Default: TRUE

runModac <- function(inputFile,
                     comparisonsFile,
                     type,
                     outdir = "results",
                     project_title = "Report",
                     project_subtitle = "by",
                     heatmap_color_scale = c("blue", "black", "yellow"),
                     scriptPath,
                     samplesAreRows = TRUE,
                     sampleIDRow = 2,
                     min_signal = NULL,
                     replaceNA = NULL,
                     replaceZeros = NULL,
                     onlyNorm = FALSE,
                     # --- Settings (previously in comparisonsFile settings sheet) ---
                     normalization = "none",
                     input_data_space = "log2",
                     differential_analysis_space = "log2",
                     QC_row = NA,
                     cv_cutoff_internal_standard = 0.25,
                     ISTD_columns = NULL,       # named vector, e.g. c(AminoAcids = 25)
                     padj_method = "fdr",
                     padj_cutoff = 0.25,
                     linear_fc_cutoff = 1.5,
                     compute_log2fc = TRUE) {

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

  # Set the number of cores to use
  num_cores <- 4

  # Register the parallel backend
  cl <- makeCluster(num_cores, outfile="")
  registerDoParallel(cl)

  # Convert the loop to parallel using foreach
  foreach(i = 1) %dopar% {

    # setting CRAN mirror
    options(repos = c(CRAN = "https://cran.r-project.org"))

    # installing required packages
    list.of.packages <- c("ggplot2", "readxl", "openxlsx","tidyr","foreach","doParallel")
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
                            outdir = outdir,
                            scriptPath = scriptPath,
                            normalization = normalization,
                            input_data_space = input_data_space,
                            differential_analysis_space = differential_analysis_space,
                            QC_row = QC_row,
                            cv_cutoff_internal_standard = cv_cutoff_internal_standard,
                            ISTD_columns = ISTD_columns,
                            samplesAreRows = FALSE,
                            replaceNA = NULL,
                            replaceZeros = NULL)
    } else {
      exprsdf <- preProcess(inputFile = inputFile,
                            outdir = outdir,
                            scriptPath = scriptPath,
                            normalization = normalization,
                            input_data_space = input_data_space,
                            differential_analysis_space = differential_analysis_space,
                            QC_row = QC_row,
                            cv_cutoff_internal_standard = cv_cutoff_internal_standard,
                            ISTD_columns = ISTD_columns,
                            samplesAreRows = samplesAreRows,
                            replaceNA = NULL,
                            replaceZeros = NULL)
    }
    print(cat("##### Preprocessing complete #####\n\n"))

    # Boxplot of all samples — skip when column is all NA
    tryCatch({
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
    }, error = function(e) {
      cat("##### Warning: could not generate all-samples boxplot:", conditionMessage(e), "#####\n")
    })

    # stop here if only normalization requested
    if (onlyNorm) {
      print("##### Normalization complete. No further downstream analysis requested. #####")
    } else {

      # ---- Read new-format comparisons file (2 sheets) ----
      comparisons_df <- suppressMessages(data.frame(
        read_excel(comparisonsFile, sheet = "comparisons", trim_ws = TRUE),
        check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE
      ))

      samples_df <- suppressMessages(data.frame(
        read_excel(comparisonsFile, sheet = "samples", trim_ws = TRUE),
        check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE
      ))
      rownames(samples_df) <- as.character(samples_df[, 1])
      samples_df <- samples_df[, -1, drop = FALSE]

      all.comparison.labels <- c()

      # Collect all unique group levels for color palette
      all.groups <- c()
      for (i in seq_len(nrow(comparisons_df))) {
        mygroupvar_i <- as.character(comparisons_df$GroupingVariable[i])
        if (tolower(as.character(comparisons_df$TestType[i])) %in% c("t-test", "limma")) {
          all.groups <- c(all.groups,
                          as.character(comparisons_df$Control[i]),
                          as.character(comparisons_df$Test[i]))
        } else {
          grp_vals <- samples_df[, mygroupvar_i]
          all.groups <- c(all.groups, unique(grp_vals[!is.na(grp_vals)]))
        }
      }
      all.groups <- unique(all.groups)

      # Expanding color set to color 36 groups
      colors_set <- c("#FD001C", "#16FF00", "#0000FF", "#63422A", "#FE16D4", "#00F7FF",
                      "#F1E500", "#D7B0FD", "#FD8E3D", "#008E3B", "#CD1660", "#000099",
                      "#CB0DF7", "#225C75", "#FC9DBD", "#22AAFE", "#920D83", "#FCCB95",
                      "#512EA3", "#ACF25F", "#00FFC4", "#94352A", "#858B22", "#A2D1F5",
                      "#FC83FD", "#946A8E", "#22A593", "#F2532A", "#FDBB0D", "#D63B9A",
                      "#ABEFAE", "#B37AFF", "#FCD5F2", "#005332", "#1CE26D", "#FC7C7D")
      col_list <- list(group = colors_set[seq_along(all.groups)])
      names(col_list$group) <- all.groups

      source(paste0(scriptPath,"/myStatTest.R"))
      source(paste0(scriptPath,"/plotPCA.R"))


###Begin loop through comparisons -------------------------------------------
      for (i in seq_len(nrow(comparisons_df))) {
        mycomparison <- as.character(comparisons_df$ComparisonName[i])
        mytest       <- as.character(comparisons_df$TestType[i])
        mygroupvar   <- as.character(comparisons_df$GroupingVariable[i])

        if (tolower(mytest) %in% c("t-test", "limma")) {
          mytest_group <- as.character(comparisons_df$Test[i])
          myctrl_group <- as.character(comparisons_df$Control[i])
          mygroups <- c(myctrl_group, mytest_group)  # control first, then test
        } else if (tolower(mytest) == "limma-full") {
          mytest_group <- as.character(comparisons_df$Test[i])
          myctrl_group <- as.character(comparisons_df$Control[i])
          # Include all group levels so lmFit sees the full sample set
          grp_vals <- samples_df[, mygroupvar]
          mygroups <- unique(grp_vals[!is.na(grp_vals)])
        } else {
          # ANOVA: use all levels in order of appearance in samples sheet
          grp_vals <- samples_df[, mygroupvar]
          mygroups <- unique(grp_vals[!is.na(grp_vals)])
          mytest_group <- NULL
          myctrl_group <- NULL
        }

        print(cat("##### Performing", mytest, "for comparison:", mycomparison, "#####\n"))

        # Build metadata from samples sheet
        mymeta <- data.frame(group = samples_df[, mygroupvar],
                             row.names = rownames(samples_df),
                             stringsAsFactors = FALSE)
        mymeta <- mymeta[!is.na(mymeta$group) & mymeta$group %in% mygroups, , drop = FALSE]

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

        myexprs.raw  <- exprsdf[["raw"]][commonsamples, , drop=F]
        myexprs.norm <- exprsdf[["norm"]][commonsamples, , drop=F]

        dropZeroVarArrays <- function(inputdf, checkRows = F, checkCols = T){
          if (checkCols) {
            flag <- colSums(inputdf, na.rm = T) == 0 | colSums(!is.na(inputdf)) == 0
            if (any(flag, na.rm = T)) {
              print(cat("##### Dropping features with no expression or missing values across all samples for this comparison #####\n"))
              print(cat("##### Number of features dropped:", length(flag), "#####\n"))
              print(cat("##### Exact features dropped:", colnames(inputdf)[!flag], "#####\n"))
              inputdf <- inputdf[,!flag, drop=F]
            }
          }

          if (checkRows) {
            flag <- rowSums(inputdf, na.rm = T) == 0 | rowSums(!is.na(inputdf)) == 0
            if (any(flag, na.rm = T)) {
              print(cat("##### Dropping samples with no expression or missing values across all features for this comparison #####\n"))
              print(cat("##### Number of samples dropped:", length(flag), "#####\n"))
              print(cat("##### Exact samples dropped:", rownames(inputdf)[!flag], "#####\n"))
              inputdf <- inputdf[!flag, , drop=F]
            }
          }

          return(inputdf)
        }

        myexprs.raw  <- dropZeroVarArrays(myexprs.raw,  checkRows = T, checkCols = F)
        myexprs.norm <- dropZeroVarArrays(myexprs.norm, checkRows = T, checkCols = F)

        # Extract pair IDs for paired testing (NULL when Pairing column is absent or empty)
        mypair_col <- NULL
        if ("Pairing" %in% colnames(comparisons_df)) {
          pval_raw <- as.character(comparisons_df$Pairing[i])
          if (!is.na(pval_raw) && nzchar(trimws(pval_raw))) {
            mypair_col <- trimws(pval_raw)
          }
        }
        mypair_id <- NULL
        if (!is.null(mypair_col)) {
          if (!mypair_col %in% colnames(samples_df)) {
            warning(paste0("Pairing column '", mypair_col, "' not found in samples sheet ",
                           "for comparison '", mycomparison, "'. Running unpaired."))
          } else {
            mypair_id <- samples_df[rownames(myexprs.norm), mypair_col]
            names(mypair_id) <- rownames(myexprs.norm)
          }
        }

        # Build unique output label: ComparisonName_TestType[_paired_by_PairingVar]
        mycomparison_full <- paste0(mycomparison, "_", mytest)
        if (!is.null(mypair_col)) {
          mycomparison_full <- paste0(mycomparison_full, "_paired_by_", mypair_col)
        }

        # PCA plot
        if (ncol(myexprs.raw) > 1) {
          plotPCA(exprs = myexprs.raw,
                  meta = mymeta,
                  outdir = paste0(outdir, "/pca"),
                  suffix = paste0("_raw_",mycomparison_full),
                  groupColors = col_list$group[mygroups],
                  samplesAreRows = T)
        }

        if (ncol(myexprs.norm) > 1) {
          plotPCA(exprs = myexprs.norm,
                  meta = mymeta,
                  outdir = paste0(outdir, "/pca"),
                  suffix = paste0("_norm_",mycomparison_full),
                  groupColors = col_list$group[mygroups],
                  samplesAreRows = T)
        }

        # limma requires log2 space; skip comparison with a warning if space is wrong
        if (tolower(mytest) %in% c("limma", "limma-full") && differential_analysis_space != "log2") {
          warning(paste0("Comparison '", mycomparison, "' skipped: ", mytest, " requires ",
                         "differential_analysis_space = 'log2' but current space is '",
                         differential_analysis_space, "'. Set differential_analysis_space = 'log2' to run limma."))
          next
        }

        all.comparison.labels <- c(all.comparison.labels, mycomparison_full)

        myStatTest(exprs = myexprs.norm,
                   meta = mymeta,
                   test = mytest,
                   request.type = type,
                   comparison = mycomparison_full,
                   group.ctrl.test = mygroups,
                   group.colors = col_list$group[mygroups],
                   diff.space = differential_analysis_space,
                   compute.log2fc = compute_log2fc,
                   samplesAreRows = T,
                   test.group = mytest_group,
                   ctrl.group = myctrl_group,
                   pair_id = mypair_id,
                   outdir = paste0(outdir, "/report/"))

        if (type == "rppa") {
          fullreport <- read.xlsx(paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison_full,".xlsx"), colNames = T, rowNames = F)
          tempdf <- read.xlsx(paste0(outdir,"/report/full_aggregate_data.xlsx"), sheet = "Norm_Median", rowNames = F)
          geneSymbols <- unlist(unname(sapply(fullreport$ID[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]})))
          finalreport <- data.frame(GeneSymbol = c(NA,geneSymbols), fullreport)
          names(finalreport)[2] <- "AB_name"
          write.xlsx(finalreport,
                     paste0(outdir,"/report/FullReport_",mytest,"_",mycomparison_full,".xlsx"),
                     rowNames = F, overwrite = T)
        }

        if (tolower(mytest) %in% c("t-test", "limma", "limma-full")) {
          createDir(paste0(outdir,"/Rnk"))
          createDir(paste0(outdir,"/Signature"))

          source(paste0(scriptPath,"/createRNK.R"))
          createRNK(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison_full,".xlsx"),
                    outFile = paste0(outdir,"/Rnk/Rnk_",mycomparison_full,".rnk"))

          if (type == "biocrates") {
            cat("#### No Signatures for Biocrates analysis since we are using difference of means and not linear or log2 fold change. #### \n")
          } else {
            if (compute_log2fc | tolower(mytest) %in% c("limma", "limma-full")) {
              temp_fc_cutoff <- log2(as.numeric(linear_fc_cutoff))
              outFileName <- paste0(mycomparison_full,"_log2FC",round(temp_fc_cutoff,2),"_",padj_method,padj_cutoff)
            } else {
              temp_fc_cutoff <- as.numeric(linear_fc_cutoff)
              outFileName <- paste0(mycomparison_full,"_linearFC",round(temp_fc_cutoff,2),"_",padj_method,padj_cutoff)
            }

            source(paste0(scriptPath,"/createSignature.R"))
            createSignature(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison_full,".xlsx"),
                            outFileName = outFileName,
                            outDir = paste0(outdir,"/Signature/"),
                            fcCutoff = temp_fc_cutoff,
                            statType = padj_method,
                            statCutoff = as.numeric(padj_cutoff))
          }

          if (type == "rppa") {
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

            sigdf$X1[-1] <- unname(sapply(sigdf$X1[-1], function(x){tempdf$GeneSymbol[tempdf$AB_name == x]}))
            sigdf <- sigdf[!is.na(sigdf[,1]), ,drop=F]
            tempFileName <- paste0(outdir,"/Signature/sig.",outFileName,".txt")
            write_tsv(x = sigdf,
                      file = tempFileName,
                      col_names = FALSE)
          }
        }

        if (type == "biocrates") {
          cat("#### No Volcano Plots for Biocrates analysis since we are using difference of means and not linear or log2 fold change. #### \n")
        } else {
          source(paste0(scriptPath,"/plotVolcano.R"))
          if (tolower(mytest) %in% c("limma", "limma-full") | (compute_log2fc & tolower(mytest) == "t-test")) {
            plotVolcano(reportFile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison_full,".xlsx"),
                        myComparison = mycomparison_full,
                        outDir = paste0(outdir, "/volcano_plots/"),
                        fcCutoff = as.numeric(linear_fc_cutoff),
                        padjCutoff = as.numeric(padj_cutoff),
                        padjMethod = padj_method)
          }
        }

        # heatmaps
        source(paste0(scriptPath,"/plotHeatmap.R"))
        plotHeatmap(exprs = myexprs.norm,
                    meta = mymeta,
                    test = mytest,
                    comparison = mycomparison_full,
                    outdir = paste0(outdir, "/heatmaps/"),
                    groupOrder = mygroups,
                    groupColors = col_list$group[mygroups],
                    heatmapColorScale = heatmap_color_scale,
                    reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison_full,".xlsx"),
                    cutoffStat = padj_method,
                    cutoff = as.numeric(padj_cutoff),
                    samplesAreRows = T)
        plotHeatmap(exprs = myexprs.norm,
                    meta = mymeta,
                    test = mytest,
                    comparison = mycomparison_full,
                    outdir = paste0(outdir, "/heatmaps/"),
                    groupOrder = mygroups,
                    groupColors = col_list$group[mygroups],
                    heatmapColorScale = heatmap_color_scale,
                    reportfile = paste0(outdir,"/report/Report_",mytest,"_",mycomparison_full,".xlsx"),
                    cutoffStat = padj_method,
                    cutoff = 1,
                    samplesAreRows = T)
      }

      source(paste0(scriptPath,"/createPptx.R"))
      template_pptx_path <- paste0(gsub("/R$|/R/$","/data",scriptPath),"/Project_Report_Template.pptx")
      createPptx(project_title = project_title,
                 project_subtitle = project_subtitle,
                 exprsdf = exprsdf,
                 normalization = normalization,
                 input_data_space = input_data_space,
                 differential_analysis_space = differential_analysis_space,
                 padj_method = padj_method,
                 padj_cutoff = padj_cutoff,
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

  # Stop the parallel backend
  stopCluster(cl)

  # End log
  sink(NULL)
}
