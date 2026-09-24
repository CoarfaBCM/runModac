source("~/Box/runModac/R/runModac.R")

# ---------------------------------------------------------------------------
# comparisonsFile format (v2): 2-sheet Excel
#   Sheet 1 "samples"     — SampleID | GroupingVariable1 | GroupingVariable2 ...
#   Sheet 2 "comparisons" — ComparisonName | TestType | Test | Control | GroupingVariable | GroupOrder (optional) | Pairing (optional)
#
# All analysis settings are now passed directly as runModac() arguments.
# GroupOrder (optional column in "comparisons"): comma-separated group names
# controlling heatmap/boxplot display order only (Test/Control designation
# and PCA plots are unaffected). Falls back to the default derived order
# when absent or empty for a row.
# ---------------------------------------------------------------------------

# example1-metabolomics
runModac(inputFile    = "~/Box/runModac/examples/input-metabolomics.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-metabolomics-v2.xlsx",
         type         = "metabolomics",
         outdir       = "~/Box/runModac/examples/results-metabolomics_v2",
         project_title    = "Project Report",
         project_subtitle = "Metabolomics Analysis",
         scriptPath   = "~/Box/runModac/R/",
         # Data orientation
         samplesAreRows = TRUE,
         # Normalization settings
         normalization              = "istd",
         input_data_space           = "linear",
         differential_analysis_space = "log2",
         QC_row                     = 41,
         cv_cutoff_internal_standard = 0.25,
         ISTD_columns               = c(AminoAcids = 25),
         # Statistical settings
         padj_method     = "fdr",
         padj_cutoff     = 0.25,
         linear_fc_cutoff = 1,
         compute_log2fc  = TRUE)

# example2-iqr
runModac(inputFile = "~/Box/runModac/examples/input-iqr.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-iqr-v2.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-iqr_v2",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [IQR]",
         heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/",
         samplesAreRows = TRUE,
         normalization = "iqr",
         input_data_space = "linear",
         differential_analysis_space = "log2",
         QC_row = 41,
         padj_method = "fdr",
         padj_cutoff = 0.25,
         linear_fc_cutoff = 1,
         compute_log2fc = TRUE)

# example3-prenormalized
runModac(inputFile = "~/Box/runModac/examples/input-prenormalized.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-prenormalized-v2.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-prenormalized_v2",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [prenormalized]",
         scriptPath = "~/Box/runModac/R/",
         samplesAreRows = TRUE,
         normalization = "none",
         input_data_space = "log2",
         differential_analysis_space = "log2",
         padj_method = "fdr",
         padj_cutoff = 0.25,
         linear_fc_cutoff = 1,
         compute_log2fc = TRUE)

# example4-transposed
runModac(inputFile = "~/Box/runModac/examples/input-transposed.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-transposed-v2.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-transposed_v2",
         project_title = "Project Report",
         project_subtitle = "Example w/ Transposed Data",
         scriptPath = "~/Box/runModac/R/",
         samplesAreRows = FALSE,
         normalization = "none",
         input_data_space = "log2",
         differential_analysis_space = "log2",
         padj_method = "fdr",
         padj_cutoff = 0.25,
         linear_fc_cutoff = 1,
         compute_log2fc = TRUE)

# example5-rppa
runModac(inputFile = "~/Box/runModac/examples/input-rppa.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-rppa-v2.xlsx",
         type = "rppa",
         outdir = "~/Box/runModac/examples/results-rppa_v2",
         project_title = "Project Report",
         project_subtitle = "RPPA Analysis",
         scriptPath = "~/Box/runModac/R/",
         samplesAreRows = FALSE,
         sampleIDRow = 1,
         normalization = "none",
         input_data_space = "linear",
         differential_analysis_space = "log2",
         padj_method = "pval",
         padj_cutoff = 0.05,
         linear_fc_cutoff = 1,
         compute_log2fc = TRUE)
