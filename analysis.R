source("/Users/amritkoirala/Library/CloudStorage/Box-Box/Amrit_Projects_Box/GitHub_codeNnotes/runModac/R/runModac.R")

# ---------------------------------------------------------------------------
# New comparisonsFile format (v2): 2-sheet Excel
#   Sheet 1 "samples"     — SampleID | GroupingVariable1 | GroupingVariable2 ...
#   Sheet 2 "comparisons" — ComparisonName | TestType | Test | Control | GroupingVariable
#
# All analysis settings are now passed directly as runModac() arguments.
# ---------------------------------------------------------------------------

# example1-metabolomics (new comparisons format)
runModac(inputFile    = "./examples/input-metabolomics.xlsx",
         comparisonsFile = "./examples/comparisons-metabolomics-v2.xlsx",
         type         = "metabolomics",
         outdir       = "./examples/results-metabolomics_v2",
         project_title    = "Project Report",
         project_subtitle = "Metabolomics Analysis",
         scriptPath   = "/Users/amritkoirala/Library/CloudStorage/Box-Box/Amrit_Projects_Box/GitHub_codeNnotes/runModac/R/",
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

# ---------------------------------------------------------------------------
# Examples 2-5 below use comparisonsFile in the OLD format (pre-v2).
# They need their comparisons files converted to the new 2-sheet format
# (sheet "samples" + sheet "comparisons") before they will run.
# Settings that were in the comparisonsFile settings sheet must now be
# supplied as arguments to runModac() — see example1 above for reference.
# ---------------------------------------------------------------------------

# example2-iqr (comparisons file needs updating to v2 format)
# runModac(inputFile = "~/Box/runModac/examples/input-iqr.xlsx",
#          comparisonsFile = "~/Box/runModac/examples/comparisons-iqr-v2.xlsx",
#          type = "metabolomics",
#          outdir = "~/Box/runModac/examples/results-iqr",
#          project_title = "Project Report",
#          project_subtitle = "Metabolomics Analysis [IQR]",
#          heatmap_color_scale = c("blue","white","red"),
#          scriptPath = "~/Box/runModac/R/",
#          normalization = "iqr",
#          input_data_space = "linear",
#          differential_analysis_space = "log2",
#          padj_method = "fdr",
#          padj_cutoff = 0.25,
#          linear_fc_cutoff = 1.5,
#          compute_log2fc = TRUE)

# example3-prenormalized (comparisons file needs updating to v2 format)
# runModac(inputFile = "~/Box/runModac/examples/input-prenormalized.xlsx",
#          comparisonsFile = "~/Box/runModac/examples/comparisons-prenormalized-v2.xlsx",
#          type = "none",
#          outdir = "~/Box/runModac/examples/results-prenormalized",
#          project_title = "Project Report",
#          project_subtitle = "Metabolomics Analysis [prenormalized]",
#          scriptPath = "~/Box/runModac/R/",
#          normalization = "none",
#          input_data_space = "log2",
#          differential_analysis_space = "log2",
#          padj_method = "fdr",
#          padj_cutoff = 0.25,
#          linear_fc_cutoff = 1.5,
#          compute_log2fc = TRUE)

# example4-transposed (comparisons file needs updating to v2 format)
# runModac(inputFile = "~/Box/runModac/examples/input-transposed.xlsx",
#          comparisonsFile = "~/Box/runModac/examples/comparisons-transposed-v2.xlsx",
#          type = "none",
#          outdir = "~/Box/runModac/examples/results-transposed",
#          project_title = "Project Report",
#          project_subtitle = "Example w/ Transposed Data",
#          scriptPath = "~/Box/runModac/R/",
#          samplesAreRows = FALSE,
#          normalization = "none",
#          input_data_space = "log2",
#          differential_analysis_space = "log2",
#          padj_method = "fdr",
#          padj_cutoff = 0.25,
#          linear_fc_cutoff = 1.5,
#          compute_log2fc = TRUE)

# example5-rppa (comparisons file needs updating to v2 format)
# runModac(inputFile = "~/Box/runModac/examples/input-rppa.xlsx",
#          comparisonsFile = "~/Box/runModac/examples/comparisons-rppa-v2.xlsx",
#          type = "rppa",
#          outdir = "~/Box/runModac/examples/results-rppa",
#          project_title = "Project Report",
#          project_subtitle = "RPPA Analysis",
#          scriptPath = "~/Box/runModac/R/",
#          samplesAreRows = FALSE,
#          sampleIDRow = 1,
#          normalization = "none",
#          input_data_space = "log2",
#          differential_analysis_space = "log2",
#          padj_method = "fdr",
#          padj_cutoff = 0.25,
#          linear_fc_cutoff = 1.5,
#          compute_log2fc = TRUE)
