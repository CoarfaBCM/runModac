source("~/Box/runModac/R/runModac.R")

# example1-metabolomics
runModac(inputFile = "~/Box/runModac/examples/input-metabolomics.xlsx", # path to input data Excel file
         comparisonsFile = "~/Box/runModac/examples/comparisons-metabolomics.xlsx", # path to Excel file defining comparisons, groups, and settings
         type = "metabolomics", # analysis type: "metabolomics", "rppa", "biocrates", or "none"
         outdir = "~/Box/runModac/examples/results-metabolomics", # directory where all outputs are written
         project_title = "Project Report", # title shown on the PPTX report
         project_subtitle = "Metabolomics Analysis", # subtitle shown on the PPTX report
         # heatmap_color_scale = c("blue", "black", "yellow"), # low/mid/high colors for heatmaps
         # samplesAreRows = T, # TRUE if input sheets have samples as rows (FALSE if transposed)
         # sampleIDRow = 2, # which metadata row/column holds sample IDs (used for RPPA aggregation)
         # max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
         # min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
         # replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
         # replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
         # onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
         scriptPath = "~/Box/runModac/R/") # directory containing the pipeline's R scripts

# example2-iqr
runModac(inputFile = "~/Box/runModac/examples/input-iqr.xlsx", # path to input data Excel file
         comparisonsFile = "~/Box/runModac/examples/comparisons-iqr.xlsx", # path to Excel file defining comparisons, groups, and settings
         type = "metabolomics", # analysis type: "metabolomics", "rppa", "biocrates", or "none"
         outdir = "~/Box/runModac/examples/results-iqr", # directory where all outputs are written
         project_title = "Project Report", # title shown on the PPTX report
         project_subtitle = "Metabolomics Analysis [IQR]", # subtitle shown on the PPTX report
         heatmap_color_scale = c("blue","white","red"), # low/mid/high colors for heatmaps
         # samplesAreRows = T, # TRUE if input sheets have samples as rows (FALSE if transposed)
         # sampleIDRow = 2, # which metadata row/column holds sample IDs (used for RPPA aggregation)
         # max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
         # min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
         # replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
         # replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
         # onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
         scriptPath = "~/Box/runModac/R/") # directory containing the pipeline's R scripts

# example3-prenormalized
runModac(inputFile = "~/Box/runModac/examples/input-prenormalized.xlsx", # path to input data Excel file
         comparisonsFile = "~/Box/runModac/examples/comparisons-prenormalized.xlsx", # path to Excel file defining comparisons, groups, and settings
         type = "none", # analysis type: "metabolomics", "rppa", "biocrates", or "none"
         outdir = "~/Box/runModac/examples/results-prenormalized", # directory where all outputs are written
         project_title = "Project Report", # title shown on the PPTX report
         project_subtitle = "Metabolomics Analysis [prenormalized]", # subtitle shown on the PPTX report
         # heatmap_color_scale = c("blue","white","red"), # low/mid/high colors for heatmaps
         # samplesAreRows = T, # TRUE if input sheets have samples as rows (FALSE if transposed)
         # sampleIDRow = 2, # which metadata row/column holds sample IDs (used for RPPA aggregation)
         # max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
         # min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
         # replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
         # replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
         # onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
         scriptPath = "~/Box/runModac/R/") # directory containing the pipeline's R scripts

# example4-transposed
runModac(inputFile = "~/Box/runModac/examples/input-transposed.xlsx", # path to input data Excel file
         comparisonsFile = "~/Box/runModac/examples/comparisons-transposed.xlsx", # path to Excel file defining comparisons, groups, and settings
         type = "none", # analysis type: "metabolomics", "rppa", "biocrates", or "none"
         outdir = "~/Box/runModac/examples/results-transposed", # directory where all outputs are written
         project_title = "Project Report", # title shown on the PPTX report
         project_subtitle = "Example w/ Transposed Data", # subtitle shown on the PPTX report
         # heatmap_color_scale = c("blue","white","red"), # low/mid/high colors for heatmaps
         samplesAreRows = F, # TRUE if input sheets have samples as rows (FALSE if transposed)
         # sampleIDRow = 2, # which metadata row/column holds sample IDs (used for RPPA aggregation)
         # max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
         # min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
         # replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
         # replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
         # onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
         scriptPath = "~/Box/runModac/R/") # directory containing the pipeline's R scripts

# example5-rppa
runModac(inputFile = "~/Box/runModac/examples/input-rppa.xlsx", # path to input data Excel file
         comparisonsFile = "~/Box/runModac/examples/comparisons-rppa.xlsx", # path to Excel file defining comparisons, groups, and settings
         type = "rppa", # analysis type: "metabolomics", "rppa", "biocrates", or "none"
         outdir = "~/Box/runModac/examples/results-rppa", # directory where all outputs are written
         project_title = "Project Report", # title shown on the PPTX report
         project_subtitle = "RPPA Analysis", # subtitle shown on the PPTX report
         # heatmap_color_scale = c("blue","white","red"), # low/mid/high colors for heatmaps
         samplesAreRows = F, # TRUE if input sheets have samples as rows (FALSE if transposed)
         sampleIDRow = 1, # which metadata row/column holds sample IDs (used for RPPA aggregation)
         # max_missing_frac = 1, # missingness filter: drop a feature if >(max_missing_frac*100)% missing in ANY comparison group (1 = off). A per-group "at least 2 non-missing" rule is always enforced.
         # min_signal = NULL, # drop a feature if no sample (raw) exceeds this signal cutoff (NULL = off)
         # replaceNA = NULL, # value to replace NAs with before normalization (NULL = leave as NA)
         # replaceZeros = NULL, # value to replace zeros with before normalization (NULL = leave as 0; set to 1 for non-biocrates metabolomics core requests as per core request)
         # onlyNorm = F, # if TRUE, stop after normalization/QC (skip stats, plots, reports)
         scriptPath = "~/Box/runModac/R/") # directory containing the pipeline's R scripts
