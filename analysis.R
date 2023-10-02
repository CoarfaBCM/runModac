source("~/Box/runModac/R/runModac.R")

# example1-metabolomics
runModac(inputFile = "~/Box/runModac/examples/input-metabolomics.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-metabolomics.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-metabolomics",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         # samplesAreRows = T,
         # sampleIDRow = 2,
         # min_signal = NULL,
         scriptPath = "~/Box/runModac/R/")

# example2-iqr
runModac(inputFile = "~/Box/runModac/examples/input-iqr.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-iqr.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-iqr",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [IQR]",
         heatmap_color_scale = c("blue","white","red"),
         # samplesAreRows = T,
         # sampleIDRow = 2,
         # min_signal = NULL,
         scriptPath = "~/Box/runModac/R/")

# example3-prenormalized
runModac(inputFile = "~/Box/runModac/examples/input-prenormalized.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-prenormalized.xlsx",
         type = "none",
         outdir = "~/Box/runModac/examples/results-prenormalized",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [prenormalized]",
         # heatmap_color_scale = c("blue","white","red"),
         # samplesAreRows = T,
         # sampleIDRow = 2,
         # min_signal = NULL,
         scriptPath = "~/Box/runModac/R/")

# example4-transposed
runModac(inputFile = "~/Box/runModac/examples/input-transposed.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-transposed.xlsx",
         type = "none",
         outdir = "~/Box/runModac/examples/results-transposed",
         project_title = "Project Report",
         project_subtitle = "Example w/ Transposed Data",
         # heatmap_color_scale = c("blue","white","red"),
         samplesAreRows = F,
         # sampleIDRow = 2,
         # min_signal = NULL,
         scriptPath = "~/Box/runModac/R/")

# example5-rppa
runModac(inputFile = "~/Box/runModac/examples/input-rppa.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-rppa.xlsx",
         type = "rppa",
         outdir = "~/Box/runModac/examples/results-rppa",
         project_title = "Project Report",
         project_subtitle = "RPPA Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         samplesAreRows = F,
         sampleIDRow = 1,
         # min_signal = NULL,
         scriptPath = "~/Box/runModac/R/")
