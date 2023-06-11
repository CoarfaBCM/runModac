source("~/Box/runModac/R/runModac.R")

# example1-metabolomics
runModac(inputFile = "~/Box/runModac/example1-metabolomics/input/input-A.xlsx",
         comparisonsFile = "~/Box/runModac/example1-metabolomics/input/comparisons-A.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/example1-metabolomics/results-A",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example2-prenormalized-dummy
runModac(inputFile = "~/Box/runModac/example2-prenormalized-dummy/input/input.xlsx",
         comparisonsFile = "~/Box/runModac/example2-prenormalized-dummy/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/example2-prenormalized-dummy/results",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [prenormalized w/ dummy QC]",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example3-iqr
runModac(inputFile = "~/Box/runModac/example3-iqr/input/input.xlsx",
         comparisonsFile = "~/Box/runModac/example3-iqr/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/example3-iqr/results",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [IQR]",
         heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example4-prenormalized-nodummy
runModac(inputFile = "~/Box/runModac/example4-prenormalized-nodummy/input/input.xlsx",
         comparisonsFile = "~/Box/runModac/example4-prenormalized-nodummy/input/comparisons.xlsx",
         type = "none",
         outdir = "~/Box/runModac/example4-prenormalized-nodummy/results",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [prenormalized]",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example5-transposed
runModac(inputFile = "~/Box/runModac/example5-transposed/input/input.xlsx",
         comparisonsFile = "~/Box/runModac/example5-transposed/input/comparisons.xlsx",
         type = "metabolomics",
         samplesAreRows = F,
         outdir = "~/Box/runModac/example5-transposed/results",
         project_title = "Project Report",
         project_subtitle = "Example w/ Transposed Data",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example6-rppa
runModac(inputFile = "~/Box/runModac/example6-rppa/input/input.xlsx",
         comparisonsFile = "~/Box/runModac/example6-rppa/input/comparisons.xlsx",
         type = "rppa",
         outdir = "~/Box/runModac/example6-rppa/results",
         samplesAreRows = F,
         sampleIDRow = 1,
         project_title = "Project Report",
         project_subtitle = "RPPA Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")
