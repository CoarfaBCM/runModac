source("~/Box/runModac/R/runModac.R")

# example1-metabolomics
runModac(inputFile = "~/Box/runModac/examples/input-ex1.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-ex1.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-ex1",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example2-iqr
runModac(inputFile = "~/Box/runModac/examples/input-ex2.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-ex2.xlsx",
         type = "metabolomics",
         outdir = "~/Box/runModac/examples/results-ex2",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [IQR]",
         heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example3-prenormalized
runModac(inputFile = "~/Box/runModac/examples/input-ex3.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-ex3.xlsx",
         type = "none",
         outdir = "~/Box/runModac/examples/results-ex3",
         project_title = "Project Report",
         project_subtitle = "Metabolomics Analysis [prenormalized]",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example4-transposed
runModac(inputFile = "~/Box/runModac/examples/input-ex4.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-ex4.xlsx",
         type = "none",
         samplesAreRows = F,
         outdir = "~/Box/runModac/examples/results-ex4",
         project_title = "Project Report",
         project_subtitle = "Example w/ Transposed Data",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")

# example5-rppa
runModac(inputFile = "~/Box/runModac/examples/input-ex5.xlsx",
         comparisonsFile = "~/Box/runModac/examples/comparisons-ex5.xlsx",
         type = "rppa",
         outdir = "~/Box/runModac/examples/results-ex5",
         samplesAreRows = F,
         sampleIDRow = 1,
         project_title = "Project Report",
         project_subtitle = "RPPA Analysis",
         # heatmap_color_scale = c("blue","white","red"),
         scriptPath = "~/Box/runModac/R/")
