source("R/runModac.R")

# example1-metabolomics
runModac(inputFile = "example1-metabolomics/input/input-A.xlsx",
         comparisonsFile = "example1-metabolomics/input/comparisons-A.xlsx",
         type = "metabolomics",
         outdir = "example1-metabolomics/results-A",
         scriptPath = "~/Box Sync/runModac/R/")

# example2-prenormalized
runModac(inputFile = "example2-prenormalized/input/input.xlsx",
         comparisonsFile = "example2-prenormalized/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "example2-prenormalized/results",
         scriptPath = "~/Box Sync/runModac/R/")
