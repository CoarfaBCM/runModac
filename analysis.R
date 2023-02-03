source("R/runModac.R")

# example1-metabolomics
runModac(inputFile = "example1-metabolomics/input/input-A.xlsx",
         comparisonsFile = "example1-metabolomics/input/comparisons-A.xlsx",
         type = "metabolomics",
         outdir = "example1-metabolomics/results-A",
         scriptPath = "~/Box Sync/runModac/R/")

# example2-prenormalized-dummy
runModac(inputFile = "example2-prenormalized-dummy/input/input.xlsx",
         comparisonsFile = "example2-prenormalized-dummy/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "example2-prenormalized-dummy/results",
         scriptPath = "~/Box Sync/runModac/R/")

# example3-iqr
runModac(inputFile = "example3-iqr/input/input.xlsx",
         comparisonsFile = "example3-iqr/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "example3-iqr/results",
         scriptPath = "~/Box Sync/runModac/R/")

# example4-prenormalized-nodummy
runModac(inputFile = "example4-prenormalized-nodummy/input/input.xlsx",
         comparisonsFile = "example4-prenormalized-nodummy/input/comparisons.xlsx",
         type = "metabolomics",
         outdir = "example4-prenormalized-nodummy/results",
         scriptPath = "~/Box Sync/runModac/R/")

# example5-transposed
runModac(inputFile = "example5-transposed/input/input.xlsx",
         comparisonsFile = "example5-transposed/input/comparisons.xlsx",
         type = "metabolomics",
         samplesAreRows = F,
         outdir = "example5-transposed/results",
         scriptPath = "~/Box Sync/runModac/R/")

# example6-rppa
runModac(inputFile = "example6-rppa/input/input.xlsx",
         comparisonsFile = "example6-rppa/input/comparisons.xlsx",
         type = "rppa",
         outdir = "example6-rppa/results",
         samplesAreRows = F,
         sampleIDRow = 1,
         scriptPath = "~/Box Sync/runModac/R/")
