
runModac <- function(inputFile, comparisonsFile, type = "metabolomics", outdir = "results") {
  library(readxl)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  
  source("R/createDir.R")
  source("R/preprocessMetab.R")
  
  exprsdf <- preprocessMetab(inputFile = inputFile,
                             comparisonsFile = comparisonsFile,
                             outdir = outdir)
  
  # reading in all comparisons
  all.comparisons <- excel_sheets(comparisonsFile)[-1]
  
  # PCA plot
  ## read in groups
  mygroups <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = all.comparisons[1], skip = 3), row.names = 2, check.rows = F)
  
  source("R/plotPCA.R")
  plotPCA(exprs = exprsdf[["raw"]], meta = mygroups, outdir = paste0(outdir, "/pca"), suffix = "_raw", transpose = T)
  plotPCA(exprs = exprsdf[["norm"]], meta = mygroups, outdir = paste0(outdir, "/pca"), suffix = "_norm", transpose = T)
  
  source("R/myStatTest.R")
  for (i in seq_along(all.comparisons)) {
    mycomparison <- all.comparisons[i]
    settings <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 8, col_names = F), row.names = 1, check.rows = F)
    myttest(exprs = exprsdf[["norm"]],
            meta = mygroups,
            comparison = mycomparison,
            outdir = paste0(outdir, "/report/"))
  }
}

runModac(inputFile = "scratch/input/input-A.xlsx",
         comparisonsFile = "scratch/input/comparisons-A.xlsx",
         outdir = "scratch/results-A")
