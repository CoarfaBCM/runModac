source("R/runModac.R")

# scratch1 - metabolomics
runModac(inputFile = "example1-metabolomics/input/input-A.xlsx",
         comparisonsFile = "example1-metabolomics/input/comparisons-A.xlsx",
         outdir = "scratch/results-ex1",
         samplesAreRows = T,
         scriptPath = "~/Box Sync/runModac/R/",
         heatmap_color_scale = c("blue","white","red"))
# 
# # scratch2 - prenormalized
# runModac(inputFile = "example2-prenormalized-dummy/input/input.xlsx",
#          comparisonsFile = "example2-prenormalized-dummy/input/comparisons.xlsx",
#          type = "metabolomics",
#          outdir = "scratch/results-ex2",
#          scriptPath = "~/Box Sync/runModac/R/")
#
# scratch6 - rppa
runModac(inputFile = "example6-rppa/input/input.xlsx",
         comparisonsFile = "example6-rppa/input/comparisons.xlsx",
         type = "rppa",
         outdir = "scratch/results-ex6",
         samplesAreRows = F,
         sampleIDRow = 1,
         scriptPath = "~/Box Sync/runModac/R/")
# 
# # input4
# runModac(inputFile = "scratch/input4/input.xlsx",
#          comparisonsFile = "scratch/input4/comparisons.xlsx",
#          type = "metabolomics",
#          outdir = "scratch/results4",
#          scriptPath = "~/Box Sync/runModac/R/")
# 
# # input3
# runModac(inputFile = "scratch/input3/input.xlsx",
#          comparisonsFile = "scratch/input3/comparisons.xlsx",
#          type = "metabolomics",
#          outdir = "scratch/results3",
#          scriptPath = "~/Box Sync/runModac/R/")
# 
# # input2
# runModac(inputFile = "scratch/input2/input.xlsx",
#          comparisonsFile = "scratch/input2/comparisons.xlsx",
#          type = "metabolomics",
#          outdir = "scratch/results2",
#          scriptPath = "~/Box Sync/runModac/R/")
# 
# # input1
# runModac(inputFile = "scratch/input1/input-A.xlsx",
#          comparisonsFile = "scratch/input1/comparisons-A.xlsx",
#          type = "metabolomics",
#          outdir = "scratch/results1",
#          scriptPath = "~/Box Sync/runModac/R/")
# 