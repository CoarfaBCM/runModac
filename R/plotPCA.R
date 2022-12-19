plotPCA <- function(exprs, meta, outdir, suffix = NULL, label = T, labSize = 2, transpose = F) {
  library(ggfortify)
  library(openxlsx)
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
    meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  }
  
  if (transpose) {
    pcaRaw <- prcomp(exprs)
  } else {
      pcaRaw <- prcomp(t(exprs))
    }
  
  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
  }
  
  createDir(outdir)
  
  saveFile <- paste0("pca",suffix,".pdf")
  write.xlsx(pcaRaw, paste(outdir,gsub("pdf","xlsx",saveFile),sep="/"), overwrite = T)
  pdf(paste(outdir,saveFile,sep="/"))
  print(autoplot(pcaRaw, data = meta, colour = colnames(meta)[1]))
  if (label) {
    print(autoplot(pcaRaw,
                   data = meta,
                   colour = colnames(meta)[1],
                   label = T,
                   label.size = labSize))
  }
  dev.off()
}