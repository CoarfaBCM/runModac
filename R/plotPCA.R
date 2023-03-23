plotPCA <- function(exprs, meta, outdir, suffix = NULL, label = T, labSize = 2, samplesAreRows = F) {
  
  # Loading required packages and installing ones not present
  list.of.packages <- c("ggfortify", "openxlsx")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
  }
  
  if (!samplesAreRows) {
    exprs <- t(exprs)
  }
  
  keep <- intersect(rownames(exprs), rownames(meta))
  exprs <- exprs[keep,,drop=F]
  meta <- meta[keep,,drop=F]
  if (!(identical(rownames(exprs), rownames(meta)))) {
    stop(cat("##### Error matching samples between expression and meta data. #####"))
  }
  # meta <- meta[match(colnames(exprs), rownames(meta)), , drop=F]
  
  pcaRaw <- prcomp(exprs)
  
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
  
  saveFile <- paste0("pca",suffix,".jpg")
  jpeg(paste(outdir,saveFile,sep="/"), quality = 1)
  print(autoplot(pcaRaw, data = meta, colour = colnames(meta)[1]))
  dev.off()
}