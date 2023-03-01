
rppaAggr <- function(inputFile, outdir, geneNames = F, sampleIDRow = 2) {
  library(openxlsx)
  mydf <- openxlsx::read.xlsx(inputFile)
  mygroups <- data.frame(colnames(mydf)[-c(1:5)], t(mydf[1,-c(1:5)]))
  names(mygroups) <- c("ID", "Sample")
  
  mydf <- data.frame(mydf[-1,-c(1,3,5)], row.names = NULL, check.rows = F, check.names = F)
  mydf[,1] <- gsub("_R_V","",mydf[,1])
  mydf[,2] <- unname(sapply(mydf[,2], function(x){strsplit(x, ",")[[1]][1]}))
  
  mydf[,-c(1,2)] <- sapply(mydf[,-c(1,2)], as.numeric)
  mydf[is.na(mydf)] <- 1
  
  tempdf <- data.frame(Sample = mygroups$Sample, t(data.frame(mydf[,-2], row.names = 1, check.rows = F, check.names = F)), check.names = F, check.rows = F)
  newdf <- aggregate(. ~ Sample, data = tempdf, FUN = median)
  cvdf <- aggregate(. ~ Sample, data = tempdf, FUN = function(x){sd(x)/mean(x)})
  newdf.1 <- newdf[,-1]
  cvdf.1 <- cvdf[,-1]
  newdf.1[cvdf.1 > 0.25] <- 1
  newdf[,-1] <- newdf.1
  
  if (sampleIDRow == 1) {
    uniqueNames <- data.frame(ID = unname(sapply(mygroups$ID, function(x){strsplit(x,"[.]")[[1]][1]})),
                              sample = mygroups$Sample)
    uniqueNames <- uniqueNames[!duplicated(uniqueNames$sample),]
    uniqueNames <- uniqueNames[match(newdf$Sample, uniqueNames$sample),]
    if (identical(newdf$Sample, uniqueNames$sample)) {
      newdf$Sample <- uniqueNames$ID 
    }
  }
  openxlsx::write.xlsx(list(rppa = newdf), paste0(outdir, "/aggregate_data.xlsx"), rowNames = F) 
  
  newdf <- data.frame(newdf, row.names = 1, check.rows = F, check.names = F)
  
  finaldf <- data.frame(GeneSymbol = mydf[,2], AB_name = colnames(newdf), t(newdf), check.rows = F, check.names = F)
  openxlsx::write.xlsx(list(rppa = finaldf), paste0(outdir, "/full_aggregate_data.xlsx"), rowNames = F)
}
