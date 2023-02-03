
aggregateRPPA <- function(inputFile, outFile, geneNames = F, sampleIDRow = 2) {
  
  library(openxlsx)
  mydf <- openxlsx::read.xlsx(inputFile)
  mygroups <- data.frame(colnames(mydf)[-c(1:5)], t(mydf[1,-c(1:5)]))
  names(mygroups) <- c("ID", "Sample")
  if (geneNames) {
    mydf <- data.frame(mydf[-1,-c(1:3,5)], row.names = NULL, check.rows = F, check.names = F)
    mydf[,1] <- unname(sapply(mydf[,1], function(x){strsplit(x, ",")[[1]][1]}))
  } else {
    mydf <- data.frame(mydf[-1,-c(1,3:5)], row.names = NULL, check.rows = F, check.names = F)
    mydf[,1] <- gsub("_R_V","",mydf[,1])
  }
  
  mydf[,-1] <- sapply(mydf[,-1], as.numeric)
  mydf[is.na(mydf)] <- 1
  
  # keeping only antibodies or genes with the highest expression out of duplicates
  tempdf <- data.frame(table(mydf[,1]))
  dups <- as.character(tempdf$Var1[tempdf$Freq > 1])
  
  for (x in dups) {
    idx <- which(mydf[,1] == x)
    tempsums <- rowSums(mydf[idx, -1])
    tempmax <- max(tempsums)
    dropidx <- idx[which(tempsums != tempmax)]
    mydf <- mydf[-dropidx,]
  }
  
  mydf <- data.frame(Sample = mygroups$Sample, t(data.frame(mydf, row.names = 1, check.rows = F, check.names = F)), check.names = F, check.rows = F)
  newdf <- aggregate(. ~ Sample, data = mydf, FUN = median)
  
  if (sampleIDRow == 1) {
    uniqueNames <- data.frame(ID = unname(sapply(mygroups$ID, function(x){strsplit(x,"[.]")[[1]][1]})),
                              sample = mygroups$Sample)
    uniqueNames <- uniqueNames[!duplicated(uniqueNames$sample),]
    uniqueNames <- uniqueNames[match(newdf$Sample, uniqueNames$sample),]
    if (identical(newdf$Sample, uniqueNames$sample)) {
      newdf$Sample <- uniqueNames$ID 
    }
  }
  
  openxlsx::write.xlsx(list(rppa = newdf), outFile, rowNames = F) 
}
