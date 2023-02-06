createRNK <- function(reportFile, outFile) {
  
  data <- read.csv(reportFile, check.names = F)
  data <- data[order(data[,grepl("fc",colnames(data))], decreasing = T),]
  
  write.table(data[,c(1, which(grepl("fc",colnames(data))))], outFile, quote = F, sep = "\t", row.names = F)
}