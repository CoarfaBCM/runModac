createRNK <- function(reportFile, outFile) {
  data <- read.xlsx(reportFile, check.names = F)
  if (any(grepl("log2",colnames(data), ignore.case = T))) {
    data <- data[order(data[,which(grepl("log2",colnames(data)))], decreasing = T),]
    write.table(data[,c(1, which(grepl("log2",colnames(data))))], outFile, quote = F, sep = "\t", row.names = F)
  } else {
    data <- data[order(data[,which(grepl("fc",colnames(data)))], decreasing = T),]
    write.table(data[,c(1, which(grepl("fc",colnames(data))))], outFile, quote = F, sep = "\t", row.names = F)
  }
}