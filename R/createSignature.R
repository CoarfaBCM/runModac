createSignature <- function(reportFile, outFileName, outDir, fcCutoff, statType, statCutoff) {
  # Loading required packages and installing ones not present
  list.of.packages <- c("tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  data <- read.csv(reportFile, check.names = F)
  
  temp.idx.statType <- which(grepl(statType, colnames(data), ignore.case = T))
  if (any(grepl("log2",colnames(data), ignore.case = T))){
    temp.idx.fc <- which(grepl("log2", colnames(data), ignore.case = T))
  } else {
    temp.idx.fc <- which(grepl("fc", colnames(data), ignore.case = T))
  }
  data <- data[order(data[,temp.idx.fc], decreasing = T),]
  
  finaldf <- data %>% 
    filter(.[[temp.idx.statType]] <= as.numeric(statCutoff)) %>% 
    filter(abs(.[[temp.idx.fc]]) >= as.numeric(fcCutoff)) %>% 
    select(1, all_of(temp.idx.fc))
  
  # tempname <- basename(outFile) %>% sub("[.]txt$","",.) %>% sub("sig.","",.,ignore.case = T)
  # finaldf <- rbind(c(tempname,""),finaldf)
  # write_tsv(x = finaldf, 
  #           file = outFile,
  #           col_names = FALSE)
  
  finaldf <- rbind(c(outFileName,""),finaldf)
  write_tsv(x = finaldf,
            file = paste0(outDir,"/Sig_",outFileName,".txt"),
            col_names = FALSE,
            na = "")
}
