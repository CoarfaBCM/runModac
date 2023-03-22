createSignature <- function(reportFile, outFile, fcCutoff, statType, statCutoff) {
  # Loading required packages and installing ones not present
  list.of.packages <- c("tidyverse")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  data <- read.csv(reportFile, check.names = F)
  data <- data[order(data[,grepl("fc",colnames(data))], decreasing = T),]
  
  temp.idx.statType <- which(grepl(statType, colnames(data), ignore.case = T))
  temp.idx.fc <- which(grepl("fc", colnames(data), ignore.case = T))
  
  finaldf <- data %>% 
    filter(.[[temp.idx.statType]] <= as.numeric(statCutoff)) %>% 
    filter(abs(.[[temp.idx.fc]]) >= as.numeric(fcCutoff)) %>% 
    select(1, all_of(temp.idx.fc))
  write_tsv(x = finaldf, 
            file = outFile,
            col_names = FALSE)
}
