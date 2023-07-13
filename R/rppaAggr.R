
rppaAggr <- function(inputFile,
                     outdir,
                     cv_cutoff = 0.25,
                     sampleIDRow = 2,
                     replacement = 1) {
  # Loading required packages and installing ones not present
  list.of.packages <- c("openxlsx")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  mydf <- read.xlsx(inputFile)
  
  mygroups <- data.frame(colnames(mydf)[-c(1:5)], t(mydf[1,-c(1:5)]))
  names(mygroups) <- c("ID", "Sample")
  
  tempIDs <- mydf$X1[-1]
  mydf <- data.frame(mydf[-1,-c(1,3,5)], row.names = NULL, check.rows = F, check.names = F)
  mydf[,1] <- gsub("_R_V","",mydf[,1])
  mydf[,2] <- unname(sapply(mydf[,2], function(x){strsplit(x, ",")[[1]][1]}))
  
  # Function to add ID to duplicate names
  makeUniqueNames <- function(IDs, names) {
    unique_names <- character(length(names))  # Initialize a vector to store unique names
    
    # Store the indices of duplicates
    duplicate_indices <- which(duplicated(names))
    
    # Iterate over each name
    for (i in seq_along(names)) {
      # Check if the name is a duplicate
      if (i %in% duplicate_indices) {
        # Find the index of the first instance of the current name
        first_instance_index <- min(which(names == names[i]))
        
        # Append the specific ID to both the first instance and the duplicates
        unique_names[first_instance_index] <- paste(names[first_instance_index], IDs[first_instance_index], sep = "_")
        unique_names[i] <- paste(names[i], IDs[i], sep = "_")
      } else {
        unique_names[i] <- names[i]  # Name is unique, no need to modify it
      }
    }
    
    return(unique_names)
  }
  
  mydf$AB_name <- makeUniqueNames(IDs = tempIDs, names = mydf$AB_name)
  
  mydf[,-c(1,2)] <- sapply(mydf[,-c(1,2)], as.numeric)
  mydf[is.na(mydf)] <- 1
  
  tempdf <- data.frame(Sample = mygroups$Sample, t(data.frame(mydf[,-2], row.names = 1, check.rows = F, check.names = F)), check.names = F, check.rows = F)
  newdf <- aggregate(. ~ Sample, data = tempdf, FUN = median)
  cvdf <- aggregate(. ~ Sample, data = tempdf, FUN = function(x){sd(x)/mean(x)})
  cvdf.1 <- cvdf[,-1]
  if (!(is.null(replacement))) {
    newdf.1 <- newdf[,-1]
    newdf.1[cvdf.1 > cv_cutoff] <- replacement
    newdf[,-1] <- newdf.1
  }
  
  if (sampleIDRow == 1) {
    uniqueNames <- data.frame(ID = unname(sapply(mygroups$ID, function(x){strsplit(x,"[.]")[[1]][1]})),
                              sample = mygroups$Sample)
    uniqueNames <- uniqueNames[!duplicated(uniqueNames$sample),]
    uniqueNames <- uniqueNames[match(newdf$Sample, uniqueNames$sample),]
    if (identical(newdf$Sample, uniqueNames$sample)) {
      newdf$Sample <- uniqueNames$ID 
    }
  }
  
  write.xlsx(list(rppa = newdf), paste0(outdir, "/aggregate_data.xlsx"), rowNames = F)
  newdf <- data.frame(newdf, row.names = 1, check.rows = F, check.names = F)
  # rownames(cvdf.1) <- rownames(newdf)
  
  finaldf <- data.frame(GeneSymbol = mydf[,2], AB_name = colnames(newdf), t(newdf), check.rows = F, check.names = F)
  # finaldf.cv <- data.frame(GeneSymbol = mydf[,2], AB_name = colnames(cvdf.1), t(cvdf.1), check.rows = F, check.names = F)
  write.xlsx(list(rppa = finaldf), paste0(outdir, "/full_aggregate_data.xlsx"), rowNames = F)
}
