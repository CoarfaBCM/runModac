
rppaAggr <- function(inputFile,
                     outputFile,
                     cv_cutoff = NULL,
                     sampleIDRow,
                     replace.cvcutoff = NULL,
                     replace.na = NULL) {
  # Loading required packages and installing ones not present
  list.of.packages <- c("readxl","openxlsx","dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)}
  lapply(list.of.packages, require, character.only = TRUE)

  tempfunc <- function(inputFile,
                       wb,
                       sheetName,
                       cv_cutoff,
                       sampleIDRow,
                       replace.cvcutoff = NULL,
                       replace.na = NULL) {
    mydf <- read.xlsx(inputFile, sheet = sheetName)
    
    if (!("AB_species" %in% names(mydf))) {
      rawdf <- mydf
      
      # trimming any trailing whitespaces from antibody IDs, antibody names and gene symbols
      rawdf[,c(1:5)] <- apply(rawdf[,c(1:5)], 2, trimws)
      
      # removing trailing characters other than antibody name
      rawdf$AB_name <- gsub("_V$|_$","", rawdf$AB_name)
      
      # adding antibody species column
      AB_species <- rep("",nrow(rawdf))
      AB_species[1] <- "AB_species"
      AB_species[grepl("_R$", rawdf$AB_name)] <- "rabbit"
      AB_species[grepl("_M$", rawdf$AB_name)] <- "mouse"
      
      rawdf$AB_name <- gsub("_R$","",rawdf$AB_name)
      rawdf$AB_name <- gsub("_M$","",rawdf$AB_name)
      # rawdf$Gene_ID <- unname(sapply(rawdf$Gene_ID, function(x){strsplit(x, ",")[[1]][1]}))
      
      # trimming any trailing whitespaces from antibody IDs, antibody names and gene symbols
      rawdf[,c(1:5)] <- apply(rawdf[,c(1:5)], 2, trimws)
      
      colnames(rawdf)[1] <- "AB_ID"
      rawdf$AB_species <- AB_species
      mycols1 <- c("AB_ID","AB_name","AB_species","Slide_file","Gene_ID", "Swiss_ID")
      mycols2 <- setdiff(names(rawdf), mycols1)
      rawdf <- rawdf[, c(mycols1, mycols2)]
      writeData(wb,sheetName,rawdf,rowNames = FALSE)
      # bold the first column
      addStyle(wb = wb,
               sheet = sheetName,
               style = createStyle(textDecoration = "bold"),
               cols = 1,
               rows = 1:nrow(rawdf),
               gridExpand = T)
      
      # bold the first row
      addStyle(wb = wb,
               sheet = sheetName,
               style = createStyle(textDecoration = "bold"),
               cols = 1:ncol(rawdf),
               rows = 1,
               gridExpand = T)
      
      mydf <- rawdf 
    }
    
    mygroups <- data.frame(ID=unname(sapply(colnames(mydf)[-c(1:6)], function(x){strsplit(x,"[.]")[[1]][1]})),
                           Sample=t(mydf[1,-c(1:6)]))
    colnames(mygroups) <- c("ID","Sample")
    
    AB_species <- mydf$AB_species[-1]
    mydf <- data.frame(mydf[-1,-c(3,4,6)], row.names = NULL, check.rows = F, check.names = F)
    
    # cleaning up gene IDs
    mydf$Gene_ID <- unname(sapply(mydf$Gene_ID, function(x){strsplit(x, ",")[[1]][1]}))
    
    # trimming any trailing whitespaces from antibody IDs, antibody names and gene symbols
    mydf[,c(1:3)] <- apply(mydf[,c(1:3)], 2, trimws)
    
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
    
    mydf$AB_name <- makeUniqueNames(IDs = mydf$AB_ID, names = mydf$AB_name)
    
    mydf[,-c(1:3)] <- sapply(mydf[,-c(1:3)], as.numeric)
    
    if (!(is.null(replace.na))) {
      mydf[,-c(1:3)][is.na(mydf[,-c(1:3)])] <- replace.na
    }
    
    tempdf <- data.frame(Sample = mygroups[,sampleIDRow], t(data.frame(mydf[,-c(1,3)], row.names = 1, check.rows = F, check.names = F)), check.names = F, check.rows = F)
    
    # custom median and cv functions because the default R
    # functions drop a sample if it has all NAs for even one feature
    custom_median <- function(x) {
      if(all(is.na(x))) {
        NA
      } else {
        median(x, na.rm = T)
      }
    }
    
    custom_cv <- function(x) {
      if(all(is.na(x))) {
        NA
      } else {
        sd(x, na.rm = T)/mean(x, na.rm = T)
      }
    }
    
    newdf <- tempdf %>%
      group_by(Sample) %>%
      summarise(across(.cols = everything(), .fns = custom_median))
    cvdf <- tempdf %>%
      group_by(Sample) %>%
      summarise(across(.cols = everything(), .fns = custom_cv))
    
    # The commented out lines below drop a sample if all it's replicates have NAs for even just 1 feature
    # newdf <- aggregate(. ~ Sample, data = tempdf, FUN = function(x){median(x, na.rm = T)})
    # cvdf <- aggregate(. ~ Sample, data = tempdf, FUN = function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})
    
    newdf <- data.frame(newdf, row.names = 1, check.rows = F, check.names = F)
    newdf <- newdf[unique(mygroups[,sampleIDRow]),]
    
    cvdf <- data.frame(cvdf, row.names = 1, check.rows = F, check.names = F)
    cvdf <- cvdf[unique(mygroups[,sampleIDRow]),]
    
    if (!(is.null(replace.cvcutoff)) & !(is.null(cv_cutoff))) {
      newdf[cvdf > cv_cutoff] <- replace.cvcutoff
    }
    
    finaldf <- data.frame(AB_ID = as.numeric(mydf[,1]), AB_name = colnames(newdf), AB_species = AB_species, GeneSymbol = mydf[,3], t(newdf), check.rows = F, check.names = F)
    finaldf.cv <- data.frame(AB_ID = as.numeric(mydf[,1]), AB_name = colnames(cvdf), AB_species = AB_species, GeneSymbol = mydf[,3], t(cvdf), check.rows = F, check.names = F)
    # write.xlsx(list(rppa = finaldf), paste0(outdir, "/full_aggregate_data.xlsx"), rowNames = F)
    
    addWorksheet(wb,paste0(sheetName,"_Median"))
    writeData(wb,paste0(sheetName,"_Median"),finaldf,rowNames = FALSE)
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_Median"),
             style = createStyle(numFmt = "0.00"),
             cols = 5:ncol(finaldf.cv),
             rows = 2:(nrow(finaldf.cv)+1),
             gridExpand = T)
    
    # bold the first column
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_Median"),
             style = createStyle(textDecoration = "bold"),
             cols = 1,
             rows = 1:(nrow(finaldf.cv)+1),
             gridExpand = T)
    
    # bold the first row
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_Median"),
             style = createStyle(textDecoration = "bold"),
             cols = 1:ncol(finaldf.cv),
             rows = 1,
             gridExpand = T)
    
    addWorksheet(wb,paste0(sheetName,"_CV"))
    writeData(wb,paste0(sheetName,"_CV"),finaldf.cv,rowNames = FALSE)
    mystyle <- createStyle(fontColour = "#000000", bgFill = "#FFFF00")
    conditionalFormatting(wb,
                          paste0(sheetName,"_CV"),
                          cols = 5:ncol(finaldf.cv),
                          rows = 2:(nrow(finaldf.cv)+1),
                          rule = paste0(">",cv_cutoff),
                          style = mystyle)
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_CV"),
             style = createStyle(numFmt = "PERCENTAGE"),
             cols = 5:ncol(finaldf.cv),
             rows = 2:(nrow(finaldf.cv)+1),
             gridExpand = T)
    
    # bold the first column
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_CV"),
             style = createStyle(textDecoration = "bold"),
             cols = 1,
             rows = 1:(nrow(finaldf.cv)+1),
             gridExpand = T)
    
    # bold the first row
    addStyle(wb = wb,
             sheet = paste0(sheetName,"_CV"),
             style = createStyle(textDecoration = "bold"),
             cols = 1:ncol(finaldf.cv),
             rows = 1,
             gridExpand = T)
    
    return(wb)
  }
  
  all.sheets <- excel_sheets(inputFile)
  wb <- loadWorkbook(inputFile)
  if(!("Norm_Median" %in% all.sheets)) {
    wb1 <- tempfunc(inputFile = inputFile,
                    wb = wb,
                    sheetName = "Norm",
                    cv_cutoff = cv_cutoff,
                    sampleIDRow = sampleIDRow,
                    replace.cvcutoff = replace.cvcutoff,
                    replace.na = replace.na)
    saveWorkbook(wb1,outputFile,overwrite = TRUE)
  } else {
    saveWorkbook(wb,outputFile,overwrite = TRUE)
  }
}
