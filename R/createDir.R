# function for creating folders
createDir <- function(folder) {
  if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
}