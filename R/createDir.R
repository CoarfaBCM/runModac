# function for creating folders
createDir <- function(folder) {
  if (!dir.exists(folder)){dir.create(folder, recursive = TRUE)}       
}