normIQR <- function(inputdf, comparisonsFile) {
  
  settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 7, col_names = F), row.names = 1, check.rows = F))
  log2tr <- as.logical(settings["log2transform",1])
  
  # log2 transform
  if(log2tr){
    inputdf <- log2(inputdf)
  }
  
  all.medians <- unname(apply(inputdf,1,median))
  all.upper.quartile <- unname(apply(inputdf,1,function(x){quantile(x, 0.75)}))
  all.lower.quartile <- unname(apply(inputdf,1,function(x){quantile(x, 0.25)}))
  newdf <- data.frame(t(sapply(1:nrow(inputdf),function(x){
    sapply(1:ncol(inputdf),function(y){
      inputdf[x,y] <- (inputdf[x,y] - all.medians[x])/(all.upper.quartile[x] - all.lower.quartile[x])
    })
  })), row.names = rownames(inputdf))
  colnames(newdf) <- colnames(inputdf)
  
  return(newdf)
}