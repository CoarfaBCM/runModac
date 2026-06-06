normIQR <- function(inputdf, input_space, diff_space, qc_row = NA) {

  log2TF <- input_space == "linear" & diff_space == "log2"
  linearTF <- input_space == "log2" & diff_space == "linear"

  qc.idx <- as.numeric(qc_row) - 1

  # log2 transform
  if(log2TF){
    inputdf <- log2(inputdf)
  }

  #linear transform
  if(linearTF){
    inputdf <- 2^inputdf
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

  if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
    return_qc_num <- 0
  } else {
    return_qc_num <- length(qc.idx:nrow(inputdf))
  }

  return(list(data = newdf,
              qc_num = return_qc_num
              )
         )
}
