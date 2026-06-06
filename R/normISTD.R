normISTD <- function(inputdf, input_space, diff_space, qc_row = NA, cv_cutoff, istd_column, myoutdir, mymethod) {

  log2TF <- input_space == "linear" & diff_space == "log2"
  linearTF <- input_space == "log2" & diff_space == "linear"

  # saving qc row start idx
  qc.idx <- as.numeric(qc_row) - 1
  if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
    inputdf.noqc <- inputdf
  } else {
    inputdf.noqc <- inputdf[-c(qc.idx:nrow(inputdf)),]
  }

  # calculating CV
  istd.idx <- istd_column - 1

  if(ncol(inputdf) == istd.idx) {
    cv <- sd(inputdf.noqc[,istd.idx])/mean(inputdf.noqc[,istd.idx])
  } else {
    all.istd.idx <- istd.idx:ncol(inputdf)
    cv <- sapply(all.istd.idx, function(x){sd(inputdf.noqc[,x])/mean(inputdf.noqc[,x])})
    keep <- all.istd.idx[which(cv == min(cv))]
    inputdf <- inputdf[,c(1:(istd.idx-1),keep)]
    cv <- min(cv)
  }

  # checking if cv meets user specified cv cutoff
  if(cv > cv_cutoff) {
    stop(paste0("For ",mymethod,", ",colnames(inputdf)[ncol(inputdf)]," CV = ", round(cv,2), " > CV cutoff = ", cv_cutoff))
  }

  # plotting ISTD
  temp <- inputdf[,ncol(inputdf)]

  # log2 transform
  if(log2TF){
    temp <- log2(temp)
  }

  #linear transform
  if(linearTF){
    temp <- 2^temp
  }

  pdf(paste(myoutdir,paste0("dIS_dist_",mymethod,".pdf"),sep = "/"))
  plot(temp,
       ylim = c(ceiling(mean(temp))-4, ceiling(mean(temp))+4),
       xlim = c(0,nrow(inputdf)),
       xlab = "Samples",
       ylab = "Log2 ISTD",
       main = paste0(mymethod, " ISTD distribution (",colnames(inputdf)[ncol(inputdf)],") [CV = ",round(cv,2),"]"),
       cex.main = 0.8,
       col = "blue",
       pch = 19)
  dev.off()

  jpeg(paste(myoutdir,paste0("dIS_dist_",mymethod,".jpg"),sep = "/"))
  plot(temp,
       ylim = c(ceiling(mean(temp))-4, ceiling(mean(temp))+4),
       xlim = c(0,nrow(inputdf)),
       xlab = "Samples",
       ylab = "Log2 ISTD",
       main = paste0(mymethod, " ISTD distribution (",colnames(inputdf)[ncol(inputdf)],") [CV = ",round(cv,2),"]"),
       cex.main = 0.8,
       col = "blue",
       pch = 19)
  dev.off()

  # ISTD normalization
  if(input_space == "linear") {
    inputdf[1:(ncol(inputdf)-1)] <- apply(inputdf[1:(ncol(inputdf)-1)], 2, function(x){x/inputdf[,ncol(inputdf)]})
  } else if(input_space == "log2") {
    inputdf[1:(ncol(inputdf)-1)] <- apply(inputdf[1:(ncol(inputdf)-1)], 2, function(x){x - inputdf[,ncol(inputdf)]})
  }

  # log2 transform
  if(log2TF){
    inputdf <- log2(inputdf)
  }

  #linear transform
  if(linearTF){
    inputdf <- 2^inputdf
  }

  # dropping istd
  if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
    return_qc_num <- 0
  } else {
    return_qc_num <- length(qc.idx:nrow(inputdf))
  }
  return(list(data = inputdf[,-ncol(inputdf)],
              istd = colnames(inputdf)[ncol(inputdf)],
              cv = round(cv,2),
              qc_num = return_qc_num
              )
         )
}
