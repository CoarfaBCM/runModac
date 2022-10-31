normISTD <- function(inputdf, qc.idx, comparisonsFile, myoutdir, mymethod) {
  istd.info <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", skip = 7), row.names = 1, check.rows = F) 
  
  # calculating CV
  istd.idx <- istd.info[mymethod, "ISTD_column"]
  
  if(ncol(inputdf) == istd.idx) {
    cv <- sd(inputdf[-c(qc.idx:nrow(inputdf)),istd.idx])/mean(inputdf[,istd.idx]) 
  } else {
    all.istd.idx <- istd.idx:ncol(inputdf)
    cv <- sapply(all.istd.idx, function(x){sd(inputdf[-c(qc.idx:nrow(inputdf)),x])/mean(inputdf[-c(qc.idx:nrow(inputdf)),x])})
    keep <- all.istd.idx[which(cv == min(cv))]
    inputdf <- inputdf[,c(1:(istd.idx-1),keep)]
    cv <- min(cv)
  }
  
  # Plotting ISTD
  temp <- inputdf[,ncol(inputdf)]
  
  pdf(paste(myoutdir,paste0("dIS_dist_",mymethod,".pdf"),sep = "/"))
  plot(log2(temp),
       ylim = c(ceiling(mean(log2(temp)))-4, ceiling(mean(log2(temp)))+4),
       xlim = c(0,nrow(inputdf)),
       xlab = "Samples",
       ylab = "Log2 ISTD",
       main = paste0(mymethod, " ISTD distribution (",colnames(inputdf)[ncol(inputdf)],") [CV = ",round(cv,2),"]"),
       cex.main = 0.8,
       col = "blue",
       pch = 19)
  dev.off()
  
  # ISTD normalization
  inputdf[1:(ncol(inputdf)-1)] <- apply(inputdf[1:(ncol(inputdf)-1)], 2, function(x){x/inputdf[,ncol(inputdf)]})
  
  # dropping istd
  inputdf <- inputdf[,-ncol(inputdf)]
  
  return(inputdf)
}