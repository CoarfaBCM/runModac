normISTD <- function(inputdf, comparisonsFile, myoutdir, mymethod) {
  istd.info <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", skip = 8), check.rows = F)
  rownames(istd.info) <- istd.info$tab
  istd.info <- istd.info[,-1, drop=F]
  settings <- data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 8, col_names = F), row.names = 1, check.rows = F)
  log2tr <- as.logical(settings["log2transform",1])
  
  # saving qc row start idx
  qc.idx <- as.numeric(settings["QC_row",1])
  
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
  
  # checking if cv meets user specified cv cutoff
  if(cv > settings["cv_cutoff_internal_standard",1]) {
    stop(paste0("cv = ", cv, "< cv cuttoff = ", settings["cv_cutoff_internal_standard",1]))
  }
  
  # plotting ISTD
  temp <- inputdf[,ncol(inputdf)]
  
  if(log2tr){
    temp <- log2(temp)
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
  
  # ISTD normalization
  inputdf[1:(ncol(inputdf)-1)] <- apply(inputdf[1:(ncol(inputdf)-1)], 2, function(x){x/inputdf[,ncol(inputdf)]})
  
  # dropping istd
  inputdf <- inputdf[,-ncol(inputdf)]
  
  return(inputdf)
}