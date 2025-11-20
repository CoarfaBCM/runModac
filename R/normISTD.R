normISTD <- function(inputdf, comparisonsFile, myoutdir, mymethod) {
  
  istd.info <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", skip = 9), check.rows = F))
  rownames(istd.info) <- istd.info$tab
  istd.info <- istd.info[,-1, drop=F]
  settings <- suppressMessages(data.frame(read_excel(comparisonsFile, trim_ws = T, sheet = "settings", n_max = 9, col_names = F), row.names = 1, check.rows = F))
  input_space <- settings["input_data_space",1]
  diff_space <- settings["differential_analysis_space",1]
  
  log2TF <- input_space == "linear" & diff_space == "log2"
  linearTF <- input_space == "log2" & diff_space == "linear"
  
  # saving qc row start idx
  qc.idx <- as.numeric(settings["QC_row",1])-1
  if (is.na(qc.idx) | qc.idx > nrow(inputdf)) {
    inputdf.noqc <- inputdf
  } else {
    inputdf.noqc <- inputdf[-c(qc.idx:nrow(inputdf)),]
  }
  
  # calculating CV
  istd.idx <- istd.info[mymethod, "ISTD_column"]-1
  
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
  if(cv > settings["cv_cutoff_internal_standard",1]) {
    stop(paste0("For ",mymethod,", ",colnames(inputdf)[ncol(inputdf)]," CV = ", round(cv,2), " > CV cuttoff = ", settings["cv_cutoff_internal_standard",1]))
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