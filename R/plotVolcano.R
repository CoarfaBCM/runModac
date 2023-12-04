# Volcano Plot
plotVolcano <- function(reportFile, myComparison, outDir, fcCutoff = 1.5, padjMethod = "fdr", padjCutoff = 0.25, label = T, suffix = NULL) {
  
  list.of.packages <- c("ggplot2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
  }
  createDir(outDir)
  
  # data <- read.csv(reportFile, check.names = F)
  # data$DE <- "No"
  # data$DE[(data[,"log2_fc"] > log2(fcCutoff) | data[,"log2_fc"] < -log2(fcCutoff)) & data[,padjMethod] < padjCutoff] <- "DE"
  
  data <- read.csv(reportFile, check.names = F)
  
  data$direction <- "NS"
  data$direction[data[,"log2_fc"] > log2(fcCutoff) & data[,padjMethod] < padjCutoff] <- "up"
  data$direction[data[,"log2_fc"] < -log2(fcCutoff) & data[,padjMethod] < padjCutoff] <- "down"
  
  mycolors <- c("red", "grey30", "blue")
  names(mycolors) <- c("up","NS", "down")
  
  mylabels <- data$ID
  mylabels[data$direction == "NS"] <- ""
  
  # myvolcanoplot <- ggplot(data=data, aes(x=log2_fc, y=-log10(get(padjMethod)), col=DE)) + 
  #   geom_point(show.legend = T) + 
  #   theme_minimal() + 
  #   scale_color_manual(values=c("red", "black")) +
  #   geom_hline(yintercept = -log10(padjCutoff)) +
  #   labs(x = 'log2FC', y = paste0("-log10(",padjMethod,")")) +
  #   {if(fcCutoff != 1) geom_vline(xintercept = c(-log2(fcCutoff), log2(fcCutoff)))}
  
  myvolcanoplot1 <- ggplot2::ggplot(data=data, aes(x=log2_fc, y=-log10(get(padjMethod)), col=direction)) +
    ggplot2::geom_point(alpha = 0.5, size = 2) +
    ggplot2::scale_color_manual(values=mycolors) +
    ggplot2::scale_y_continuous(n.breaks = 10, labels = function(x) {
      scales::label_number(accuracy = 0.05)(x)
    }) +
    {if(any(-log10(data[[padjMethod]]) > -log10(padjCutoff))) {
      ggplot2::geom_hline(yintercept = -log10(padjCutoff), linetype = "dashed", colour = "goldenrod", size = 0.75, alpha = 0.5)
    }} +
    {if(any(data$log2_fc > log2(fcCutoff))) {
      ggplot2::geom_vline(xintercept = log2(fcCutoff), linetype = "dashed", colour = "red", size = 0.75, alpha = 0.5)
    }} +
    {if(any(data$log2_fc < -log2(fcCutoff))) {
      ggplot2::geom_vline(xintercept = -log2(fcCutoff), linetype = "dashed", colour = "blue", size = 0.75, alpha = 0.5)
    }} +
    ggplot2::labs(
      caption = stringr::str_interp('${sum(data$direction != "NS")}/${nrow(data)} signficant; ${sum(data$direction == "down")} down, ${sum(data$direction == "up")} up\n${padjMethod}: ${padjCutoff}, log2FC: ${round(log2(fcCutoff), 2)}, linearFC: ${round(fcCutoff, 2)}'),
      x = "log2FC",
      y = paste0("-log10(",padjMethod,")")
    ) +
    ggplot2::theme(legend.position = "bottom")
  
  if (label) {
    myvolcanoplot2 <- myvolcanoplot1 +
      ggrepel::geom_text_repel(ggplot2::aes(label = mylabels), size = 3, show.legend = FALSE)
  }
  
  pdf(paste0(outDir,"/VolcanoPlot_",myComparison,"_FC",fcCutoff,"_",padjMethod,padjCutoff,suffix,".pdf"),
      height = 7, width = 7)
  print(myvolcanoplot1)
  if (label) {print(myvolcanoplot2)}
  dev.off()
  
  jpeg(paste0(outDir,"/VolcanoPlot_",myComparison,"_FC",fcCutoff,"_",padjMethod,padjCutoff,suffix,".jpg"))
  print(myvolcanoplot1)
  dev.off()
  
  write.csv(myvolcanoplot1$data,
            paste0(outDir,"/VolcanoPlot_",myComparison,"_FC",fcCutoff,"_",padjMethod,padjCutoff,suffix,".csv"),
            quote = F,
            row.names = F)
}
