library(dplyr)
library(circlize)
library(ComplexHeatmap)

plotHeatmap <- function(exprs, meta, test, comparison, outdir, groupOrder, reportfile, cutoffStat = "fdr", cutoff = 0.25) {
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
    meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  }
  
  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
  }
  createDir(outdir)
  
  mygrouping <- factor(meta[,1], levels = groupOrder)
  exprdf <- data.frame(t(exprs), check.rows = F, check.names = F)
  
  if (cutoffStat == "fdr") {
    sigFeatures <- read.csv(reportfile) %>% filter(fdr < cutoff) %>% select(1)
  } else if (cutoffStat == "pval") {
    sigFeatures <- read.csv(reportfile) %>% filter(pval < cutoff) %>% select(1)
  }
  if (nrow(sigFeatures) == 0) {
    return(cat("No significant features at ", cutoffStat, " < ", cutoff))
  } else {
    scaled_df <- as.data.frame(t(scale(t(exprdf[sigFeatures[,1],]))))
    
    tempdf <- cbind(c("", "Metabolite", rownames(scaled_df)),rbind(meta[,1], colnames(scaled_df), scaled_df))
    names(tempdf) <- NULL
    write.csv(tempdf, paste0(outdir,"/heatmap_",test,"_",comparison,"_",cutoffStat,cutoff,".csv"), row.names	= F)
    browser()
    x<-abs(max(scaled_df))
    y<-abs(min(scaled_df))
    scale<-max(c(x,y))
    # scale <- 4
    col_heatmap <- colorRamp2(c((-scale)*.75, 0,(scale)*.75), c("BLUE2", "black", "yellow"))
    
    ht <- Heatmap(matrix = scaled_df,
                  col=col_heatmap,
                  name = "z-score",
                  rect_gp = gpar(col = NA, lty = 1, lwd = 1),
                  cluster_rows = T,
                  row_names_gp = gpar(fontsize=5),
                  cluster_columns = F,
                  column_split = mygrouping,
                  column_title = NULL,
                  top_annotation = HeatmapAnnotation(Group = anno_block(gp = gpar(fill=c("red","steelblue","green", "orange", "pink","aquamarine","purple","grey","black","maroon","khaki")),
                                                                        labels = levels(mygrouping),
                                                                        labels_gp = gpar(font=2, fontsize=7))),
                  show_column_names = F,
                  border = T,
                  show_heatmap_legend = T,
                  heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop"))
    
    pdf(paste0(outdir,"/heatmap_",test,"_",comparison,"_",cutoffStat,cutoff,".pdf"), height = 14, width = 10)
    draw(ht, heatmap_legend_side="bottom")
    dev.off()
  }
}