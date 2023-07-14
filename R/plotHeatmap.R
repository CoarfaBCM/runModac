plotHeatmap <- function(exprs,
                        meta,
                        test,
                        comparison,
                        samplesAreRows = F,
                        outdir,
                        groupOrder,
                        groupColors,
                        heatmapColorScale,
                        reportfile,
                        cutoffStat = "fdr",
                        cutoff = 0.25) {
  # Loading required packages and installing ones not present
  list.of.packages <- c("dplyr", "circlize", "RColorBrewer")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
  }
  
  # Ensuring data has rows of samples and columns of features
  if (!samplesAreRows) {
    exprs <- t(exprs)
  }
  
  meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  
  meta <- meta[order(match(meta$group, groupOrder)), , drop=F]
  exprs <- exprs[match(rownames(meta), rownames(exprs)),]
  
  # function for creating folders
  createDir <- function(folder) {
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
  }
  createDir(outdir)
  
  mygrouping <- factor(meta[,1], levels = groupOrder)
  exprdf <- data.frame(t(exprs), check.rows = F, check.names = F)
  
  if (grepl("fdr", cutoffStat, ignore.case = T)) {
    sigFeatures <- read.csv(reportfile) %>% filter(fdr <= cutoff) %>% select(1)
  } else if (grepl("pval", cutoffStat, ignore.case = T)) {
    sigFeatures <- read.csv(reportfile) %>% filter(pval <= cutoff) %>% select(1)
  }
  if (nrow(sigFeatures) == 0) {
    print(cat("No significant features at", cutoffStat, "<", cutoff,"\n"))
  } else {
    scaled_df <- as.data.frame(t(scale(t(exprdf[sigFeatures[,1],]))))
    scaled_df <- na.omit(scaled_df)
    
    tempdf <- cbind(c("", "Metabolite", rownames(scaled_df)),rbind(meta[,1], colnames(scaled_df), scaled_df))
    names(tempdf) <- NULL
    write.csv(tempdf, paste0(outdir,"/heatmap_",test,"_",comparison,"_",cutoffStat,cutoff,".csv"), row.names	= F)
    
    x<-abs(max(scaled_df))
    y<-abs(min(scaled_df))
    scale<-max(c(x,y))
    # if (heatmap_color_scale == 1) {
    #   heatmap_color_scale <- c("blue", "black", "yellow")
    # } else {
    #   heatmap_color_scale <- c("blue", "white", "red")
    # }
    col_heatmap <- colorRamp2(c((-scale)*.75, 0,(scale)*.75), heatmapColorScale)
    col_list <- list(group=groupColors)
    names(col_list$group) <- groupOrder
    
    ha1 <- HeatmapAnnotation(df = meta, show_annotation_name = T, annotation_height = .25,
                            col = col_list)
    
    ht<-Heatmap(matrix = as.matrix(scaled_df),
                top_annotation = ha1,
                clustering_distance_rows = "euclidean",
                name = "z-score",
                column_names_side = "bottom",
                col=col_heatmap,
                cluster_columns = F,
                cluster_rows = T,
                row_names_gp = gpar(fontsize=6),
                rect_gp = gpar(col = NA, lty = 1, lwd = 1),
                column_names_gp = gpar(fontsize=6),
                show_row_dend = T,
                km=1,
                heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(1, "in")))
    
    pdf(paste0(outdir,"/heatmap_",test,"_",comparison,"_",cutoffStat,cutoff,".pdf"), height = 14, width = 10)
    draw(ht, heatmap_legend_side="bottom")
    dev.off()
    
    jpeg(paste0(outdir,"/heatmap_",test,"_",comparison,"_",cutoffStat,cutoff,".jpg"), height = 14, width = 10, units = "in", quality = 1, res = 300)
    draw(ht, heatmap_legend_side="bottom")
    dev.off()
  }
}