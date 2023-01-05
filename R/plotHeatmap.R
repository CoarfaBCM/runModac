library(dplyr)
library(circlize)
library(ComplexHeatmap)

runHeatmap <- function(exprs, meta, reportfile, outdir, cutoffStat, cutoff) {
  rawdf <- openxlsx::read.xlsx(inputfile, check.names = F)
  
  if (is.character(exprs)) {
    exprs <- read.csv(exprs, row.names = 1, check.names = F)
  }
  
  if (is.character(meta)) {
    meta <- read.csv(meta, row.names = 1, check.names = F)
    meta <- meta[match(rownames(exprs), rownames(meta)), , drop=F]
  }
  
  mygrouping <- as.factor(mymeta[,1])
  exprdf <- t(data.frame(rawdf[,-1], row.names = 1, check.rows = F, check.names = F))
  
  sigFeatures <- read.csv(reportfile) %>% filter(fdr < fdrCutoff) %>% select(1)
  
  scaled_df <- as.data.frame(t(scale(t(exprdf[sigFeatures[,1],]))))
  
  x<-abs(max(scaled_df))
  y<-abs(min(scaled_df))
  scale<-max(c(x,y))
  # scale <- 4
  col_heatmap <- colorRamp2(c((-scale)*.75, 0,(scale)*.75), c("BLUE2", "white", "RED"))
  
  ht <- Heatmap(matrix = scaled_df,
                col=col_heatmap,
                name = "z-score",
                rect_gp = gpar(col = NA, lty = 1, lwd = 1),
                cluster_rows = T,
                row_names_gp = gpar(fontsize=6),
                cluster_columns = F,
                column_split = mygrouping,
                column_title = NULL,
                top_annotation = HeatmapAnnotation(Group = anno_block(gp = gpar(fill=c("green", "orange", "pink","aquamarine","purple")),
                                                                      labels = levels(mygrouping),
                                                                      labels_gp = gpar(font=2, fontsize=12))),
                show_column_names = F,
                border = T,
                show_heatmap_legend = T,
                heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop"))
  
  pdf(paste0(outdir,"/heatmap-fdr",fdrCutoff,".pdf"), height = 14, width = 10)
  draw(ht, heatmap_legend_side="bottom")
  dev.off()
}