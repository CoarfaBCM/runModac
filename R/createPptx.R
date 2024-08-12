createPptx <- function(project_title = "Project Report",
                       project_subtitle = "Metabolomics Analysis",
                       exprsdf = exprsdf,
                       settings = settings,
                       all.comparison.labels = all.comparison.labels,
                       outdir = outdir,
                       template_pptx_path) {
  
  # Loading required packages and installing ones not present
  list.of.packages <- c("officer", "tidyverse", "jpeg")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)}
  lapply(list.of.packages, require, character.only = TRUE)
  
  normalization <- settings["normalization",1]
  input_space <- settings["input_data_space",1]
  diff_space <- settings["differential_analysis_space",1]
  log2TF <- input_space == "linear" & diff_space == "log2"
  linearTF <- input_space == "log2" & diff_space == "linear"
  
  sample_number<-nrow(exprsdf[["norm"]]) #Not including the quality control samples
  qc_sample_number<-exprsdf$qc_num #"Quality Control" samples
  number_of_species<-ncol(exprsdf[["norm"]]) #metabolites or lipids, not including the internal standards
  number_of_methods<-length(exprsdf$methods)
  
  if (normalization == "istd") {
    norm_method<-sapply(1:number_of_methods,
                        function(x){
                          paste0(exprsdf$methods[x],
                                 " was normalized by Internal Standard ",
                                 exprsdf$istd[x],
                                 " [CV = ",
                                 exprsdf$cv[x],
                                 "]")
                        }) 
  } else if (normalization == "iqr") {
    norm_method<-sapply(1:number_of_methods,
                        function(x){
                          paste0(exprsdf$methods[x],
                                 " was normalized by IQR")
                        })
  } else if (normalization == "none") {
    if (log2TF) {
      norm_method<-c("Data was prenormalized and transformed to log2 space")
    } else if (linearTF) {
      norm_method<-c("Data was prenormalized and transformed to linear space") 
    } else {
      norm_method<-c("Data was prenormalized") 
    }
  }
  
  if (all(grepl("_over_", all.comparison.labels, ignore.case = T))) {
    comparison <- c("t-test (2 tailed, assuming unequal variance)",
                    all.comparison.labels)
    comparison_levels <- c(2, rep(3, length(all.comparison.labels)))
  } else {
    myflag <- grepl("_over_", all.comparison.labels, ignore.case = T)
    comparison <- c("t-test (2 tailed, assuming unequal variance)",
                    all.comparison.labels[myflag],
                    "anova (analysis of variance)",
                    all.comparison.labels[!(myflag)])
    comparison_levels <- c(2, rep(3, length(all.comparison.labels[myflag])),
                           2, rep(3, length(all.comparison.labels[!(myflag)])))
  }
  
  #define where the template powerpoint is located
  report <- read_pptx(template_pptx_path)
  
  # Functions ---------------------------------------------------------------
  ctrl_ppt<-function(report,files,title,names) {
    for(i in seq_along(files)){
      img<-readJPEG(files[i])
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      
      report<-report%>%add_slide(layout="Title and Content", master="Custom Design")
      report<-report%>%ph_with(value=str_c(title,names[i]), location=ph_location_type("title"))
      report<-report%>%ph_with(value=external_img(files[i]),ph_location(top=1,left=0.5,width=width,height=height))
    }
  }
  
  # Report Preparation ------------------------------------------------------
  # Title Slide -------------------------------------------------------------
  report<-report%>%add_slide(layout="Title Slide", master="Custom Design")
  report<-report%>%ph_with(location=ph_location_type("ctrTitle"),value=project_title)
  report<-report%>%ph_with(location=ph_location_type("subTitle"),value=project_subtitle)
  # Project Overview --------------------------------------------------------
  report<-report%>%add_slide(layout="Title and Content", master="Custom Design")
  report<-report%>%ph_with(location=ph_location_type("title"),value="Project Overview")
  if (qc_sample_number == 0) {
    ul_levels <- c(1,2,
                   1,2,
                   1,comparison_levels)
    ul_str <- c("Samples",
                paste0(sample_number," experimental samples"),
                paste0("Total of ",number_of_species," features from ",number_of_methods," methods"),
                norm_method,
                "Comparisons performed",
                comparison)
  } else {
    ul_levels <- c(1,2,2,
                   1,rep(2,number_of_methods),
                   1,comparison_levels)
    ul_str <- c("Samples",
                paste0(sample_number," experimental samples"),
                paste0(qc_sample_number," quality control samples"),
                paste0("Total of ",number_of_species," features from ",number_of_methods," methods"),
                norm_method,
                "Comparisons performed",
                comparison)
  }
  ul <- unordered_list(
    level_list = ul_levels,
    str_list = ul_str
  )
  report <- report%>%ph_with(location = ph_location_type(type = "body"), value = ul)
  
  # QA Plots ----------------------------------------------------------------
  qa<-list.files(path = paste0(outdir,"/QA_plots/"),
                 pattern="*.jpg")
  
  qa_plots<-qa[grepl("boxplot_",qa,ignore.case = T)]
  qa_names<-str_replace(qa_plots,"boxplot_","")%>%str_replace(".jpg","")
  qc_plots<-qa[grepl("qc_",qa,ignore.case = T)]
  qc_names<-str_replace(qc_plots,"qc_dist_","")%>%str_replace(".jpg","")
  IS_plots<-qa[grepl("dIS",qa)]
  IS_names<-str_replace(IS_plots,"dIS_dist_","")%>%str_replace(".jpg","")
  
  
  if(!is_empty(qa_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",qa_plots),title=c("QA plots: "),names=qa_names)}
  if(!is_empty(qc_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",qc_plots),title=c("Quality Control: "),names=qc_names)}
  if(!is_empty(IS_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",IS_plots),title=c("Internal Standards Distribution: "),names=IS_names)}
  
  # PCA ----------------------------------------------------------------
  for(i in seq_along(all.comparison.labels)){
    pcas_norm<-paste0(outdir,"/pca/pca_norm_",all.comparison.labels[i],".jpg")
    pcas_raw<-paste0(outdir,"/pca/pca_raw_",all.comparison.labels[i],".jpg")
    
    report<-report%>%add_slide(layout="Comparison",master="Custom Design")
    report<-report%>%ph_with(location=ph_location_type("title"),value=str_c("Comparison ",all.comparison.labels[i]))
    
    if (normalization == "none") {
      if (log2TF) {
        temp1 <- "Log2 Transformed PCA"
      } else if (linearTF) {
        temp1 <- "Linear Transformed PCA"
      } else {
        temp1 <- "Prenormalized PCA"
      }
      temp2 <- "Prenormalized PCA"
    } else {
      temp1 <- "Normalized PCA"
      temp2 <- "Raw PCA"
    }
    report<-report%>%ph_with(location=ph_location_type("body", id=3),value=temp1)
    report<-report%>%ph_with(location=ph_location_type("body", id=1),value=temp2)
    
    if(file.exists(pcas_norm)){
      img<-readJPEG(pcas_norm)
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      report<-report%>%ph_with(value=external_img(pcas_norm),ph_location(left = 5.08,top = 1.56,width = 4.5,height = 4.25/width*height))
    }
    if(file.exists(pcas_raw)){
      img<-readJPEG(pcas_raw)
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      report<-report%>%ph_with(value=external_img(pcas_raw),ph_location(left = 0.5,top = 1.56,width = 4.5,height = 4.25/width*height))
    }
  }
  
  # Volcano plots ----------------------------------------------------------------
  volcanoplots<-list.files(path = paste0(outdir, "/volcano_plots/"), pattern="*.jpg")
  volcanoplots_jpg<-comparison[grepl("_over_", comparison, ignore.case = T)]
  volcanoplots_names<-volcanoplots_jpg
  volcanoplots_jpg<-gsub(x=volcanoplots_jpg,pattern = "\\+","\\\\\\+")
  volcanoplots_jpg<-gsub(x=volcanoplots_jpg,pattern = "\\-","\\\\\\-")
  
  for(i in seq_along(volcanoplots_jpg)){
    volcanoplots_sub <- volcanoplots[grepl(paste0("VolcanoPlot_",volcanoplots_jpg[i]), volcanoplots, ignore.case = T)]
    
    if(!is_empty(volcanoplots_sub)){
      ctrl_ppt(report=report,files=paste0(outdir,"/volcano_plots/",volcanoplots_sub),title=c("Comparison "),names=volcanoplots_names[i])
    }
  }
  
  # Heatmaps ----------------------------------------------------------------
  heatmaps<-list.files(path = paste0(outdir, "/heatmaps/"), pattern="*.jpg")
  heatmaps_jpg<-comparison[!grepl("anova|t-test", comparison, ignore.case = T)]
  heatmaps_names<-heatmaps_jpg
  heatmaps_jpg<-gsub(x=heatmaps_jpg,pattern = "\\+","\\\\\\+")
  heatmaps_jpg<-gsub(x=heatmaps_jpg,pattern = "\\-","\\\\\\-")
  
  for(i in seq_along(heatmaps_jpg)){
    
    heatmaps_sub_0<-heatmaps[grepl(str_c("anova_",heatmaps_jpg[i],"_",settings["padj_method",1],settings["padj_cutoff",1],"|","t-test","_",heatmaps_jpg[i],"_",settings["padj_method",1],settings["padj_cutoff",1]), heatmaps, ignore.case = T)]
    heatmaps_sub_1<-heatmaps[grepl(str_c("anova_",heatmaps_jpg[i],"_",settings["padj_method",1],1,"|","t-test","_",heatmaps_jpg[i],"_",settings["padj_method",1],1), heatmaps, ignore.case = T)]
    
    report<-report%>%add_slide(layout="Comparison",master="Custom Design")
    report<-report%>%ph_with(location=ph_location_type("title"),value=str_c("Comparison ",heatmaps_names[i]))
    report<-report%>%ph_with(location=ph_location_type("body", id=3),value=paste0("Differential ",settings["padj_method",1],"<",settings["padj_cutoff",1]))
    report<-report%>%ph_with(location=ph_location_type("body", id=1),value=paste0("Complete Set ",settings["padj_method",1],"=1"))
    
    if(!is_empty(heatmaps_sub_0)){
      img<-readJPEG(paste0(outdir,"/heatmaps/",heatmaps_sub_0))
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      report<-report%>%ph_with(value=external_img(paste0(outdir,"/heatmaps/",heatmaps_sub_0)),ph_location(left = 5.08,top = 1.56,width = 4.5,height = 4.25/width*height))
    }
    if(!is_empty(heatmaps_sub_1)){
      img<-readJPEG(paste0(outdir,"/heatmaps/",heatmaps_sub_1))
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      report<-report%>%ph_with(value=external_img(paste0(outdir,"/heatmaps/",heatmaps_sub_1)),ph_location(left = 0.5,top = 1.56,width = 4.5,height = 4.25/width*height))
    }
  }
  
  # Save Report -------------------------------------------------------------
  print(report,paste0(outdir,"/report.pptx"))
}