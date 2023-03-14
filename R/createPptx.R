createPptx <- function(project_title = "Project Report",
                       project_subtitle = "Metabolomics Analysis",
                       exprsdf = exprsdf,
                       settings = settings,
                       all.comparisons = all.comparisons,
                       outdir = outdir) {
  
  library(officer)
  library(tidyverse)
  library(jpeg)
  
  normalization <- settings["normalization",1]
  log2tr <- as.logical(settings["log2transform",1])
  
  sample_number<-nrow(exprsdf[["norm"]]) #Not including the quality control samples
  liver_sample_number<-exprsdf$qc_num #"Quality Control" samples
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
    if (log2tr) {
      norm_method<-c("Data was log2 transformed")
    } else {
      norm_method<-c("Data was prenormalized") 
    }
  }
  
  if (all(grepl("_over_", all.comparisons, ignore.case = T))) {
    comparison <- c("t-test",
                    all.comparisons)
    comparison_levels <- c(2, rep(3, length(all.comparisons)))
  } else {
    myflag <- grepl("_over_", all.comparisons, ignore.case = T)
    comparison <- c("t-test",
                    all.comparisons[myflag],
                    "anova",
                    all.comparisons[!(myflag)])
    comparison_levels <- c(2, rep(3, length(all.comparisons[myflag])),
                           2, rep(3, length(all.comparisons[!(myflag)])))
  }
  
  #define where the template powerpoint is located
  report <- read_pptx("data/Project_Report_Template.pptx")
  
  # Functions ---------------------------------------------------------------
  ctrl_ppt<-function(report,files,title,names) {
    for(i in seq_along(files)){
      img<-readJPEG(files[i])
      height<-dim(img)[1]/96
      width<-dim(img)[2]/96
      
      report<-report%>%add_slide(layout="Title and Content", master="Custom Design")
      report<-report%>%ph_with(value=str_c(title,": ",names[i]), location=ph_location_type("title"))
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
  if (normalization == "istd") {
    ul_levels <- c(1,2,2,
                   1,rep(2,number_of_methods),
                   1,comparison_levels)
    ul_str <- c("Samples",
                paste0(sample_number," experimental samples"),
                paste0(liver_sample_number," liver control samples"),
                paste0("Total of ",number_of_species," features from ",number_of_methods," methods"),
                norm_method,
                "Comparisons performed",
                comparison)
  } else {
    ul_levels <- c(1,2,
                   1,rep(2,number_of_methods),
                   1,comparison_levels)
    ul_str <- c("Samples",
                paste0(sample_number," experimental samples"),
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
  liver_plots<-qa[grepl("liver_",qa,ignore.case = T)]
  liver_names<-str_replace(liver_plots,"liver_dist_","")%>%str_replace(".jpg","")
  IS_plots<-qa[grepl("dIS",qa)]
  IS_names<-str_replace(IS_plots,"dIS_dist_","")%>%str_replace(".jpg","")
  
  
  if(!is_empty(qa_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",qa_plots),title=c("QA plots"),names=qa_names)}
  if(!is_empty(liver_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",liver_plots),title=c("Liver Control"),names=liver_names)}
  if(!is_empty(IS_plots)){ctrl_ppt(report=report,files=paste0(outdir,"/QA_plots/",IS_plots),title=c("Internal Standards Distribution"),names=IS_names)}
  
  
  # Heatmaps ----------------------------------------------------------------
  heatmaps<-list.files(path = paste0(outdir, "/heatmaps/"), pattern="*.jpg")
  heatmaps_jpg<-comparison[!grepl("anova|t-test", comparison, ignore.case = T)]
  heatmaps_names<-heatmaps_jpg
  heatmaps_jpg<-gsub(x=heatmaps_jpg,pattern = "\\+","\\\\\\+")
  heatmaps_jpg<-gsub(x=heatmaps_jpg,pattern = "\\-","\\\\\\-")
  
  for(i in seq_along(heatmaps_jpg)){
    
    heatmaps_sub_0<-heatmaps[grepl(str_c("anova_",heatmaps_jpg[i],"_",settings["statistic_selector",1],settings["statistic_cutoff",1],"|","t-test","_",heatmaps_jpg[i],"_",settings["statistic_selector",1],settings["statistic_cutoff",1]), heatmaps, ignore.case = T)]
    heatmaps_sub_1<-heatmaps[grepl(str_c("anova_",heatmaps_jpg[i],"_",settings["statistic_selector",1],1,"|","t-test","_",heatmaps_jpg[i],"_",settings["statistic_selector",1],1), heatmaps, ignore.case = T)]
    
    report<-report%>%add_slide(layout="Comparison",master="Custom Design")
    report<-report%>%ph_with(location=ph_location_type("title"),value=str_c("Comparison ",heatmaps_names[i]))
    report<-report%>%ph_with(location=ph_location_type("body", id=3),value="Differential FDR<0.25")
    report<-report%>%ph_with(location=ph_location_type("body", id=1),value="Complete Set FDR=1")
    
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