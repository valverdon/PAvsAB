
source("code/Evalmetrics.R")
for (PAAB in c("AB")){
  for (GtoM in c("PR")){
    print(GtoM)
    # for (algo in c("GAM","RF","GLM","GBM")){
    
    for (algo in c("RF")){
      Fit_list<-c()
      #algo="RF";PAAB="AB";GtoM="PR"
      print(algo)
      files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_data/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
      files<- gsub("Fit_data","",files)
      for (file in files) { #for each of the files file=files[2]
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_data/Fit_data",file)) #load it
        Metrics<-lapply(Fit_data_list,function(X){Evalmetrics(obs=X$obs[!(is.na(X$fit))],pred=X$fit[!(is.na(X$fit))],PAAB="AB")})
        Metrics<-do.call("rbind", Metrics)
        Fit_list<-rbind(Fit_list,Metrics)  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")

        # if(nrow(explVar_list)!=nrow(Eval_list)){ #nrow(Eval_met_mat) ; nrow(Fit_met_mat)
        #   print(paste0("pb with",file))
        # }
      }
      # nrow(explVar_list) # nrow(Eval_list)  ; rownames(Eval_met_mat) ; rownames(Fit_met_mat)
      fitmet<-as_tibble(Fit_list,.name_repair = "unique")
      fitmet$OTU<-as.numeric(gsub("OTU","",rownames(Fit_list)))
      fitmet <- dplyr::arrange(fitmet,fitmet$OTU)
      save(fitmet,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
    }
  }
}
