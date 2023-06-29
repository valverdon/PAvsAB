# Extract Varimp from RF and GBM models
PAAB = "AB"
GtoM="BA"
algo="GBM"
library(RRF)
library(gbm)
library(tidyverse)
# pathtowork<-"../../work/FAC/FBM/DEE/aguisan/sometalp/vverdon/"
# pathtowork<-""



# for (PAAB in c("PA","AB")){
  # for(GtoM in c("PR","FU","BA")){
for (algo in c("RF","GBM")){
  
  files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Models/Models"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
  files <- gsub("Models","",files)
  
#RRF
load(paste0(PAAB,"/",GtoM,"/data/OTUdata.Rda"))
importance_list<-list()

for (file in files) { #for each of the files
# file=files[5]
  print(file)
  load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Models/Models/Models",file)) #load it
  print(length(Model_list))
  for (i in 1:length(names(Model_list))){
    if (algo=="GBM"){
      if(all(is.na(Model_list[[names(Model_list[i])]]))){
        importance_list[[names(Model_list[i])]]<-NA
      }else{
    importance_list[[names(Model_list[i])]]<-summary.gbm(Model_list[[names(Model_list[i])]],plotit=FALSE)
    }}
    if(algo=="RF"){
      if(all(is.na(Model_list[[names(Model_list[i])]]))){
        importance_list[[names(Model_list[i])]]<-NA
      }else{
      importance_list[[names(Model_list[i])]]<-RRF::importance(Model_list[[names(Model_list[i])]])
      }
    }
    # i=1
}}
# names(importance_list)
#all have 1 
# all(table(as.numeric(gsub("OTU","",names(importance_list))))==1)

save(importance_list,file=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"importance_list.Rda"))

    }