library(tidyverse)
library(gridExtra)
library(stringi)
#Gathering
PAAB = "PA"
GtoM="BA"
algo="GAM"

#TODEBUG TUESDAY
# Fitlist<-list()
# for (i in 1:length(files)) { #for each of the files file=files[2]
#   load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_Met/",files[i])) #load it
#   Fitlist[[i]]<-unlist(Fit_met_mat)  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")
# }
# which(lapply(Fitlist, ncol)==16)
# Evallist<-list()
# for (i in 1:length(files)) { #for each of the files file=files[2]
#   load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/",files[i])) #load it
#   Evallist[[i]]<-unlist(Eval_met_mat)  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")
# }

#checking Fit metrics

# files <-list.files(path=paste0("Abundance/",GtoM,"/Outputs/",algo,"/Fit_met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
# files <-list(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
# files <- gsub("VarSel","",files)
# 
# for (PAAB in c("PA","AB")){ 
#   print(PAAB)
#   for (GtoM in c("PR","BA","FU")){ 
#     print(GtoM)
#     for (algo in c("GAM","GLM","RF","GBM")){
#       print(algo)
#       # load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",files[1])) #array 1
#       load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel",files[1])) #array 1
#      print(length(VarSel$ranking))
#      print(length(VarSel$preselection))
#       # load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",files[3])) #array 1
#       load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel",files[20])) #array 10
#       print(length(VarSel$ranking))
#       print(length(VarSel$preselection))
#     }
#   }
# }
# #conclusion : preselection gives 15 best variables
# #ranking gives the value of best ones
# #ranking fail : IF OTU"X" does not exist, it had passed the preselection step but no variable were selected for the following


for (PAAB in c("PA","AB")){# PAAB="AB"
print(PAAB)
  for (GtoM in c("PR","BA","FU")){#GtoM="FU"
    # load(paste0(PAAB,"/PR/data/ENVdata.Rda"))
    # colnames(ENVdata)
    # load(paste0(PAAB,"/BA/data/ENVdata.Rda"))
    # colnames(ENVdata)
    # load(paste0(PAAB,"/FU/data/ENVdata.Rda"))
    # colnames(ENVdata)
    print(GtoM)
    load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
    load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
    #patch
    colnames(ENVdata)[which(colnames(ENVdata)=="Elevation")]<-"Altitude"
    
    for (algo in c("GAM","GLM")){
      # algo="GAM"
      seqs<-colnames(OTUdata)

      table_to_fill_rankings<-as.data.frame(matrix(0,ncol=length(ENVdata),nrow=length(seqs)))
      table_to_fill_rankings2<-as.data.frame(matrix(NA,ncol=length(ENVdata),nrow=length(seqs)))
      colnames(table_to_fill_rankings2)<-colnames(table_to_fill_rankings)<-colnames(ENVdata)
      rownames(table_to_fill_rankings2)<-rownames(table_to_fill_rankings)<-paste0("OTU",1:length(seqs))

      table_to_fill_rankings2$seq<-table_to_fill_rankings$seq<-colnames(OTUdata)
      table_to_fill_presels<-table_to_fill_rankings
      
      print(algo)
      files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
      files<- gsub("VarSel","",files)
      #reglage PR champi perdu
      if(GtoM=="PR"){
      load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel_temp244.Rda")) #load last file to check Number of OTU
        if(as.numeric(gsub("OTU","",names(VarSel$preselection)[length(VarSel$preselection)]))!=length(seqs)){
          seqs<-c(colnames(OTUdata)[1:3018],"ChampiPerdu",colnames(OTUdata)[3020:3161])
          table_to_fill_rankings<-as.data.frame(matrix(0,ncol=length(ENVdata),nrow=length(seqs)))
          table_to_fill_rankings2<-as.data.frame(matrix(NA,ncol=length(ENVdata),nrow=length(seqs)))
          colnames(table_to_fill_rankings2)<-colnames(table_to_fill_rankings)<-colnames(ENVdata)
          rownames(table_to_fill_rankings2)<-rownames(table_to_fill_rankings)<-paste0("OTU",1:length(seqs))
          
          table_to_fill_rankings2$seq<-table_to_fill_rankings$seq<-seqs
          table_to_fill_presels<-table_to_fill_rankings

        }
      }
      for (file in files) { #for each of the files 
        # file=files[10]
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel",file)) #load it
          for (i in 1:length(names(VarSel$ranking))){
            #i=12
            # VarSel$ranking[[names(VarSel$ranking)[i]]]["var"]
            linetofill<-names(VarSel$ranking)[i]
            if(algo=="GLM"){

              coltofill<-match(unlist(VarSel$ranking[[names(VarSel$ranking)[i]]]["var"]),colnames(table_to_fill_rankings))


              if(!(all(is.na(coltofill)))){
                if(any(VarSel$ranking[[names(VarSel$ranking)[i]]]["var"]=="Elevation")){
                  VarSel$ranking[[names(VarSel$ranking)[i]]][VarSel$ranking[[names(VarSel$ranking)[i]]]["var"]=="Elevation","var"]<-"Altitude"
                  col_NAvar<-which(colnames(table_to_fill_rankings)%in%c("Altitude","Elevation"))
                  coltofill[which(is.na(coltofill))]<-col_NAvar
                  }
                table_to_fill_rankings[linetofill,coltofill]<-unlist(VarSel$ranking[[names(VarSel$ranking)[i]]]["coef"])
                table_to_fill_rankings2[linetofill,coltofill]<-unlist(VarSel$ranking[[names(VarSel$ranking)[i]]]["rank"])
              }
            }
            if (algo=="GAM"){
              if(!any(is.na(VarSel$ranking[[names(VarSel$ranking)[i]]]))){
                smooths<-VarSel$ranking[[names(VarSel$ranking)[i]]][["smooth"]][VarSel$ranking[[names(VarSel$ranking)[i]]][["smooth"]]["p.value"]<0.05,]
                colnames(smooths)[4]<-"test_val"
                params<-data.frame(VarSel$ranking[[names(VarSel$ranking)[i]]][["param"]],variable=rownames(VarSel$ranking[[names(VarSel$ranking)[i]]][["param"]]))
                params_sel<-params[VarSel$ranking[[names(VarSel$ranking)[i]]][["param"]][,4]<0.05 & names(VarSel$ranking[[names(VarSel$ranking)[i]]][["param"]][,4]<0.05)!="(Intercept)",]

                extract_pvalues<-rbind(smooths[,c(1,4,5)],data.frame(var=rownames(params_sel),test_val=abs(params_sel[,3]),p.value=params_sel[,4]))
                extract_pvalues<-extract_pvalues[order(extract_pvalues$p.value),]
                if (nrow(extract_pvalues)!=0){
                  extract_pvalues$rank<-1:nrow(extract_pvalues)
                  coltofill<-match(extract_pvalues$var,colnames(table_to_fill_rankings))
                  if(any(is.na(coltofill))){
                    prob_var<-extract_pvalues$var[is.na(coltofill)]
                    prob_var_corr<-gsub("1","",prob_var)#pH.1 correction
                    col_NAvar<-match(prob_var_corr,colnames(table_to_fill_rankings))
                    coltofill[which(is.na(coltofill))]<-col_NAvar
                  }
                  if(!all(is.na(coltofill))){
                    table_to_fill_rankings[linetofill,coltofill]<-extract_pvalues$test_val
                    table_to_fill_rankings2[linetofill,coltofill]<-extract_pvalues$rank

                    #if use ranks
                    # table_to_fill_rankings[linetofill,coltofill]<-unlist(VarSel$ranking[[names(VarSel$ranking)[i]]]["rank"])
                  }
                }

            }
          }
        }
        for (i in 1:length(names(VarSel$preselection))){
          #i=1
          VarSel$preselection[[names(VarSel$preselection)[i]]]
          linetofill2<-names(VarSel$preselection)[i]
          coltofill2<-match(VarSel$preselection[[names(VarSel$preselection)[i]]],colnames(table_to_fill_presels))
          if(!all(is.na(VarSel$preselection[[names(VarSel$preselection)[i]]]))){
            col_NAvar2<-which(colnames(table_to_fill_rankings)%in%c("Altitude","Elevation"))
            coltofill2[which(is.na(coltofill2))]<-col_NAvar2
            table_to_fill_presels[linetofill2,coltofill2]<-1
          }
        }
        #   print(paste0("pb with",file))
        # }
      }
      Variable_importance<-table_to_fill_rankings
      # table_to_fill_rankings[758:765,]
      Variable_ranks<-table_to_fill_rankings2
      # table_to_fill_rankings2[758:765,]
      Variable_preselected<-table_to_fill_presels
      # table_to_fill_presels[758:765,]
      save(Variable_importance,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
      save(Variable_ranks,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_ranks.Rda"))
      save(Variable_preselected,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_preselected.Rda"))
    }#end algo

  }#end GtoM
}

# summary(Variable_importance)
# summary(Variable_ranks)
# Variable_preselected[is.na(Variable_preselected$bio1_t),]
# summary(Variable_preselected)
# allpresel<-apply(Variable_preselected[,-which(colnames(Variable_preselected)=="seq")],2,function(X){sum(X)/length(X)})
# allsel
# allranks<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){na.exclude(X)})
#   Var_ranks<-lapply(allranks,as.vector)
#   unlist(Var_ranks,recursive=FALSE)
#   boxplot(Var_ranks,main=paste0(PAAB,"_",GtoM,"_",algo,"_Rank_Distribution"),las = 2)

# for (PAAB in c("PA","AB")){# PAAB="AB"
#   
#   for (GtoM in c("PR","BA","FU")){#GtoM="PR"
#     # load(paste0(PAAB,"/PR/data/ENVdata.Rda"))
#     # colnames(ENVdata)
#     # load(paste0(PAAB,"/BA/data/ENVdata.Rda"))
#     # colnames(ENVdata)
#     # load(paste0(PAAB,"/FU/data/ENVdata.Rda"))
#     # colnames(ENVdata)
#     
#     load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
#     #patch
#     colnames(ENVdata)[which(colnames(ENVdata)=="Elevation")]<-"Altitude"
#     
#     for (algo in c("GAM","GLM")){
for (PAAB in c("PA","AB")){# PAAB="AB"
  print(PAAB)
  for (GtoM in c("PR","BA","FU")){#GtoM="PR"
    print(GtoM)
    load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
    colnames(ENVdata)[which(colnames(ENVdata)=="Elevation")]<-"Altitude"
    load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))


    for (algo in c("GBM","RF")){
      seqs<-colnames(OTUdata)
      # algo="GBM"
      print(algo)
      table_to_fill_rankings<-as.data.frame(matrix(0,ncol=length(ENVdata),nrow=length(seqs)))
      table_to_fill_rankings2<-as.data.frame(matrix(NA,ncol=length(ENVdata),nrow=length(seqs)))
      colnames(table_to_fill_rankings2)<-colnames(table_to_fill_rankings)<-colnames(ENVdata)
      rownames(table_to_fill_rankings2)<-rownames(table_to_fill_rankings)<-paste0("OTU",1:length(seqs))
      
      table_to_fill_rankings2$seq<-table_to_fill_rankings$seq<-colnames(OTUdata)
      table_to_fill_presels<-table_to_fill_rankings
      
      print(algo)
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"importance_list.Rda")) #load last file to check Number of OTU¨
      
      files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
      files<- gsub("VarSel","",files)
      
      if(GtoM=="PR"){
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel_temp244.Rda")) #load last file to check Number of OTU
        if(as.numeric(gsub("OTU","",names(VarSel$preselection)[length(VarSel$preselection)]))!=length(seqs)){
          seqs<-c(colnames(OTUdata)[1:3018],"ChampiPerdu",colnames(OTUdata)[3020:3161])
          table_to_fill_rankings<-as.data.frame(matrix(0,ncol=length(ENVdata),nrow=length(seqs)))
          table_to_fill_rankings2<-as.data.frame(matrix(NA,ncol=length(ENVdata),nrow=length(seqs)))
          colnames(table_to_fill_rankings2)<-colnames(table_to_fill_rankings)<-colnames(ENVdata)
          rownames(table_to_fill_rankings2)<-rownames(table_to_fill_rankings)<-paste0("OTU",1:length(seqs))
          
          table_to_fill_rankings2$seq<-table_to_fill_rankings$seq<-seqs
          table_to_fill_presels<-table_to_fill_rankings
          
        }
      }
      for (file in files) { #for each of the files 
        # file=files[10]
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel",file)) #load it
        importancetolook<-importance_list[names(VarSel$preselection)]
        
        for (i in 1:length(importancetolook)){
          #i=1
          # VarSel$ranking[[names(VarSel$ranking)[i]]]["var"]
          linetofill<-names(importancetolook)[i]
          coltofill<-match(rownames(importancetolook[[i]]),colnames(table_to_fill_rankings))
          #####################################################################################################
          if(any(is.na(coltofill))){ ############ ICI a mod
            if("Elevation"%in%rownames(importancetolook[[i]])&"Altitude"%in%colnames(table_to_fill_rankings)){
              coltofill[which(is.na(coltofill))]<-which(colnames(table_to_fill_rankings)=="Altitude")
            }
              
          }
          ###############################################################
          if(!(all(is.na(coltofill)))){
            
            
            if(algo=="RF"){
              table_to_fill_rankings[linetofill,coltofill]<-importancetolook[[i]][,1]
              ranking<-importancetolook[[i]][,1][order(importancetolook[[i]][,1],decreasing=TRUE)]
              ranking2<-1:length(ranking)
              ranking2[which(ranking==0)]<-NA
              coltofill2<-match(names(ranking),colnames(table_to_fill_rankings))
              if("Elevation"%in%names(ranking)&"Altitude"%in%colnames(table_to_fill_rankings)){
                coltofill2[which(is.na(coltofill2))]<-which(colnames(table_to_fill_rankings)=="Altitude")
              }
              table_to_fill_rankings2[linetofill,coltofill2]<-ranking2
            }
            if (algo=="GBM"){
              table_to_fill_rankings[linetofill,coltofill]<-importancetolook[[i]][,2]
              coltofill2<-match(rownames(importancetolook[[i]]),colnames(table_to_fill_rankings))
              ranking2<-1:length(importancetolook[[i]][,2])
              ranking2[which(importancetolook[[i]][,2]==0)]<-NA
              table_to_fill_rankings2[linetofill,coltofill2]<-importancetolook[[i]][,2]
            }


          }
        }
        for (i in 1:length(names(VarSel$preselection))){
          #i=1
          VarSel$preselection[[names(VarSel$preselection)[i]]]
          linetofill3<-names(VarSel$preselection)[i]
          coltofill3<-match(VarSel$preselection[[names(VarSel$preselection)[i]]],colnames(table_to_fill_presels))
          if(!(any(is.na(coltofill3)))){
            table_to_fill_presels[linetofill3,coltofill3]<-1
          }
        }
        #   print(paste0("pb with",file))
        # }
      }
  
      Variable_importance<-table_to_fill_rankings
      Variable_ranks<-table_to_fill_rankings2
      Variable_preselected<-table_to_fill_presels
      save(Variable_importance,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
      save(Variable_ranks,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_ranks.Rda"))
      save(Variable_preselected,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_preselected.Rda"))
    }
  }
}

# #checking
# for (PAAB in c("PA","AB")){
#   for (GtoM in c("PR","BA","FU")){
#     print(GtoM)
#     for (algo in c("GAM","RF","GLM","GBM")){
#       print(algo)
#       load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
#       print(nrow(fitmet))
#       load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
#       print(nrow(evalmet))
#     }}}
