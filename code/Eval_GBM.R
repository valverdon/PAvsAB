
Eval_GBM<-function(PAAB,arrayID,GtoM,ENVdata,OTUdata,OTUstoRun, Model_list, NbRunEval=100, DataSplit=0.8, validation.method = "bootstrap",CCV=CCV,test=FALSE){
  #Model_list<-Modlist
  #for gam
  library(gbm)
  library(parallel)
  library(doParallel)
  # 
  # arrayID=1
  # GtoM="NULL_PR"
  # PAAB <- "PA"
  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  # load(paste0(PAAB,"/",GtoM,"/Outputs/GBM/Models/Models_temp",arrayID,".Rda"))
  # NbRunEval=100
  # validation.method = "bootstrap"
  
  binvar <- c()
  contvar <- c()
  for (covar in names(ENVdata)){#i =names(ENVdata) [31]
    if(length(unique(ENVdata[,covar]))>2){
      contvar <- c(contvar,covar)
    } else {
      binvar <- c(binvar,covar)
      ENVdata[,covar]<-as.factor(ENVdata[,covar])
    }
  }
  
  #functions I wrote that will be used
  #Function I wrote that compute some metrics of model evaluation
  #part of the P-A metrics code are from Antoine A.'s work
  source("code/Evalmetrics.R")
  source("code/create_DataSplitTable.R") #function to split dataset for the CV procedure.
  
  nthreads <- 2
  registerDoParallel(cores = nthreads)
  #test=TRUE
  #endthisloop=5
  Eval_met_list <- list()
  Pred_data_list <- list()
  Eval_values_list <-list()
  if(test){
    time1<-Sys.time()
  }
  
  for (OTUtoRun in OTUstoRun){#OTUtoRun=OTUstoRun[1]
    # for (i in 1:2){ print(i)
    # OTUtoRun <- as.numeric(substring(names(Model_list)[i],first=4))
    if(test){
      print(OTUtoRun)
    }
    Mod <- Model_list[[paste0("OTU",OTUtoRun)]]
    
    # if(validation.method=="LOO"){
    #   NbRunEval <- nrow(OTUdata)
    # }
    if(CCV){ #if community cross validation, get common CV setup, else create CV setup for this OTU.
      dataSplits<-readRDS(paste0("data/splitOTUdata_",validation.method,"_",NbRunEval,"splits_",group,".Rds"))
    }else{
      dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = validation.method, DataSplit=DataSplit)
    }

    
    #OTUtoRun<-3
    if(PAAB=="PA"){
      if(all(is.na(Mod))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        pred_GBM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        VarSel<-readRDS(paste0("data/Outputs/GBM/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_",arrayID,".rds"))
        preselected<-VarSel$preselection[[paste0("OTU",OTUtoRun)]]

        ##############################################################################################################
        #### Part1 : preselection of covariates #########################################################################
        # Method : Fit  univariate polynomial2 GLMs and take covariates from the 15 best models (avoiding colinear covariates)
        
        #not sure if using weights is good or not yet
        
        # for each covariate, fit an univariate model (polynomial 2nd degree)
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          # weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
          ENVeval[binvar] <- lapply(ENVeval[binvar],factor)
          ENVeval[contvar] <- lapply(ENVeval[contvar],as.numeric)
          
          formCV <- as.formula(paste0("OTUtrain~ ", paste(preselected,collapse= " + ")))
          #################################

          ModEval <- tryCatch(
            { gbm(formCV, data=as.data.frame(ENVtrain), distribution="bernoulli", n.cores= nthreads, n.trees=Mod$n.trees,shrinkage=Mod$shrinkage,bag.fraction=Mod$bag.fraction)
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")"
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          if(!(is.na(ModEval[1]))){
            predCV <- predict.gbm(ModEval, n.trees=Mod$n.trees, newdata = ENVeval,type="response")
            pred_GBM[rownames(ENVeval), paste0("pred",f)] <- predCV

            if(validation.method!="LOO"){
              tryCatch(
                { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUeval,pred=predCV,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
                }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                  message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                  message(paste("boot will be ignored"))
                  return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
                }
              )
            }
          }
        }#end that split
        if(validation.method=="LOO"){
          pred_GLM2<-apply(pred_GLM,1,function(X){mean(X,na.rm=TRUE)})
          tryCatch(
            { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[,OTUtoRun],pred=pred_GLM2,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
              message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
              message(paste("boot will be ignored"))
              return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
            }
          )
        }
        
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GBM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }#end else
    }#PA
    
    # if(PAAB=="AB"){
    #   load(paste0(PAAB,"/",GtoM,"/data/dataTotSeqSum.Rda")) #total read count data (to go to relative abundance)
    #   
    #   
    #   if (is.na(Mod[1])){
    #     pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
    #     
    #     Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GAM)
    #     Eval_met_list[[paste0("OTU",OTUtoRun)]] <- rep(NA, 15)
    #     next
    #   } else {
    #     
    #     pred_GBM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
    #     Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion","R2_scaled","MAE_scaled", "RMSE_scaled", "Dspear_scaled", "Dpear_scaled", "Pdisp_scaled")))
    #     
    #     load(paste0(PAAB,"/",GtoM,"/Outputs/GBM/VarSel/VarSel_temp",arrayID,".Rda"))
    #     preselected<-VarSel$preselection[[OTUtoRun_name]]
    #     
    #     for (f in 1:ncol(dataSplits)) { #f=2 for the different selection of traning/eval plots.
    #       # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
    #       # ENVtrain <- ENVdata[dataSplits[,f], ]
    #       # Offset <- log(TotSeqSum[dataSplits[,f]])
    #       # print(f)
    #       OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
    #       ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
    #       TotSeqSumtrain <- rep(TotSeqSum,times=dataSplits[,f])
    #       
    #       weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
    #       
    #       ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
    #       ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
    #       ENVtrain <- cbind(ENVtrain,TotSeqSumtrain=TotSeqSumtrain)
    #       
    #       ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
    #       OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
    #       TotSeqSumEval <- TotSeqSum[dataSplits[,f]==FALSE]
    #       
    #       formCV <- as.formula(paste0("OTUtrain~ ", paste(preselected,collapse= " + "),"+ offset(log(TotSeqSumtrain))"))
    # 
    #       #################################
    #       ModEval <- tryCatch(
    #         { gbm(formCV, data=ENVtrain, distribution="poisson", weights = weights, n.cores= nthreads, n.trees=Mod$n.trees,shrinkage=Mod$shrinkage,bag.fraction=Mod$bag.fraction)
    #         }, error=function(cond){
    #           message(paste0("error with boot ",f," of OTU ", OTUtoRun))
    #           message(cond)
    #           message(paste("boot will be ignored"))
    #           return(NA)
    #         }
    #       )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
    #       #                data=ENVtrain, family="poisson")
    #       # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
    #       #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
    #       
    #       if(!(is.na(ModEval[1]))){
    #         predCV <- predict.gbm(ModEval, n.trees=Mod$n.trees, newdata = ENVeval,type="response")+log(TotSeqSumEval)
    #         pred_GBM[rownames(ENVeval), paste0("pred",f)] <- predCV
    # 
    #         if(validation.method!="LOO"){
    #           tryCatch(
    #             { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUeval,pred=fitted.val,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
    #             }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
    #               message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
    #               message(paste("boot will be ignored"))
    #               return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
    #             }
    #           )
    #         }
    #       }
    #     }
    #     if(validation.method=="LOO"){
    #       pred_GLM2<-apply(pred_GBM,1,function(X){mean(X,na.rm=TRUE)})
    #       tryCatch(
    #         { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[,OTUtoRun],pred=pred_GLM2,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
    #         }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
    #           message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
    #           message(paste("boot will be ignored"))
    #           return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
    #         }
    #       )
    #     }
    #     
    #     pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
    #     
    #     Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GBM)
    #     Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
    #     Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
    #   }#end else
    # }#end AB
  }#end ASV
  
  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  
  Eval_met_mat <- do.call("rbind", Eval_met_list)

  
  if(CCV){
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Eval_met"))){
      saveRDS(Eval_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met/Eval_met_CCV_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Eval_met"))
      save(Eval_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met/Eval_met_CCV_temp_", arrayID, ".Rds"))
    }
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Pred_data"))){
      save(Pred_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Pred_data/Pred_data_CCV_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Pred_data"))
      save(Pred_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Pred_data/Pred_data_CCV_temp_", arrayID, ".Rds"))
    }
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/Eval_met_CCV_allboot_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/Eval_met_CCV_allboot_temp_", arrayID, ".Rds"))
    }
    if(test){
      return(Eval_met_mat<-Eval_met_mat)
    }
  }else{
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Eval_met"))){
      saveRDS(Eval_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met/Eval_met_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Eval_met"))
      save(Eval_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met/Eval_met_temp_", arrayID, ".Rds"))
    }
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Pred_data"))){
      save(Pred_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Pred_data/Pred_data_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Pred_data"))
      save(Pred_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Pred_data/Pred_data_temp_", arrayID, ".Rds"))
    }
    if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/Eval_met_allboot_temp_", arrayID, ".Rds"))
    } else {
      dir.create(paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("data/Outputs/GBM/",PAAB,"/Eval_met_allboot/Eval_met_allboot_temp_", arrayID, ".Rds"))
    }
    if(test){
      return(Eval_met_mat<-Eval_met_mat)
    }
  }
}