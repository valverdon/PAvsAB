
########### WORK IN PROGRESS

###################################################################################################
# Full pipe for P-A fit and evaluation. Called from cluster using arguments                       # 
# GtoM for microbial group and  Modeltype for Model type                                          #
# e.g. sbatch pipe_core.txt PR GAMnb     for protists modeled by neg binom                        #
##################################  called from PA folder     #####################################
.libPaths("/work/FAC/FBM/DEE/aguisan/default/rlib/4.0")

library(parallel)
library(mgcv)
library(tidyverse)

#for evalmetrics
library(modEvA)
library(ecospat)
library(ROCR)

#get options from slurm
args <- commandArgs(trailingOnly = TRUE)
arrayID= as.numeric(args[[1]]) #arrayID=57
GtoM <- args[[2]] #GtoM="BA"
Modeltype <- args[[3]] #Modeltype <- "GAMnb"
PAAB <- args[[4]] #PAAB <- "both"


#functions I wrote that will be used
source("code/Evalmetrics.R")#Function thqt compute some metrics of model evqluqtion
source("code/create_DataSplitTable.R") #function to prepqre cross validation framework


nthreads <- detectCores()-1
set.seed(seed=08022022) #date of the day I added a seed

load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))

binvar <- c()
contvar <- c()
for (i in names(ENVdata)){#i =names(ENVdata) [31]
  if(length(unique(ENVdata[,i]))>2){
    contvar <- c(contvar,i)
  } else {
    binvar <- c(binvar,i)
  }
}



#######################################################################################################################
################################################ PA ###################################################################
#######################################################################################################################
if(PAAB=="PA"|PAAB=="both"){
  load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
  
  endloop <- ceiling(ncol(OTUdata)/250) #endthisloop=5
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#207 on last loop
  if (endthisloop<=0) { q("no")} #already finished

  
  #############################################################################################################
  ##############################################   GAM   ######################################################
  ##############################################   nb    ######################################################
  
  if(Modeltype=="All"|Modeltype=="GAMnb"){
    
    UnivExpl_list <- list()
    preselected_list <- list()
    ranking_list <- list()
    Fit_met_list <- list()
    Fit_data_list <- list()
    Deviances_list <- list()
    Model_list <- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      
      if (sum(OTUdata[,OTUtoRun])/length(OTUdata[,OTUtoRun])>0.95) {#dont bother model it if no absence data, Predict Presence everywhere
        UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
        preselected_list[[paste0("OTU",OTUtoRun)]]  <- rep(NA,15)
        ranking_list[[paste0("OTU",OTUtoRun)]] <- list(smooth=NA,param=NA)
        Metrics <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, Metrics)
        Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- data.frame(obs=OTUdata[,OTUtoRun],fit=rep(1,length(OTUdata[,OTUtoRun])))
        Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
        Model_list[[paste0("OTU",OTUtoRun)]] <- NA
      } else{
        
        ##############################################################################################################
        #### Part1 : preselection of covariates #########################################################################
        # Method : Fit univariate models : Take covariates from the 15 best models (avoiding colinear covariates)
        
        #a) fit univariate models
        explVar_values <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
        
        for (j in names(ENVdata)) { #j=names(ENVdata)[31]
          if(j %in% contvar){
            tryCatch(
              {
                model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
              }, error=function(cond){
                message(paste0("error with variable ",j," of OTU ", OTUtoRun))
                return(NA)
              })
          } else {
            model <- gam(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
          }
          if (is.na(model[[1]])){explVar_values[j] <- 0 }
          else {
            # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
            explVar_values[j] <- summary(model)$dev.expl*100 
          }
        }
        UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values  
        
        #b) preselection for all OTUs 1 by one, take the most important and remove everything correlated, then take second most important in whats left ...
        thresh <- 0.8
        #correlation matrix of envdata predictors
        varcor <- cor(ENVdata, use = "pairwise.complete.obs")
        
        preselected <- c()
        for (i in 1:15){
          explVar_values <- sort(explVar_values,decreasing=T)
          preselected <- c(preselected,names(explVar_values)[1]) #Get the name of the most important var
          tokeep <- names(which(abs(varcor[,names(explVar_values)[1]]) < thresh)) #keep non correlated  names(which(abs(varcor[,names(explVar_values)[1]]) > 0.7))
          varcor <- varcor[tokeep,tokeep] # reduce the cor matrix
          explVar_values <- explVar_values[tokeep] #reduce the explVar_values vector
        }
        preselected_list[[paste0("OTU",OTUtoRun)]]  <- preselected
        
        
        ##############################################################################################################
        #### P2 : Model calibration ##################################################################################
        #### a: Fit the full model  #########################################################################
        # Method : Null-space penalization (shrinkage method) selection.
        
        form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste(c(paste0("s(",preselected[preselected%in%contvar],",bs='cs',k=4)"),paste0(preselected[preselected%in%binvar])),collapse=" + ")))
        Mgam  <- tryCatch(
          {
            mgcv::gam(form, data=ENVdata, family="binomial",select=TRUE,nthreads = nthreads)
          }, error=function(cond){
            message(paste0("error fitting OTU ", OTUtoRun))
            message("Original error:")
            message(cond)
            return(NA)
          }
        )
        if(is.na(Mgam[1])){ #warning message if no problem beacause objet with several elements
          ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=NA,param=NA)
          Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
          Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA))
          Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
          Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
          Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
          Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
          RespC_list[[paste0("OTU",OTUtoRun)]]<-NA
          next
        } else{
          # Extract results
          gam.beta<-data.frame(var=names(Mgam$model)[-(1:(1+length(preselected[preselected%in%binvar])))], summary(Mgam)$s.table, row.names = NULL)  #[-(1:2)] if offset in model form
          rank_list_cont <- data.frame(gam.beta[order(gam.beta$p.value),], rank = 1:nrow(gam.beta), model="gam")
          rank_list_bin <- summary(Mgam)$p.table
          ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=rank_list_cont,param=rank_list_bin)
          
          Fit_GAMnb <- summary(Mgam)$dev.expl*100 #Deviance explained percentage (close to R2)
          devs <- c(null.deviance = Mgam$null.deviance, mod.deviance = Mgam$deviance)
          Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs
          
          pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=Mgam$fitted.values)
          
          #compute evaluation metrics
          Metrics <- Evalmetrics(pred_expl$obs,pred=pred_expl$fit)
          Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=Fit_GAMnb, Metrics)
          Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
          Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
        } #end what to do with fitted model
      }#end if not enough absence
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    
    #Univariate models variance explained
    UnivExpl_mat <- do.call("rbind", UnivExpl_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/UnivExpl"))){
      save(UnivExpl_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb"))
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/UnivExpl"))
      save(UnivExpl_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
    }
    #Selected Variables and shrinkage
    VarSel <- list(preselection=preselected_list,ranking=ranking_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/VarSel"))){
      save(VarSel, file=paste0("PA/",GtoM,"/Outputs/GAMnb/VarSel/VarSel_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/VarSel"))
      save(VarSel, file=paste0("PA/",GtoM,"/Outputs/GAMnb/VarSel/VarSel_temp", arrayID, ".Rda"))
    }
    #Model Deviances
    Deviances_mat <- do.call("rbind", Deviances_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/Deviances"))){
      save(Deviances_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Deviances/Deviances_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/Deviances"))
      save(Deviances_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Deviances/Deviances_temp", arrayID, ".Rda"))
    }
    #Model calibration metrics
    Fit_met_mat <- do.call("rbind", Fit_met_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_met"))){
      save(Fit_met_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_met/Fit_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_met"))
      save(Fit_met_mat, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_met/Fit_met_temp", arrayID, ".Rda"))
    }
    #save fit and cross validation models data (if needed)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_data"))){
      save(Fit_data_list, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_data/Fit_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_data"))
      save(Fit_data_list, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Fit_data/Fit_data_temp", arrayID, ".Rda"))
    }
    #Save model fit
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAMnb/Models"))){
      save(Model_list, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Models/Models_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAMnb/Models"))
      save(Model_list, file=paste0("PA/",GtoM,"/Outputs/GAMnb/Models/Models_temp", arrayID, ".Rda"))
    }
  }#end GAMnb
}#end if PAAB

#######################################################################################################################
################################################ AB ###################################################################
#######################################################################################################################
if(PAAB=="AB"|PAAB=="both"){
  load(paste0("Abundance/",GtoM,"/data/OTUdata.Rda"))
  load(paste0("Abundance/",GtoM,"/data/dataTotSeqSum.Rda"))
  endloop <- ceiling(ncol(OTUdata)/250)
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#207 on last loop
  if (endthisloop<=0) { q("no")} #already finished
  
  #############################################################################################################
  ##############################################   GAM   ######################################################
  ##############################################   nb    ######################################################
  if(Modeltype=="All"|Modeltype=="GAMnb"){
    UnivExpl_list <- list()
    preselected_list <- list()
    ranking_list <- list()
    Fit_met_list <- list()
    Fit_data_list <- list()
    Deviances_list <- list()
    Model_list <- list()
    
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      print(i)
      #From 04_
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      
      ##############################################################################################################
      #### Part1 : preselection of covariates #########################################################################
      # Method : Fit univariate models : Take covariates from the 15 best models (avoiding colinear covariates)
      
      #a) fit univariate models
      explVar_values <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
      
      for (j in names(ENVdata)) { #j=names(ENVdata)[31]
        if(j %in% contvar){
          tryCatch(
            {
              model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4) + offset(log(TotSeqSum)), family="nb", na.action="na.omit", control=list(nthreads=nthreads))
            }, error=function(cond){
              message(paste0("error with variable ",j," of OTU ", OTUtoRun))
              return(NA)
            })
        } else {
          model <- gam(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]) + offset(log(TotSeqSum)), family="nb", na.action="na.omit", control=list(nthreads=nthreads))
        }
        if (is.na(model[[1]])){explVar_values[j] <- 0 }#if cant create the model, bad covariate to not select
        else {
          # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
          explVar_values[j] <- summary(model)$dev.expl*100 
        }
      }
      UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values  
      
      #b ) Fit univariate models 
      
      thresh <- 0.8
      #correlation matrix of envdata predictors
      varcor <- cor(ENVdata, use = "pairwise.complete.obs")
      
      preselected <- c()
      #preselection for all OTUs 1 by one, take the most important and remove everything correlated, then take second most important in whats left ...
      for (i in 1:15){
        explVar_values <- sort(explVar_values,decreasing=T)
        preselected <- c(preselected,names(explVar_values)[1]) #Get the name of the most important var
        tokeep <- names(which(abs(varcor[,names(explVar_values)[1]]) < thresh)) #keep non correlated  names(which(abs(varcor[,names(explVar_values)[1]]) > 0.7))
        varcor <- varcor[tokeep,tokeep] # reduce the cor matrix
        explVar_values <- explVar_values[tokeep] #reduce the explVar_values vector
      }
      preselected_list[[paste0("OTU",OTUtoRun)]]  <- preselected
      
      
      ##############################################################################################################
      #### P2 : Model calibration ##################################################################################
      
      #### a: Fit the model  #########################################################################
      
      # Method : Null-space penalization (shrinkage method) selection.
      
      form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste(c(paste0("s(",preselected[preselected%in%contvar],",bs='cs',k=4)"),paste0(preselected[preselected%in%binvar])),collapse=" + "),"+ offset(log(TotSeqSum))"))
      Mgam  <- tryCatch(
        {
          mgcv::gam(form, data=ENVdata, family="nb",select=TRUE,nthreads = nthreads)
        }, error=function(cond){
          message(paste0("error fitting OTU ", OTUtoRun))
          message("Original error:")
          message(cond)
          return(NA)
        }
      )
      if(is.na(Mgam[1])){ #warning message if no problem beacause objet with several elements
        ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=NA,param=NA)
        Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
        Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA))
        Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
        Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        RespC_list[[paste0("OTU",OTUtoRun)]]<-NA
        next
      } else{
        # Extract results
        gam.beta<-data.frame(var=names(Mgam$model)[-(1:(1+1+length(preselected[preselected%in%binvar])))], summary(Mgam)$s.table, row.names = NULL)  #[-(1+1)]for response and offset
        rank_list_cont <- data.frame(gam.beta[order(gam.beta$p.value),], rank = 1:nrow(gam.beta), model="gam")
        rank_list_bin <- summary(Mgam)$p.table
        ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=rank_list_cont,param=rank_list_bin)
        
        Fit_GAMnb <- summary(Mgam)$dev.expl*100 #Deviance explained percentage (close to R2)
        devs <- c(null.deviance = Mgam$null.deviance, mod.deviance = Mgam$deviance)
        Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs
        
        pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=Mgam$fitted.values)
        
        #compute evaluation metrics
        Metrics <- Evalmetrics(pred_expl$obs,pred=pred_expl$fit)
        Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=Fit_GAMnb, Metrics)
        Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
        Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
      }
    }# end ASV loop
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    
    #Univariate models variance explained
    UnivExpl_mat <- do.call("rbind", UnivExpl_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/UnivExpl"))){
      save(UnivExpl_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb"))
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/UnivExpl"))
      save(UnivExpl_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
    }
    #Selected Variables and shrinkage
    VarSel <- list(preselection=preselected_list,ranking=ranking_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/VarSel"))){
      save(VarSel, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/VarSel/VarSel_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/VarSel"))
      save(VarSel, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/VarSel/VarSel_temp", arrayID, ".Rda"))
    }
    #Model Deviances
    Deviances_mat <- do.call("rbind", Deviances_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Deviances"))){
      save(Deviances_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Deviances/Deviances_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Deviances"))
      save(Deviances_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Deviances/Deviances_temp", arrayID, ".Rda"))
    }
    #Model calibration metrics
    Fit_met_mat <- do.call("rbind", Fit_met_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_met"))){
      save(Fit_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_met/Fit_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_met"))
      save(Fit_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_met/Fit_met_temp", arrayID, ".Rda"))
    }
    #save fit and cross validation models data (if needed)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_data"))){
      save(Fit_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_data/Fit_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_data"))
      save(Fit_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Fit_data/Fit_data_temp", arrayID, ".Rda"))
    }
    #Save model fit
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Models"))){
      save(Model_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Models/Models_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Models"))
      save(Model_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/Models/Models_temp", arrayID, ".Rda"))
    }
  }#End GAMnb
}#end if PAAB

q("no")