
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


#for GLM
library(data.table)
library(stringi)
library(glmnet)
library(MASS)
library(plyr)
library(randomForest)
library(RRF)
# 
# arrayID=1
# GtoM="PR"
# Modeltype <- "All"
# PAAB <- "both"
#get options from slurm
args <- commandArgs(trailingOnly = TRUE)
arrayID= as.numeric(args[[1]]) 
GtoM <- args[[2]] 
Modeltype <- args[[3]]
PAAB <- args[[4]] 

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

RepCpoints <- 100

ENVmedian <- as.data.frame(matrix(apply(ENVdata,2,function(x){median(as.numeric(as.character(x)))})
                               ,nrow=RepCpoints, ncol=ncol(ENVdata),byrow=TRUE))
colnames(ENVmedian) <- c(colnames(ENVdata))

#######################################################################################################################################
#######################################################################################################################################
########################################################   P - A   ####################################################################
#######################################################################################################################################
#######################################################################################################################################

if(PAAB=="PA"|PAAB=="both"){
  load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))

  endloop <- ceiling(ncol(OTUdata)/250) #endthisloop=10

  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#207 on last loop
  if (endthisloop<=0) { q("no")} #already finished
  

  ###############################################################################################################################
  ####################################################### GAM ###################################################################
  ###############################################################################################################################
  # Time1=Sys.time()
  if(Modeltype=="All"|Modeltype=="GAM"){
    load(paste0("PA/",GtoM,"/Outputs/GAM/Models/Models_temp",arrayID,".Rda"))
    
    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=1
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      Mgam <- Model_list[[paste0("OTU",OTUtoRun)]]

      if(is.na(Mgam[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################
        
        # save(VarSel, file=paste0(GtoM,"/Outputs/VarSel_OTU", OTUtoRun, ".Rda")) #testing line
        
        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(attr(Mgam$terms,"term.labels")))#pb for NA models
        colnames(RespC) <- attr(Mgam$terms,"term.labels")
        for (covariate in colnames(RespC)){ #covariate=colnames(RespC)[1]
          ENVmed <- ENVmedian
          #create a gradient covering the whole environment for that covariate
          if(covariate %in% contvar){#covariate continuous
            gradient <- seq(min(ENVdata[,colnames(ENVdata)==covariate]),max(ENVdata[,colnames(ENVdata)==covariate]),length.out=RepCpoints)
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict.gam(Mgam, newdata=ENVmed, type="response")
          } else {#covariate binary
            gradient <- c(rep(levels(ENVdata[,colnames(ENVdata)==covariate])[1],RepCpoints/2),rep(levels(ENVdata[,colnames(ENVdata)==covariate])[2],RepCpoints/2))
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict.gam(Mgam, newdata=ENVmed, type="response")
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))){
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
    
  }#End GAM
  
  ###############################################################################################################################
  ####################################################### GLM ###################################################################
  ###############################################################################################################################
  if(Modeltype=="All"|Modeltype=="GLM"){
    load(paste0("PA/",GtoM,"/Outputs/GLM/Models/Models_temp",arrayID,".Rda"))
    
    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=1
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      Mglm <- Model_list[[paste0("OTU",OTUtoRun)]]

      
      if(is.na(Mglm[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################
        covariate_temp<-gsub(", 2\\)","",gsub("poly\\(","",rownames(Mglm$glmnet.fit$beta[-1,])))
        covariates<-unique(substring(covariate_temp,1,nchar(covariate_temp)-1))
        form<-as.formula(paste0("OTUdata[,OTUtoRun] ~ ", paste0(c(paste0("poly(",covariates[covariates%in%contvar],",2)"),paste0(covariates[covariates%in%binvar])),collapse=" + ")))

        ModMat <- model.matrix(form,data = ENVdata)
        ModMatmedian <- matrix(apply(ModMat,2,function(x){median(as.numeric(as.character(x)))}),nrow=RepCpoints,ncol=ncol(ModMat), byrow = TRUE)
        colnames(ModMatmedian) <- colnames(ModMat)
        
        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(covariates))
        colnames(RespC) <- covariates
        for (covariate in 1:length(covariates)){ #covariate=1
          ModMatR <- ModMatmedian#remove intercept
          #create a gradient covering the whole environment for that covariate
          if(covariates[covariate] %in% contvar){#covariate continuous 
            ModMatR[,(covariate*2):(covariate*2+1)]<-cbind(seq(min(ModMat[,(covariate*2)]),max(ModMat[,(covariate*2)]),length.out=RepCpoints),seq(min(ModMat[,(covariate*2+1)]),max(ModMat[,(covariate*2+1)]),length.out=RepCpoints))

            #predict on those fake point
            if(Mglm$lambda[1]!=Mglm$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.1se, type="response")
            } else{#take lambda that minimize error
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.min, type="response")
            }
          } else {#covariate binary
            ModMatR[,paste0(covariates[covariate],"1")]<-c(rep(min(ModMat[,paste0(covariates[covariate],"1")]),RepCpoints/2),rep(max(ModMat[,paste0(covariates[covariate],"1")]),RepCpoints/2))
            if(Mglm$lambda[1]!=Mglm$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.1se, type="response")
            } else{#take lambda that minimize error
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.min, type="response")
            }
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))){
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
    
  }#end GLM
  ###############################################################################################################################
  #######################################################  RF  ###################################################################
  ###############################################################################################################################
  if(Modeltype=="All"|Modeltype=="RF"){
    load(paste0("PA/",GtoM,"/Outputs/RF/Models/Models_temp",arrayID,".Rda"))
    
    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=4
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      MRF <- Model_list[[paste0("OTU",OTUtoRun)]]

      if(is.na(MRF[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################
        covariates<-attributes(MRF$terms)$term.labels

        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(covariates))
        colnames(RespC) <- covariates
        
        for (covariate in colnames(RespC)){ #covariate=colnames(RespC)[1]
          ENVmed <- ENVmedian
          #create a gradient covering the whole environment for that covariate
          if(covariate%in% contvar){#covariate continuous 
            gradient <- seq(min(ENVdata[,colnames(ENVdata)==covariate]),max(ENVdata[,colnames(ENVdata)==covariate]),length.out=RepCpoints)
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict(MRF, newdata = ENVmed, type = "prob")[,2]
          
          } else {#covariate binary
            gradient <- c(rep(levels(ENVdata[,colnames(ENVdata)==covariate])[1],RepCpoints/2),rep(levels(ENVdata[,colnames(ENVdata)==covariate])[2],RepCpoints/2))
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict(MRF, newdata = ENVmed, type = "prob")[,2]
            
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))){
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
    
  }
}#end PA


#######################################################################################################################################
#######################################################################################################################################
########################################################  ABUNDANCE ###################################################################
#######################################################################################################################################
#######################################################################################################################################

if(PAAB=="AB"|PAAB=="both"){
  load(paste0("Abundance/",GtoM,"/data/OTUdata.Rda"))
  load(paste0("Abundance/",GtoM,"/data/dataTotSeqSum.Rda"))
  
  endloop <- ceiling(ncol(OTUdata)/250) #endthisloop=5
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#207 on last loop
  if (endthisloop<=0) { q("no")} #already finished


  
  ###############################################################################################################################
  ####################################################### GAM ###################################################################
  ###############################################################################################################################
  
  if(Modeltype=="All"|Modeltype=="GAMnb"){
    load(paste0("Abundance/",GtoM,"/Outputs/GAMnb/Models/Models_temp",arrayID,".Rda"))
    
    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      Mgam <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(is.na(Mgam[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################

        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(attr(Mgam$terms,"term.labels")))#pb for NA models
        colnames(RespC) <- attr(Mgam$terms,"term.labels")
        ENVmed <- as.data.frame(matrix(c(apply(ENVdata,2,function(x){median(as.numeric(as.character(x)))}),median(TotSeqSum))
                                       ,nrow=RepCpoints, ncol=ncol(ENVdata)+1,byrow=TRUE))
        colnames(ENVmed) <- c(colnames(ENVdata),"TotSeqSum")
        for (covariate in colnames(RespC)){ #covariate=colnames(RespC)[1]
          #create a gradient covering the whole environment for that covariate
          
          if(covariate %in% contvar){
            gradient <- seq(min(ENVdata[,colnames(ENVdata)==covariate]),max(ENVdata[,colnames(ENVdata)==covariate]),length.out=RepCpoints)
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict.gam(Mgam, newdata=ENVmed, type="response")
          } else {
            gradient <- c(rep(levels(ENVdata[,colnames(ENVdata)==covariate])[1],RepCpoints/2),rep(levels(ENVdata[,colnames(ENVdata)==covariate])[2],RepCpoints/2))
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict.gam(Mgam, newdata=ENVmed, type="response")
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAMnb/RespCs"))){
      save(RespC_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAMnb/RespCs"))
      save(RespC_list, file=paste0("Abundance/",GtoM,"/Outputs/GAMnb/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
    
  }#End GAM
  
  ###############################################################################################################################
  ####################################################### GLM ###################################################################
  ###############################################################################################################################
  if(Modeltype=="All"|Modeltype=="GLM"){
    load(paste0("Abundance/",GtoM,"/Outputs/GLM/Models/Models_temp",arrayID,".Rda"))

    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=1
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      Mglm <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(is.na(Mglm[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################
        covariate_temp<-gsub(", 2\\)","",gsub("poly\\(","",rownames(Mglm$glmnet.fit$beta[-1,])))
        covariates<-unique(substring(covariate_temp,1,nchar(covariate_temp)-1))
        form<-as.formula(paste0("OTUdata[,OTUtoRun] ~ ", paste0(c(paste0("poly(",covariates[covariates%in%contvar],",2)"),paste0(covariates[covariates%in%binvar])),collapse=" + ")))
        
        ModMat <- model.matrix(form,data = ENVdata)
        ModMatmedian <- matrix(apply(ModMat,2,function(x){median(as.numeric(as.character(x)))}),nrow=RepCpoints,ncol=ncol(ModMat), byrow = TRUE)
        colnames(ModMatmedian) <- colnames(ModMat)
        TotSeqSumMed <- median(TotSeqSum)
        
        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(covariates))
        colnames(RespC) <- covariates
        
        for (covariate in 1:length(covariates)){ #covariate=1
          #create a gradient covering the whole environment for that covariate
          ModMatR <- ModMatmedian#remove intercept
          
          if(covariates[covariate] %in% contvar){
            ModMatR[,(covariate*2):(covariate*2+1)]<-cbind(seq(min(ModMat[,(covariate*2)]),max(ModMat[,(covariate*2)]),length.out=RepCpoints),seq(min(ModMat[,(covariate*2+1)]),max(ModMat[,(covariate*2+1)]),length.out=RepCpoints))
            #predict on those fake point
            if(Mglm$lambda[1]!=Mglm$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.1se, type="response", newoffset=log(TotSeqSumMed))
            } else{#take lambda that minimize error
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.min, type="response", newoffset=log(TotSeqSumMed))
            }
          } else {
            ModMatR[,paste0(covariates[covariate],"1")]<-c(rep(min(ModMat[,paste0(covariates[covariate],"1")]),RepCpoints/2),rep(max(ModMat[,paste0(covariates[covariate],"1")]),RepCpoints/2))
            if(Mglm$lambda[1]!=Mglm$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.1se, type="response", newoffset=log(TotSeqSumMed))
            } else{#take lambda that minimize error
              RespC[,covariates[covariate]] <- predict(Mglm,newx=ModMatR, s=Mglm$lambda.min, type="response", newoffset=log(TotSeqSumMed))
            }
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GLM/RespCs"))){
      save(RespC_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GLM/RespCs"))
      save(RespC_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
  }#End GLM 
  
  ####################################################################################################################################
  ###########################################################  RF  ###################################################################
  ####################################################################################################################################
  if(Modeltype=="All"|Modeltype=="RF"){
    load(paste0("Abundance/",GtoM,"/Outputs/RF/Models/Models_temp",arrayID,".Rda"))
    
    RespC_list<- list()
    
    #for each species
    for (i in 1:endthisloop){ #i=1
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      MRF <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(is.na(MRF[1])){#dont bother do RespC if no model
        RespC_list[[paste0("OTU",OTUtoRun)]]<- matrix(NA,nrow=RepCpoints,ncol=15)#pb for NA models
      } else{
        
        ##############################################################################################################
        #### P4 : Get response curves ##################################################################
        covariates<-attributes(MRF$terms)$term.labels
        form<-as.formula(paste0("factor(OTUtrain)~ ", paste(covariates,collapse= " + ")))
        
        RespC <- matrix(NA,nrow=RepCpoints,ncol=length(covariates))
        colnames(RespC) <- covariates
        
        
        for (covariate in 1:length(covariates)){ #covariate=1
          ENVmed <- ENVmedian
          #create a gradient covering the whole environment for that covariate
          if(covariates[covariate] %in% contvar){#covariate continuous 
            gradient <- seq(min(ENVdata[,colnames(ENVdata)==covariate]),max(ENVdata[,colnames(ENVdata)==covariate]),length.out=RepCpoints)
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict(MRFEval, newdata = ENVmed)
            
          } else {#covariate binary
            gradient <- c(rep(levels(ENVdata[,colnames(ENVdata)==covariate])[1],RepCpoints/2),rep(levels(ENVdata[,colnames(ENVdata)==covariate])[2],RepCpoints/2))
            ENVmed[,colnames(ENVmed)==covariate] <- gradient
            #predict on those fake point
            RespC[,covariate] <- predict(MRFEval, newdata = ENVmed)
            
          }
        }
        RespC_list[[paste0("OTU",OTUtoRun)]]<-RespC
      }
    }#end ASV loop
    #Save RespC results
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))){
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/RespCs"))
      save(RespC_list, file=paste0("PA/",GtoM,"/Outputs/GLM/RespCs/RespCs_temp", arrayID, ".Rda"))
    }
    
  }
}#end AB

q("no")
