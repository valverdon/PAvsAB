
Fit_GAM<-function(PAAB,arrayID,GtoM,OTUdata,ENVdata,test=FALSE){
  library(mgcv)
  #test=TRUE
  # arrayID=100
  # GtoM="PR"
  # PAAB <- "AB"
  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  # load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
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
  nthreads <- 2
  registerDoParallel(cores = nthreads)
  

  ENVdataforcor<-ENVdata
  ENVdataforcor[binvar]<-apply(ENVdataforcor[binvar], 2, function(x) as.numeric(as.character(x)))
  var_cor <- cor(ENVdataforcor, use = "pairwise.complete.obs")
  
  endloop <- ceiling(ncol(OTUdata)/250)
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop} 
  if (endthisloop<=0) { q("no")} #if the array don't contain anything (finished), quit R
  #endthisloop=5
  
  UnivExpl_list <- list()
  preselected_list <- list()
  ranking_list <- list()
  Fit_met_list <- list()
  Fit_data_list <- list()
  Deviances_list <- list()
  Model_list <- list()

  loop<-1:endthisloop
  #set lists that will be filled with results.
  if(test){
    if(GtoM!="PR"){
      loop<-sample(1:endthisloop,10)
    }
    time1<-Sys.time()
  }
for (i in loop){ #i=1
  # for (i in 1:2){ print(i)

  OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  
  if(test){
    print(OTUtoRun)
  }
  weights <- sapply(OTUdata[,OTUtoRun],function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
  #weight of present = nsite/nsite_where_present.
  #weight of absent = nsite/nsite_where_absent.
  explVar_values <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
  
  #OTUtoRun<-3
  if(PAAB=="PA"){
  if (sum(OTUdata[,OTUtoRun])/length(OTUdata[,OTUtoRun])>0.95) {#dont bother model it if no absence data (more than 95% presence over all plots)
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
    # Method : Fit  univariate polynomial2 GLMs and take covariates from the 15 best models (avoiding colinear covariates)
    
    #not sure if using weights is good or not yet
    
    # for each covariate, fit an univariate model (polynomial 2nd degree)
    for (j in names(ENVdata)) { #j=names(ENVdata)[27]
      if(j %in% contvar){
        tryCatch(
          {
            # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4), family="binomial", weights = weights, na.action="na.omit", control=list(nthreads=nthreads, maxit=500))
            model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVdata[,j], 2), family="quasibinomial", weights = weights, na.action="na.omit")
          }, error=function(cond){
            message(paste0("error with variable ",j," of OTU ", OTUtoRun))
            return(NA)
          })
      } else { #handle binary variable
        model <- glm(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]), family="quasibinomial", weights = weights, na.action="na.omit")
      }
      if (all(is.na(model))){ #if fit failed, drop that variable
        explVar_values[j] <- 0 
      } else { #if everything when right, store variance explianed by the model of that variable
        # model <- bam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
        explVar_values[j] <- (summary(model)$null.deviance-summary(model)$deviance)/summary(model)$null.deviance*100 
      }
    }
    # store all var expl of univariate models for that ASV
    UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values  
    
    #b) preselection of variables, 
    # Method take the most important and remove everything correlated, then take second most important in whats left ...
    thresh <- 0.8 #correlation threshold
    
    #get a new object with correlation matrix 
    varcor <- var_cor
    preselected <- c()
    for (k in 1:15){
      explVar_values <- sort(explVar_values,decreasing=T) #sort variable by their explanatory power
      preselected <- c(preselected,names(explVar_values)[1]) #Get the name of the best var
      tokeep <- names(which(abs(varcor[,names(explVar_values)[1]]) < thresh)) #names of variable uncorrelated with selected one;  
      #names(which(abs(varcor[,names(explVar_values)[1]]) > 0.7))
      varcor <- varcor[tokeep,tokeep] # reduce the corr matrix, keeping only the variable not selected and not correlated with selected ones
      explVar_values <- explVar_values[tokeep] #reduce the explVar_values vector the same way
    }
    #store the 15 selected "best" variables for that ASV
    preselected_list[[paste0("OTU",OTUtoRun)]]  <- preselected
    
    ##############################################################################################################
    #### P2 : Model calibration ##################################################################################
    #### a: Fit the full model  #########################################################################
    # Method : Null-space penalization (shrinkage method) selection.
    
    #get model formula
    form<-as.formula(paste0 ("OTUdata[,OTUtoRun]~ ", paste(c(paste0("s(",preselected[preselected%in%contvar],",bs='cs',k=3)"),paste0(preselected[preselected%in%binvar])),collapse=" + ")))
    
    #fit model
    Mgam  <- tryCatch(
      {
        mgcv::gam(form, data=ENVdata, family="binomial", weights = weights,select=TRUE, control=list(nthreads=nthreads, maxit=500))
      }, error=function(cond){
        message(paste0("error fitting OTU ", OTUtoRun))
        message("Original error:")
        message(cond)
        return(NA)
      }
    )
    if(all(is.na(Mgam))){ #if fit failed set all results to NA
      ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=NA,param=NA)
      Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
      Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA))
      Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
      Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
      next
    } else{
      # Extract results
      #variables ranked by their "importance"
      gam.beta<-data.frame(var=names(Mgam$model)[-(c(1:(1+length(preselected[preselected%in%binvar])),17))], summary(Mgam)$s.table, row.names = NULL)  #[-(1:2)] if offset in model form
      rank_list_cont <- data.frame(gam.beta[order(gam.beta$p.value),], rank = 1:nrow(gam.beta), model="gam")
      rank_list_bin <- summary(Mgam)$p.table
      ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=rank_list_cont,param=rank_list_bin)
      
      #deviance explained
      Fit_GAM <- summary(Mgam)$dev.expl*100 #Deviance explained percentage (close to R2)
      devs <- c(null.deviance = Mgam$null.deviance, mod.deviance = Mgam$deviance)
      Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs 
      
      #table with data and prediction of the model
      fitted.val<-predict.gam(Mgam, newdata = ENVdata, type = "response")
      pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=fitted.val)
      
      #compute evaluation metrics
      Metrics <- Evalmetrics(obs=pred_expl$obs,pred=pred_expl$fit,PAAB="PA")
      Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=Fit_GAM, Metrics)
      Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
      Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
    } #end what to do with fitted model
  }#end if not enough absence
}#end PA
  
  if(PAAB=="AB"){
    load(paste0(PAAB,"/",GtoM,"/data/OTUdata.Rda"))
    load(paste0(PAAB,"/",GtoM,"/data/dataTotSeqSum.Rda")) #total read count data (to go to relative abundance)
    
    #a) fit univariate models
    
    for (j in names(ENVdata)) { #j=names(ENVdata)[27]
      if(j %in% contvar){
        tryCatch(
          {
            # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4), family="binomial", weights = weights, na.action="na.omit", control=list(nthreads=nthreads, maxit=500))
            #I added an offset of the model to convert read count into read relative abundance (log because it has to be the same scale as the response variable that is in log for poisson regression)
            model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVdata[,j], 2) + offset(log(TotSeqSum)), family="poisson", weights = weights, na.action="na.omit")
          }, error=function(cond){
            message(paste0("error with variable ",j," of OTU ", OTUtoRun))
            return(NA)
          })
      } else {
        model <- glm(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]) + offset(log(TotSeqSum)), family="poisson", weights = weights, na.action="na.omit")
      }
      if (all(is.na(model))){
        explVar_values[j] <- 0 
      } else {
        explVar_values[j] <- (summary(model)$null.deviance-summary(model)$deviance)/summary(model)$null.deviance*100 
      }
    }
    UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values 
    
    #b ) Fit univariate models 
    
    thresh <- 0.8
    
    #get a new object with correlation matrix 
    varcor <- var_cor
    preselected <- c()
    #preselection for all OTUs 1 by one, take the most important and remove everything correlated, then take second most important in whats left ...
    for (k in 1:15){
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
    #offset added to switch to rel. abundance
    form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste(c(paste0("s(",preselected[preselected%in%contvar],",bs='cs',k=3)"),paste0(preselected[preselected%in%binvar])),collapse=" + "),"+ offset(log(TotSeqSum))"))
    
    Mgam  <- tryCatch(
      {
        mgcv::gam(form, data=ENVdata, family="poisson", weights = weights, select=TRUE,control=list(nthreads = nthreads,maxit=500))
      }, error=function(cond){
        message(paste0("error fitting OTU ", OTUtoRun))
        message("Original error:")
        message(cond)
        return(NA)
      }
    )
    if(all(is.na(Mgam[1]))){ 
      ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=NA,param=NA)
      Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
      Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
      Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
      Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
      next
    } else{
      # Extract results
      gam.beta<-data.frame(var=names(Mgam$model)[-(c(1:(1+1+length(preselected[preselected%in%binvar])),18))], summary(Mgam)$s.table, row.names = NULL)  #[-(1+1)]for response and offset
      rank_list_cont <- data.frame(gam.beta[order(gam.beta$p.value),], rank = 1:nrow(gam.beta), model="gam")
      rank_list_bin <- summary(Mgam)$p.table
      ranking_list[[paste0("OTU",OTUtoRun)]]<- list(smooth=rank_list_cont,param=rank_list_bin)
      Fit_GAM <- summary(Mgam)$dev.expl*100 #Deviance explained percentage (close to R2)
      devs <- c(null.deviance = Mgam$null.deviance, mod.deviance = Mgam$deviance)
      Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs
      
      pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=10^(Mgam$fitted.values))
      
      #compute evaluation metrics
      Metrics <- Evalmetrics(pred_expl$obs,pred=pred_expl$fit,PAAB="AB")
      Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=Fit_GAM, Metrics)
      Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
      Model_list[[paste0("OTU",OTUtoRun)]] <- Mgam
    }
  }# end AB 
}#end this ASV
  if(test){
    time2<-Sys.time()
    print(time2-time1)
  }
  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  #Univariate models variance explained
  UnivExpl_mat <- do.call("rbind", UnivExpl_list)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/UnivExpl"))){
    save(UnivExpl_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM"))
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/UnivExpl"))
    save(UnivExpl_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/UnivExpl/UnivExpl_temp", arrayID, ".Rda"))
  }
  #Selected Variables and shrinkage
  VarSel <- list(preselection=preselected_list,ranking=ranking_list)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/VarSel"))){
    save(VarSel, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/VarSel/VarSel_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/VarSel"))
    save(VarSel, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/VarSel/VarSel_temp", arrayID, ".Rda"))
  }
  #Model Deviances
  Deviances_mat <- do.call("rbind", Deviances_list)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Deviances"))){
    save(Deviances_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Deviances/Deviances_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Deviances"))
    save(Deviances_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Deviances/Deviances_temp", arrayID, ".Rda"))
  }
  #Model calibration metrics
  Fit_met_mat <- do.call("rbind", Fit_met_list)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_met"))){
    save(Fit_met_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_met/Fit_met_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_met"))
    save(Fit_met_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_met/Fit_met_temp", arrayID, ".Rda"))
  }
  #save fit and cross validation models data (if needed)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_data"))){
    save(Fit_data_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_data/Fit_data_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_data"))
    save(Fit_data_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Fit_data/Fit_data_temp", arrayID, ".Rda"))
  }
  #Save model fit
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Models"))){
    save(Model_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Models/Models_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GAM/Models"))
    save(Model_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GAM/Models/Models_temp", arrayID, ".Rda"))
  }
  return(Model_list)
}
