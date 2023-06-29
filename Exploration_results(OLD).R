library(tidyverse)
library(gridExtra)
library(reshape2)
library(stringi)
PAAB = "PA"
GtoM="BA"
algo="GLM"

PrepareToPlot <- function(datatab,taxa_level,PAAB){ #datatab=FitEvalP;taxa_level="Kingdom";PAAB="PA"
  library(tidyverse)
  datatab$id <- factor(datatab$id,levels=c("Fit","Eval"))
  datatab <- datatab[!(is.na(datatab[,3])),] 

  if(taxa_level=="Phylum"){   
    datatab$labs <- 1
    labs <- aggregate(labs~Phylum,datatab,function(x){sum(x)/2})
    if (PAAB=="PA"){
      datatab<-reshape2::melt(datatab[c("id","auc","boyce","kap","TSS","SEDI","score","sensit","specif","pospredval","negpredval","accur","rmse","Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")])
      levels(datatab$variable) <- c("AUC","boyce","Kappa","TSS","SEDI","score", "sensitivity", "specificity","PPV", "NPV", "Accuracy", "RMSE")#will change names of the facets in the plot
      
    } else {
      # test<-reshape2::melt(FitEvalP[c("id","R2","D2","MAE","MAEs","MAE_scaled","RMSE","RMSEs", "RMSE_scaled","Dspear","Dpear","Pdispersion")])
      # 
      datatab<-reshape2::melt(datatab[c("id","R2","D2","MAE","MAEs","MAE_scaled","RMSE","RMSEs", "RMSE_scaled","Dspear","Dpear","Pdispersion","Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")])
      levels(datatab$variable) <- c("R2","D2","Mean Absolute Error (MAE)","MAE median_obs rescaled","MAE max_obs rescaled","RMSE", "RMSE median_obs rescaled", "RMSE max_obs rescaled", "Spearman correlation", "Pearson correlation", "P dispersion")#will change names of the facets in the plot
    }
    datatab$value <- sapply(datatab$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
    return(list(datatab,labs))
  } else if (taxa_level=="Kingdom"){
    if (PAAB=="PA"){
      datatab<-reshape2::melt(datatab[c("id","auc","boyce","kap","TSS","SEDI","score","sensit","specif","pospredval","negpredval","accur","rmse")])
      levels(datatab$variable) <- c("AUC","boyce","Kappa","TSS","SEDI","score", "sensitivity", "specificity", "PPV", "NPV", "Accuracy", "RMSE")#will change names of the facets in the plot
      
    } else {
      datatab<-reshape2::melt(datatab[c("id","R2","D2","MAE","MAEs","MAE_scaled","RMSE","RMSEs", "RMSE_scaled","Dspear","Dpear","Pdispersion")])
      levels(datatab$variable) <- c("R2","D2","Mean Absolute Error (MAE)","MAE median_obs rescaled","MAE max_obs rescaled","RMSE", "RMSE median_obs rescaled", "RMSE max_obs rescaled", "Spearman correlation", "Pearson correlation", "P dispersion")#will change names of the facets in the plot
      
    }
    datatab$value <- sapply(datatab$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
    return(datatab)
  }
  
  }

#PA
# for (GtoM in c("PR","BA","FU")){ 
for (GtoM in c("PR","BA","FU")){ 
  print(GtoM)
  load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
  for (algo in c("GLM","GAM","RF")){
    print(algo)
    load(paste0("PA/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
    load(paste0("PA/",GtoM,"/data/",algo,"/Eval_Met.Rda"))

    Fit_met_mat <- fitmet
    Eval_met_mat <- evalmet
    FitEvalP <- bind_rows(list(Fit = Fit_met_mat, Eval= Eval_met_mat), .id="id")
    if (GtoM == "BA"){
      Fit_met_mat_tax <- cbind(fitmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
      Eval_met_mat_tax <- cbind(evalmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
      FitEvalP <- bind_rows(list(Fit = Fit_met_mat_tax, Eval= Eval_met_mat_tax), .id="id")
      FitEvalP_Am<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Archaea",],taxa_level="Kingdom",PAAB="PA")
      p0 <- ggplot(FitEvalP_Am, aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="AUC",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title =element_text(hjust=0.5))
      par(mfrow=c(4,3))
      pdf(paste0("figures/PA_fitevalAR",algo,".pdf"))
      plot(p0)
      dev.off()
      
      FitEvalP_Ba<-PrepareToPlot(FitEvalP[FitEvalP$Kingdom=="Bacteria",],taxa_level="Kingdom",PAAB="PA")
      p0 <- ggplot(FitEvalP_Ba, aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="AUC",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title = element_text(hjust = 0.5))
      pdf(paste0("figures/PA_fitevalBA",algo,".pdf"))
      plot(p0)
      dev.off()
      
      FitEvalP_NA<-PrepareToPlot(FitEvalP[is.na(FitEvalP$Kingdom),],taxa_level="Kingdom",PAAB="PA")
      p0 <- ggplot(FitEvalP_NA, aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_NA[FitEvalP_NA$id=="Fit"&FitEvalP_NA$variable=="AUC",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title = element_text(hjust = 0.5))
      pdf(paste0("figures/PA_fitevalBANA",algo,".pdf"))
      plot(p0)
      dev.off()
      par(mfrow=c(1,1))

    } else{
      FitEvalP_m<-PrepareToPlot(datatab=FitEvalP,taxa_level="Kingdom",PAAB="PA")
      p0 <- ggplot(FitEvalP_m, aes(x=id, y=value))+ geom_violin(scale="width")+
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="AUC",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() +
        facet_wrap(~variable,scales="free") + theme(plot.title = element_text(hjust = 0.5))

      pdf(paste0("figures/PA_fiteval",GtoM,algo,".pdf"))
      par(mfrow=c(4,3))
      plot(p0)
      par(mfrow=c(1,1))
      dev.off()
    }

  }}
    
# Same AB
for (GtoM in c("PR","BA","FU")){ 
  print(GtoM)
  load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
  for (algo in c("GLM","GAM","RF")){
    print(algo)
    load(paste0("Abundance/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
    load(paste0("Abundance/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    
    Fit_met_mat <- fitmet
    Eval_met_mat <- evalmet
    FitEvalP <- bind_rows(list(Fit = Fit_met_mat, Eval= Eval_met_mat), .id="id")
    par(mfrow=c(4,3))
    if (GtoM == "BA"){
      Fit_met_mat_tax <- cbind(fitmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
      # Eval_met_mat_tax <- cbind(evalmet,eval(parse(text=paste0(GtoM,"taxo")))[evalmet$OTU,-1])
      Eval_met_mat_tax <- cbind(evalmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
      FitEvalP <- bind_rows(list(Fit = Fit_met_mat_tax, Eval= Eval_met_mat_tax), .id="id")
      
      FitEvalP_Am<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Archaea",],taxa_level="Kingdom",PAAB="AB")
      FitEvalP_Am$value <- sapply(FitEvalP_Am$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
      p0 <- ggplot(FitEvalP_Am[FitEvalP_Am$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="D2",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_Am[!(FitEvalP_Am$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=id, y=value))+ geom_violin(scale="width")+
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="D2",]))) + 
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + 
        theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(trans="log10") 
      pdf(paste0("figures/AB_fitevalAR",algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()

      FitEvalP_Ba<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Bacteria",],taxa_level="Kingdom",PAAB="AB")
      FitEvalP_Ba$value <- sapply(FitEvalP_Ba$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
      p0 <- ggplot(FitEvalP_Ba[FitEvalP_Ba$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="D2",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_Ba[!(FitEvalP_Ba$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="D2",]))) + 
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + 
        theme(plot.title =element_text(hjust=0.5))+ scale_y_continuous(trans="log10") 
      
      pdf(paste0("figures/AB_fitevalBA",algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()
      
      FitEvalP_NA<-PrepareToPlot(datatab=FitEvalP[is.na(FitEvalP$Kingdom),],taxa_level="Kingdom",PAAB="AB")
      FitEvalP_NA$value <- sapply(FitEvalP_NA$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
      p0 <- ggplot(FitEvalP_NA[FitEvalP_NA$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_NA[FitEvalP_NA$id=="Fit"&FitEvalP_NA$variable=="D2",]))) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + 
        facet_wrap(~variable,scales="free") + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_NA[!(FitEvalP_NA$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_NA[FitEvalP_NA$id=="Fit"&FitEvalP_NA$variable=="D2",]))) + 
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + 
        theme(plot.title =element_text(hjust=0.5))+ scale_y_continuous(trans="log10") 
      
      pdf(paste0("figures/AB_fitevalBANA",algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()
    } else{
      
      FitEvalP_m<-PrepareToPlot(datatab=FitEvalP,taxa_level="Kingdom",PAAB="AB")
      FitEvalP_m$value <- sapply(FitEvalP_m$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
      
      p0 <- ggplot(FitEvalP_m[FitEvalP_m$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="D2",]))) + 
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + 
        theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_m[!(FitEvalP_m$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=id, y=value))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="D2",])))+ 
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + 
        theme(plot.title =element_text(hjust=0.5))+ scale_y_continuous(trans="log10") 

      pdf(paste0("figures/AB_fiteval",GtoM,algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()
      par(mfrow=c(1,1))
    }
  }}

 
### TRAVAUX ###
#add taxo
for (GtoM in c("PR","BA","FU")){ #GtoM ="BA"
  print(GtoM)
  load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
  for (algo in c("GLM","GAM","RF")){#algo="GLM"
    print(algo)
    load(paste0("PA/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
    load(paste0("PA/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    #look at fit
    # nrow(fitmet) ; nrow(eval(parse(text=paste0(GtoM,"taxo")))[-1])
    # nrow(evalmet)
    Fit_met_mat <- cbind(fitmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
    Eval_met_mat <- cbind(evalmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
    FitEvalP <- bind_rows(list(Fit = Fit_met_mat, Eval= Eval_met_mat), .id="id")
    par(mfrow=c(4,3))
    if(GtoM=="BA"){
      #Archeas
      Prepa_Am<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Archaea" & !(is.na(FitEvalP$Kingdom)),],taxa_level="Phylum",PAAB="PA")
      FitEvalP_Am<-Prepa_Am[[1]]
      labs_Am<-Prepa_Am[[2]]
      p0 <- ggplot(FitEvalP_Am, aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="AUC",]))) + 
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Am, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      pdf(paste0("figures/test_phylum/PA_fitevalAR",algo,".pdf"))
      plot(p0)
      dev.off()
      
      #Bacts
      Prepa_Ba<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Bacteria" & !(is.na(FitEvalP$Kingdom)),],taxa_level="Phylum",PAAB="PA")
      FitEvalP_Ba<-Prepa_Ba[[1]]
      labs_Ba<-Prepa_Ba[[2]]
      p1 <- ggplot(FitEvalP_Ba, aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="AUC",]))) + 
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Ba, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      pdf(paste0("figures/test_phylum/PA_fitevalBA",algo,".pdf"))
      plot(p1)
      dev.off()

    } else{
    # p0 <- ggplot(allmolten, aes(x=id, y=value), group=Subphylum)+ geom_violin(scale="width") +labs(title=paste0(GtoM,algo)) +geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + facet_wrap(~variable,scales="free") + theme(axis.text.x = element_text(angle = 90))
      Prepa_m<-PrepareToPlot(datatab=FitEvalP,taxa_level="Phylum",PAAB="PA")
      FitEvalP_m<-Prepa_m[[1]]
      sum(FitEvalP_m$id=="Eval")
      labs_m<-Prepa_m[[2]]
      p1 <- ggplot(FitEvalP_m, aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width") +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="AUC",])))+ 
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_m, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      pdf(paste0("figures/test_phylum/PA_fiteval",GtoM,algo,".pdf"))
      plot(p1)
      dev.off()
    par(mfrow=c(1,1))
    }
  }}

# Same AB
for (GtoM in c("PR","BA","FU")){ 
  print(GtoM)
  load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
  for (algo in c("GLM","GAM","RF")){
    print(algo)
    load(paste0("Abundance/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
    load(paste0("Abundance/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    #look at fit
    Fit_met_mat <- cbind(fitmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])#-1 to remove sequence for readibility
    Eval_met_mat <- cbind(evalmet,eval(parse(text=paste0(GtoM,"taxo")))[-1])
    # ggplot(Fit_met_mat, aes(x=auc,y=auc)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
    # ggplot(Eval_met_mat, aes(x=auc,y=auc)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
    FitEvalP <- bind_rows(list(Fit = Fit_met_mat, Eval= Eval_met_mat), .id="id")
    par(mfrow=c(3,3))
    if(GtoM=="BA"){
      FitEvalP_Am<-FitEvalP[FitEvalP$Kingdom=="Archaea" & !(is.na(FitEvalP$Kingdom)),]
      Prepa_Am<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Archaea" & !(is.na(FitEvalP$Kingdom)),],taxa_level="Phylum",PAAB="AB")
      FitEvalP_Am<-Prepa_Am[[1]]
      labs_Am<-Prepa_Am[[2]]
      p0 <- ggplot(FitEvalP_Am[FitEvalP_Am$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,15)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Am, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_Am[!(FitEvalP_Am$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Am[FitEvalP_Am$id=="Fit"&FitEvalP_Am$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Am, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))+ 
        scale_y_continuous(trans="log10") 
      pdf(paste0("figures/test_phylum/AB_fitevalAR",algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()

      FitEvalP_Ba<-FitEvalP[FitEvalP$Kingdom=="Bacteria" & !(is.na(FitEvalP$Kingdom)),]
      Prepa_Ba<-PrepareToPlot(datatab=FitEvalP[FitEvalP$Kingdom=="Bacteria" & !(is.na(FitEvalP$Kingdom)),],taxa_level="Phylum",PAAB="AB")
      FitEvalP_Ba<-Prepa_Ba[[1]]
      labs_Ba<-Prepa_Ba[[2]]
      p0 <- ggplot(FitEvalP_Ba[FitEvalP_Ba$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,15)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Ba, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_Ba[!(FitEvalP_Ba$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_Ba[FitEvalP_Ba$id=="Fit"&FitEvalP_Ba$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_Ba, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))+ 
        scale_y_continuous(trans="log10") 
      
      pdf(paste0("figures/test_phylum/AB_fitevalBA",algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()
    } else{
      
      Prepa_m<-PrepareToPlot(datatab=FitEvalP[!(is.na(FitEvalP$Kingdom)),],taxa_level="Phylum",PAAB="AB")
      FitEvalP_m<-Prepa_m[[1]]
      labs_m<-Prepa_m[[2]]

      p0 <- ggplot(FitEvalP_m[FitEvalP_m$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation"),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,15)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_m, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))
      p1 <- ggplot(FitEvalP_m[!(FitEvalP_m$variable%in%c("R2","D2","Spearman correlation", "Pearson correlation")),], aes(x=Phylum, y=value, fill=id))+ geom_violin(scale="width")  +
        labs(title=paste0(GtoM," ",algo," quality metrics \n N ASV = ",nrow(FitEvalP_m[FitEvalP_m$id=="Fit"&FitEvalP_m$variable=="D2",]))) +
        geom_boxplot(position=position_dodge(0.9), width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + scale_x_discrete(label = function(x) {str_trunc(x,10)}) + 
        theme(axis.text.x = element_text(angle = 60,vjust=1,hjust=1)) + facet_wrap(~variable,scales="free") + 
        geom_text(aes(x=Phylum, y=Inf, vjust=1, fill=NULL,label=labs),size=3, data = labs_m, check_overlap=TRUE) + theme(plot.title =element_text(hjust=0.5))+ 
        scale_y_continuous(trans="log10") 
      
      pdf(paste0("figures/test_phylum/AB_Fiteval",GtoM,algo,".pdf"))
      plot(p0)
      plot(p1)
      dev.off()
      par(mfrow=c(1,1))
    }
  }}




#RespCtests
load(paste0("PA/",GtoM,"/Outputs/",algo,"/Models.Rda")
load("PR/Outputs/GAMnb/VarImp.Rda")
load("PR/Outputs/GAMnb/Rankings.Rda")
load("PR/data/ENVdata.Rda")
load("PR/data/OTUdata.Rda")
load("PR/data/dataTotSeqSum.Rda")
ASV_VarRank <- ranking_list$OTU1
vartoplot<- ASV_VarRank$var[which(ASV_VarRank$p.value<.05)]
ASV_VarImp <- VarImp_mat[1,]
ASV_RespCtot <- RespCs$OTU1
ABdata <- OTUdata[,"OTU1"]/TotSeqSum

par(mfrow=c(2,2))
for (i in vartoplot){#i = vartoplot[1]
  ASV_RespC<-ASV_RespCtot[,i]
  Realdata <- data.frame(Env=ENVdata[,i],Abundance = ABdata)
  plot(Realdata,ylab="relative abundance",xlab=i)
  ENVgrad <- seq(min(ENVdata[,i]),max(ENVdata[,i]),length.out=100)
  RespCtoplot <- data.frame(ENVgrad=ENVgrad,Abundance=ASV_RespC/median(TotSeqSum))
  points(RespCtoplot,type="l")
  text(min(Realdata$Env),max(Realdata$Abundance)-((max(Realdata$Abundance)-min(Realdata$Abundance))*0.1),paste0("Variable importance = ",round(VarImp_mat[1,i],2),"%\nChi.sq p value = ",round(ASV_VarRank$p.value[which(ASV_VarRank$var==i)],4)),pos=4)

}
par(mfrow=c(1,1))
par(mfrow=c(2,2))
for (i in vartoplot){#i = vartoplot[1]
  ASV_RespC<-ASV_RespCtot[,i]
  Realdata <- data.frame(Env=ENVdata[,i],Abundance = ABdata)
  
  ENVgrad <- seq(min(ENVdata[,i]),max(ENVdata[,i]),length.out=100)
  RespCtoplot <- data.frame(ENVgrad=ENVgrad,Abundance=ASV_RespC/median(TotSeqSum))
  plot(RespCtoplot,type="l",ylab="relative abundance",xlab=i)
  text(min(Realdata$Env),max(Realdata$Abundance)-((max(Realdata$Abundance)-min(Realdata$Abundance))*0.1),paste0("Variable importance = ",round(VarImp_mat[1,i],2),"%\nChi.sq p value = ",round(ASV_VarRank$p.value[which(ASV_VarRank$var==i)],4)),pos=4)
  
}
par(mfrow=c(1,1))
