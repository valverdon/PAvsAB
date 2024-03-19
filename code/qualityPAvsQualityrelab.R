#########################
#########################Check link quality P-A vs quality rel.Ab

##################################################################################
####################################PA############################################
##################################################################################

#################data################################################################
load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")#1st version of paper
#need new if rerun
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda")) #1st version of paper
# OTUdata<-readRDS(file=paste0("PA/PR/data/OTUdata_PR.Rds")) #if rerun
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.



load("PA/PR/data/GLM/Eval.Rda")
load("PA/PR/data/GLM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGLM_PR_PA<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM_PR_PA<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

# load(paste0("PA/PR/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/PR/data/GAM/Eval.Rda")
load("PA/PR/data/GAM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGAM_PR_PA<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM_PR_PA<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/GBM/Eval.Rda")
load("PA/PR/data/GBM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM_PR_PA<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM_PR_PA<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/RF/Eval.Rda")
load("PA/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF_PR_PA<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF_PR_PA<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

##################################FU##################################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("PA/FU/data/GLM/Eval.Rda")
load("PA/FU/data/GLM/Fit_Met.Rda")
fitmetGLM_FU_PA<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM_FU_PA<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
# load(paste0("PA/FU/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/FU/data/GAM/Eval.Rda")
load("PA/FU/data/GAM/Fit_Met.Rda")
fitmetGAM_FU_PA<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM_FU_PA<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/GBM/Eval.Rda")
load("PA/FU/data/GBM/Fit_Met.Rda")
fitmetGBM_FU_PA<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM_FU_PA<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/RF/Eval.Rda")
load("PA/FU/data/RF/Fit_Met.Rda")
fitmetRF_FU_PA<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF_FU_PA<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])



##################################BA##################################
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("PA/BA/data/GLM/Eval.Rda")
load("PA/BA/data/GLM/Fit_Met.Rda")

evalmetGLM_BA_PA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA_PA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR_PA<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR_PA<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/BA/data/GAM/Eval.Rda")
load("PA/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA_PA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA_PA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR_PA<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR_PA<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("PA/BA/data/GBM/Eval.Rda")
load("PA/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA_PA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA_PA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR_PA<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR_PA<-fitmet[BAtaxo$Kingdom=="Archaea",]


load("PA/BA/data/RF/Eval.Rda")
load("PA/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA_PA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA_PA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR_PA<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR_PA<-fitmet[BAtaxo$Kingdom=="Archaea",]




#####################################################################################################################################
#########################################################Same AB#####################################################################
#####################################################################################################################################

#######data#############################################################################
################################PR##############################################
load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.

# mins<-apply(OTUdata,2,min)
# maxs<-apply(OTUdata,2,max)

load("AB/PR/data/GLM/Eval_Met.Rda")
load("AB/PR/data/GLM/Fit_Met.Rda")
# evalmet<-evalmet[which(!(is.na(evalmet$Dspear))),]
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGLM_PR_AB<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM_PR_AB<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/GAM/Eval_Met.Rda")
load("AB/PR/data/GAM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGAM_PR_AB<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM_PR_AB<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])




load("AB/PR/data/GBM/Eval_Met.Rda")
load("AB/PR/data/GBM/Fit_Met.Rda")
# evalmet<-evalmet[which(!(is.na(evalmet$Dspear))),]
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM_PR_AB<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM_PR_AB<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/RF/Eval_Met.Rda")
load("AB/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF_PR_AB<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF_PR_AB<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


################################FU##############################################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("AB/FU/data/GLM/Eval_Met.Rda")
load("AB/FU/data/GLM/Fit_Met.Rda")
fitmetGLM_FU_AB<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM_FU_AB<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/FU/data/GAM/Eval_Met.Rda")
load("AB/FU/data/GAM/Fit_Met.Rda")
fitmetGAM_FU_AB<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM_FU_AB<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("AB/FU/data/GBM/Eval_Met.Rda")
load("AB/FU/data/GBM/Fit_Met.Rda")
fitmetGBM_FU_AB<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM_FU_AB<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


load("AB/FU/data/RF/Eval_Met.Rda")
load("AB/FU/data/RF/Fit_Met.Rda")
fitmetRF_FU_AB<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF_FU_AB<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


################################BA-AR###########################################


load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("AB/BA/data/GLM/Eval_Met.Rda")
load("AB/BA/data/GLM/Fit_Met.Rda")
evalmetGLM_BA_AB<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA_AB<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR_AB<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR_AB<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/BA/data/GAM/Eval_Met.Rda")
load("AB/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA_AB<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA_AB<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR_AB<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR_AB<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("AB/BA/data/GBM/Eval_Met.Rda")
load("AB/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA_AB<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA_AB<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR_AB<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR_AB<-fitmet[BAtaxo$Kingdom=="Archaea",]

load("AB/BA/data/RF/Eval_Met.Rda")
load("AB/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA_AB<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA_AB<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR_AB<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR_AB<-fitmet[BAtaxo$Kingdom=="Archaea",]



all(evalmetGLM_BA_PA$OTU==evalmetGLM_BA_AB$OTU)
quality_df_BA<-data.frame(PA_quality=evalmetGLM_BA_PA$TSS_adj,AB_quality=evalmetGLM_BA_AB$Dspear,OTUnumber=evalmetGLM_BA_PA$OTU)

all(evalmetGLM_AR_PA$OTU==evalmetGLM_AR_AB$OTU)
quality_df_AR<-data.frame(PA_quality=evalmetGLM_AR_PA$TSS_adj,AB_quality=evalmetGLM_AR_AB$Dspear,OTUnumber=evalmetGLM_AR_PA$OTU)

all(evalmetGLM_FU_PA$OTU==evalmetGLM_FU_AB$OTU)
quality_df_FU<-data.frame(PA_quality=evalmetGLM_FU_PA$TSS_adj,AB_quality=evalmetGLM_FU_AB$Dspear,OTUnumber=evalmetGLM_FU_PA$OTU)

all(evalmetGLM_PR_PA$OTU==evalmetGLM_PR_AB$OTU)
quality_df_PR<-data.frame(PA_quality=evalmetGLM_PR_PA$TSS_adj,AB_quality=evalmetGLM_PR_AB$Dspear,OTUnumber=evalmetGLM_PR_PA$OTU)

par(mfrow=c(2,2))
plot(quality_df_BA$PA_quality,quality_df_BA$AB_quality,main="bacteria")
plot(quality_df_AR$PA_quality,quality_df_AR$AB_quality,main="archaea")
plot(quality_df_FU$PA_quality,quality_df_FU$AB_quality,main="fungi")
plot(quality_df_PR$PA_quality,quality_df_PR$AB_quality,main="protist")
par(mfrow=c(1,1))
#qualities align


#################################################################################################
######################check alignment of preidcitons#############################################
#################################################################################################


Fitdata_PR_AB<-readRDS("AB/PR/Outputs/GLM/Fit_data/Fit_data_temp41.Rds")
Fitdata_PR_PA<-readRDS("PA/PR/Outputs/GLM/Fit_data/Fit_data_temp41.Rds")

Fitdata_PR_AB<-do.call(cbind, Fitdata_PR_AB)
Fitdata_PR_PA<-do.call(cbind, Fitdata_PR_PA)

Fitdata_PR_AB<-Fitdata_PR_AB[,colnames(Fitdata_PR_AB)[which(colnames(Fitdata_PR_AB)%in%colnames(Fitdata_PR_PA))]]
Fitdata_PR_PA<-Fitdata_PR_PA[,colnames(Fitdata_PR_PA)[which(colnames(Fitdata_PR_PA)%in%colnames(Fitdata_PR_PA))]]
OTUtested<-unique(as.numeric(substr(colnames(Fitdata_PR_AB),4,6)))
for (OTU in OTUtested){
  OTUqual<-c(fitmetGLM_PR_PA$TSS[which(fitmetGLM_PR_PA$OTU==OTU)],fitmetGLM_PR_AB$Dspear[which(fitmetGLM_PR_PA$OTU==OTU)])
  par(mfrow=c(1,2))
  plot(cbind(Fitdata_PR_AB[paste0("OTU",OTU,".fit")],Fitdata_PR_PA[paste0("OTU",OTU,".fit")]),xlab="fitted abundance",ylab="fitted probability of presence", main=paste0("OTU",OTU," TSS=",round(OTUqual[1],3)," Dspear=",round(OTUqual[2],3)))
  plot(cbind(Fitdata_PR_AB[paste0("OTU",OTU,".obs")],Fitdata_PR_PA[paste0("OTU",OTU,".fit")]),xlab="observed abundance",ylab="fitted probability of presence", main=paste0("OTU",OTU," TSS=",round(OTUqual[1],3)," Dspear=",round(OTUqual[2],3)))
  par(mfrow=c(1,1))
}

