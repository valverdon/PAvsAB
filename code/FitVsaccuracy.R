
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
fitmetGLM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

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
fitmetGAM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/GBM/Eval.Rda")
load("PA/PR/data/GBM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/RF/Eval.Rda")
load("PA/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

##################################FU##################################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("PA/FU/data/GLM/Eval.Rda")
load("PA/FU/data/GLM/Fit_Met.Rda")
fitmetGLM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
# load(paste0("PA/FU/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/FU/data/GAM/Eval.Rda")
load("PA/FU/data/GAM/Fit_Met.Rda")
fitmetGAM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/GBM/Eval.Rda")
load("PA/FU/data/GBM/Fit_Met.Rda")
fitmetGBM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/RF/Eval.Rda")
load("PA/FU/data/RF/Fit_Met.Rda")
fitmetRF_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])



##################################BA##################################
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("PA/BA/data/GLM/Eval.Rda")
load("PA/BA/data/GLM/Fit_Met.Rda")

evalmetGLM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/BA/data/GAM/Eval.Rda")
load("PA/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("PA/BA/data/GBM/Eval.Rda")
load("PA/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


load("PA/BA/data/RF/Eval.Rda")
load("PA/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


#########plotting#######################################################################

png(file=paste0("figures/PAAB_selection/testmodels/fitVSaccuracy_PA.png"),res=300,width=1961,height=2000)
par(mfrow=c(4,4))
plot(fitmetGLM_BA$TSS[!fitmetGLM_BA$TSS<=0],evalmetGLM_BA$TSS[!fitmetGLM_BA$TSS<=0],xlab="",ylab="model quality (TSS)",main="Bacteria",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_AR$TSS[!fitmetGLM_AR$TSS<=0],evalmetGLM_AR$TSS[!fitmetGLM_AR$TSS<=0],xlab="",ylab="",main="Archaea",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_FU$TSS[!fitmetGLM_FU$TSS<=0],evalmetGLM_FU$TSS[!fitmetGLM_FU$TSS<=0],xlab="",ylab="",main="Fungi",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_PR$TSS[!fitmetGLM_PR$TSS<=0],evalmetGLM_PR$TSS[!fitmetGLM_PR$TSS<=0],xlab="",ylab="",main="Protist",pch=16,cex=.2)
abline(0,1)


plot(fitmetGAM_BA$TSS[!fitmetGAM_BA$TSS<=0],evalmetGAM_BA$TSS[!fitmetGAM_BA$TSS<=0],xlab="",ylab="model quality (TSS)",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_AR$TSS[!fitmetGAM_AR$TSS<=0],evalmetGAM_AR$TSS[!fitmetGAM_AR$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_FU$TSS[!fitmetGAM_FU$TSS<=0],evalmetGAM_FU$TSS[!fitmetGAM_FU$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_PR$TSS[!fitmetGAM_PR$TSS<=0],evalmetGAM_PR$TSS[!fitmetGAM_PR$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)

plot(fitmetGBM_BA$TSS[!fitmetGBM_BA$TSS<=0],evalmetGBM_BA$TSS[!fitmetGBM_BA$TSS<=0],xlab="",ylab="model quality (TSS)",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_AR$TSS[!fitmetGBM_AR$TSS<=0],evalmetGBM_AR$TSS[!fitmetGBM_AR$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_FU$TSS[!fitmetGBM_FU$TSS<=0],evalmetGBM_FU$TSS[!fitmetGBM_FU$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_PR$TSS[!fitmetGBM_PR$TSS<=0],evalmetGBM_PR$TSS[!fitmetGBM_PR$TSS<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)

plot(fitmetRF_BA$TSS[!fitmetRF_BA$TSS<=0],evalmetRF_BA$TSS[!fitmetRF_BA$TSS<=0],xlab="model fit(TSS)",ylab="model quality (TSS)",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_AR$TSS[!fitmetRF_AR$TSS<=0],evalmetRF_AR$TSS[!fitmetRF_AR$TSS<=0],xlab="model fit(TSS)",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_FU$TSS[!fitmetRF_FU$TSS<=0],evalmetRF_FU$TSS[!fitmetRF_FU$TSS<=0],xlab="model fit(TSS)",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_PR$TSS[!fitmetRF_PR$TSS<=0],evalmetRF_PR$TSS[!fitmetRF_PR$TSS<=0],xlab="model fit(TSS)",ylab="",pch=16,cex=.2)
abline(0,1)
par(mfrow=c(1,1))
dev.off()


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
fitmetGLM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/GAM/Eval_Met.Rda")
load("AB/PR/data/GAM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGAM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])




load("AB/PR/data/GBM/Eval_Met.Rda")
load("AB/PR/data/GBM/Fit_Met.Rda")
# evalmet<-evalmet[which(!(is.na(evalmet$Dspear))),]
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/RF/Eval_Met.Rda")
load("AB/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF_PR<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


################################FU##############################################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("AB/FU/data/GLM/Eval_Met.Rda")
load("AB/FU/data/GLM/Fit_Met.Rda")
fitmetGLM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/FU/data/GAM/Eval_Met.Rda")
load("AB/FU/data/GAM/Fit_Met.Rda")
fitmetGAM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("AB/FU/data/GBM/Eval_Met.Rda")
load("AB/FU/data/GBM/Fit_Met.Rda")
fitmetGBM_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


load("AB/FU/data/RF/Eval_Met.Rda")
load("AB/FU/data/RF/Fit_Met.Rda")
fitmetRF_FU<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


################################BA-AR###########################################


load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("AB/BA/data/GLM/Eval_Met.Rda")
load("AB/BA/data/GLM/Fit_Met.Rda")
evalmetGLM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/BA/data/GAM/Eval_Met.Rda")
load("AB/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("AB/BA/data/GBM/Eval_Met.Rda")
load("AB/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

load("AB/BA/data/RF/Eval_Met.Rda")
load("AB/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]




#############plotting###########################################################################

png(file=paste0("figures/PAAB_selection/testmodels/fitVSaccuracy_AB.png"),res=300,width=1961,height=2000)
par(mfrow=c(4,4))
plot(fitmetGLM_BA$Dspear[!fitmetGLM_BA$Dspear<=0],evalmetGLM_BA$Dspear[!fitmetGLM_BA$Dspear<=0],xlab="",ylab="model quality (Dspear) ",main="Bacteria",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_AR$Dspear[!fitmetGLM_AR$Dspear<=0],evalmetGLM_AR$Dspear[!fitmetGLM_AR$Dspear<=0],xlab="",ylab="",main="Archaea",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_FU$Dspear[!fitmetGLM_FU$Dspear<=0],evalmetGLM_FU$Dspear[!fitmetGLM_FU$Dspear<=0],xlab="",ylab="",main="Fungi",pch=16,cex=.2)
abline(0,1)
plot(fitmetGLM_PR$Dspear[!fitmetGLM_PR$Dspear<=0],evalmetGLM_PR$Dspear[!fitmetGLM_PR$Dspear<=0],xlab="",ylab="",main="Protista",pch=16,cex=.2)
abline(0,1)


plot(fitmetGAM_BA$Dspear[!fitmetGAM_BA$Dspear<=0],evalmetGAM_BA$Dspear[!fitmetGAM_BA$Dspear<=0],xlab="",ylab="model quality (Dspear)",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_AR$Dspear[!fitmetGAM_AR$Dspear<=0],evalmetGAM_AR$Dspear[!fitmetGAM_AR$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_FU$Dspear[!fitmetGAM_FU$Dspear<=0],evalmetGAM_FU$Dspear[!fitmetGAM_FU$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGAM_PR$Dspear[!fitmetGAM_PR$Dspear<=0],evalmetGAM_PR$Dspear[!fitmetGAM_PR$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)

plot(fitmetGBM_BA$Dspear[!fitmetGBM_BA$Dspear<=0],evalmetGBM_BA$Dspear[!fitmetGBM_BA$Dspear<=0],xlab="",ylab="model quality (Dspear)",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_AR$Dspear[!fitmetGBM_AR$Dspear<=0],evalmetGBM_AR$Dspear[!fitmetGBM_AR$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_FU$Dspear[!fitmetGBM_FU$Dspear<=0],evalmetGBM_FU$Dspear[!fitmetGBM_FU$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetGBM_PR$Dspear[!fitmetGBM_PR$Dspear<=0],evalmetGBM_PR$Dspear[!fitmetGBM_PR$Dspear<=0],xlab="",ylab="",pch=16,cex=.2)
abline(0,1)

plot(fitmetRF_BA$Dspear[!fitmetRF_BA$Dspear<=0],evalmetRF_BA$Dspear[!fitmetRF_BA$Dspear<=0],xlab="model fit(Dspear)",ylab="model quality (Dspear)",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_AR$Dspear[!fitmetRF_AR$Dspear<=0],evalmetRF_AR$Dspear[!fitmetRF_AR$Dspear<=0],xlab="model fit(Dspear)",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_FU$Dspear[!fitmetRF_FU$Dspear<=0],evalmetRF_FU$Dspear[!fitmetRF_FU$Dspear<=0],xlab="model fit(Dspear)",ylab="",pch=16,cex=.2)
abline(0,1)
plot(fitmetRF_PR$Dspear[!fitmetRF_PR$Dspear<=0],evalmetRF_PR$Dspear[!fitmetRF_PR$Dspear<=0],xlab="model fit(Dspear)",ylab="",pch=16,cex=.2)
abline(0,1)
par(mfrow=c(1,1))
dev.off()





