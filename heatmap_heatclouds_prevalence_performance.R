#create graphs and  heatmaps for  TSS vs algo vs prevalence
#
tabPR <- tabPRperc <- as.data.frame(matrix(NA,nrow=4,ncol=5))
colnames(tabPR) <- colnames(tabPRperc) <- c("Good","Useful_Val", "Bad_Val","Good_Flav", "Bad_Flav")
rownames(tabPR) <- rownames(tabPRperc) <- c("GLM", "GAM", "GBM", "RF")


# which(colnames(OTUdata)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
# annoying  Champi = 3019
load("PA/PR/data/GLM/Eval.Rda")
load("PA/PR/data/GLM/Fit_Met.Rda")
if(nrow(Eval)==3161){
evalmetGLM<-Eval[-3019,]
} else {evalmetGLM <- Eval}
if(nrow(fitmet)==3161){
fitmetGLM<-fitmet[-3019,]
} else {fitmetGLM <- fitmet}
# load(paste0("PA/PR/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/PR/data/GAM/Eval.Rda")
load("PA/PR/data/GAM/Fit_Met.Rda")
if(nrow(Eval)==3161){
  evalmetGAM<-Eval[-3019,]
} else {evalmetGAM <- Eval}
if(nrow(fitmet)==3161){
  fitmetGAM<-fitmet[-3019,]
} else {fitmetGAM <- fitmet}

load("PA/PR/data/GBM/Eval.Rda")
load("PA/PR/data/GBM/Fit_Met.Rda")
if(nrow(Eval)==3161){
  evalmetGBM<-Eval[-3019,]
} else {evalmetGBM <- Eval}
if(nrow(fitmet)==3161){
  fitmetGBM<-fitmet[-3019,]
} else {fitmetGBM <- fitmet}
load("PA/PR/data/RF/Eval.Rda")
load("PA/PR/data/RF/Fit_Met.Rda")
if(nrow(Eval)==3161){
  evalmetRF<-Eval[-3019,]
} else {evalmetRF <- Eval}
if(nrow(fitmet)==3161){
  fitmetRF<-fitmet[-3019,]
} else {fitmetRF <- fitmet}



# Result table, number of model per group
#number of ASV before modelling

length(fitmetGLM$TSS)
length(fitmetGAM$TSS)
length(fitmetGLM$TSS)
length(fitmetGLM$TSS)
#number of modelled (crashed models and prev>0.95 excluded)
sum(!(is.na(fitmetGLM$TSS)))
sum(!(is.na(fitmetGAM$TSS)))
sum(!(is.na(fitmetGBM$TSS)))
sum(!(is.na(fitmetRF$TSS)))
#Number of nodel for which the fit is >0.2 (TSS)
sum(fitmetGLM$TSS>=0.2,na.rm=TRUE)
sum(fitmetGAM$TSS>=0.2,na.rm=TRUE)
sum(fitmetGBM$TSS>=0.2,na.rm=TRUE)
sum(fitmetRF$TSS>=0.2,na.rm=TRUE)
sum(evalmetGLM$TSS>=0.2,na.rm=TRUE)
sum(evalmetGAM$TSS>=0.2,na.rm=TRUE)
sum(evalmetGBM$TSS>=0.2,na.rm=TRUE)
sum(evalmetRF$TSS>=0.2,na.rm=TRUE)


tabPR[1,] <-c(length(which(evalmetGLM$TSS>=0.6)),length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6)),length(which(evalmetGLM$TSS<=0.2)),length(which(evalmetGLM$TSS>=0.4)),length(which(evalmetGLM$TSS<0.4)))
tabPR[2,] <-c(length(which(evalmetGAM$TSS>=0.6)),length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6)),length(which(evalmetGAM$TSS<=0.2)),length(which(evalmetGAM$TSS>=0.4)),length(which(evalmetGAM$TSS<0.4)))
tabPR[3,] <-c(length(which(evalmetGBM$TSS>=0.6)),length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6)),length(which(evalmetGBM$TSS<=0.2)),length(which(evalmetGBM$TSS>=0.4)),length(which(evalmetGBM$TSS<0.4)))
tabPR[4,] <-c(length(which(evalmetRF$TSS>=0.6)),length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6)),length(which(evalmetRF$TSS<=0.2)),length(which(evalmetRF$TSS>=0.4)),length(which(evalmetRF$TSS<0.4)))
tabPRperc[1,] <-c(length(which(evalmetGLM$TSS>=0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<=0.2))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.4))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<0.4))/nrow(evalmetGLM)*100)
tabPRperc[2,] <-c(length(which(evalmetGAM$TSS>=0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<=0.2))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.4))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<0.4))/nrow(evalmetGAM)*100)
tabPRperc[3,] <-c(length(which(evalmetGBM$TSS>=0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<=0.2))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.4))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<0.4))/nrow(evalmetGBM)*100)
tabPRperc[4,] <-c(length(which(evalmetRF$TSS>=0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<=0.2))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.4))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<0.4))/nrow(evalmetRF)*100)
save(tabPRperc,file="Prop_ASV_well_modelled_protists.Rda")



#best ASV modelled prot
GoodGLM<-which(evalmetGLM$TSS>=quantile(evalmetGLM$TSS,0.975,na.rm=TRUE))#top 2.5 percent
GoodGAM<-which(evalmetGAM$TSS>=quantile(evalmetGAM$TSS,0.975,na.rm=TRUE))#top 2.5 percent
GoodGBM<-which(evalmetGBM$TSS>=quantile(evalmetGBM$TSS,0.975,na.rm=TRUE))#top 2.5 percent
GoodRF<-which(evalmetRF$TSS>=quantile(evalmetRF$TSS,0.975,na.rm=TRUE))#top 2.5 percent

GoodAll<-which(evalmetGLM$TSS>=quantile(evalmetGLM$TSS,0.975,na.rm=TRUE)&evalmetGBM$TSS>=quantile(evalmetGBM$TSS,0.975,na.rm=TRUE)&evalmetRF$TSS>=quantile(evalmetRF$TSS,0.975,na.rm=TRUE)&evalmetGAM$TSS>=quantile(evalmetGAM$TSS,0.975,na.rm=TRUE))#61

EZmodel_prot <- sample(GoodAll,10,replace=FALSE)
evalmetGLM[EZmodel_prot,]$TSS
evalmetGAM[EZmodel_prot,]$TSS
evalmetGBM[EZmodel_prot,]$TSS
evalmetRF[EZmodel_prot,]$TSS
save(EZmodel_prot,file="EZmodel_prot.Rda")

BadGLM<-which(evalmetGLM$TSS<0.2)
BadAll<-which(evalmetGLM$TSS<0.4 & evalmetGLM$TSS>0 & evalmetGAM$TSS<0.4 & evalmetGAM$TSS>0 &
                evalmetGBM$TSS<0.4 & evalmetGBM$TSS>0 & evalmetRF$TSS<0.4 & evalmetRF$TSS>0)
BadAll<-which(evalmetGLM$TSS<0.2 & evalmetGLM$TSS>0.15 & evalmetGAM$TSS<0.2 & evalmetGAM$TSS>0.15 &
                evalmetGBM$TSS<0.2 & evalmetGBM$TSS>0.15 & evalmetRF$TSS<0.2 & evalmetRF$TSS>0.15)
summary(evalmetGLM[BadAll,]$TSS)
Hardmodel_prot <- BadAll
save(Hardmodel_prot,file="Hardmodel_prot.Rda")


library(ggplot2)
library(paletteer)
#color = prevalence
load(paste0("PA/PR/data/OTUdata.Rda"))
prevalence <-apply(OTUdata,2,function(X){sum(X)/length(X)})
ncol(OTUdata)
TSStab<-data.frame(GLM=evalmetGLM$TSS,GAM=evalmetGAM$TSS,GBM=evalmetGBM$TSS,RF=evalmetRF$TSS)
GAMvsGLM_TSS_protists <- ggplot(TSStab,aes(x=GAM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GAMvsGLM_TSS_protists")
GBMvsGLM_TSS_protists <- ggplot(TSStab,aes(x=GBM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGLM_TSS_protists")
RFvsGLM_TSS_protists <- ggplot(TSStab,aes(x=RF,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGLM_TSS_protists")

GBMvsGAM_TSS_protists <- ggplot(TSStab,aes(x=GBM,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGAM_TSS_protists")
RFvsGAM_TSS_protists <- ggplot(TSStab,aes(x=RF,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGAM_TSS_protists")
RFvsGBM_TSS_protists <- ggplot(TSStab,aes(x=GBM,y=RF,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsRF_TSS_protists")
library(gridExtra)
library(grid)
PrevTestPlot<-grid.arrange(GAMvsGLM_TSS_protists,GBMvsGLM_TSS_protists,RFvsGLM_TSS_protists,GBMvsGAM_TSS_protists,RFvsGAM_TSS_protists,RFvsGBM_TSS_protists)

pdf(file="figures/testmodels/Protists_prevalence_vs_algo_vs_predictability.pdf")


plot(PrevTestPlot)

dev.off()












#####################################################################################################

tabFU <- tabFUperc <- as.data.frame(matrix(NA,nrow=4,ncol=5))
colnames(tabFU) <- colnames(tabFUperc) <- c("Good","Useful_Val", "Bad_Val","Good_Flav", "Bad_Flav")
rownames(tabFU) <- rownames(tabFUperc) <- c("GLM", "GAM", "GBM", "RF")

load("PA/FU/data/GLM/Eval.Rda")
evalmetGLM<-Eval
# load(paste0("PA/FU/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien FUevalence - model quality
load("PA/FU/data/GAM/Eval.Rda")
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
evalmetGAM<-Eval

load("PA/FU/data/GBM/Eval.Rda")
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
evalmetGBM<-Eval

load("PA/FU/data/RF/Eval.Rda")
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
evalmetRF <- Eval

tabFU[1,] <-c(length(which(evalmetGLM$TSS>=0.6)),length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6)),length(which(evalmetGLM$TSS<=0.2)),length(which(evalmetGLM$TSS>=0.4)),length(which(evalmetGLM$TSS<0.4)))
tabFU[2,] <-c(length(which(evalmetGAM$TSS>=0.6)),length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6)),length(which(evalmetGAM$TSS<=0.2)),length(which(evalmetGAM$TSS>=0.4)),length(which(evalmetGAM$TSS<0.4)))
tabFU[3,] <-c(length(which(evalmetGBM$TSS>=0.6)),length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6)),length(which(evalmetGBM$TSS<=0.2)),length(which(evalmetGBM$TSS>=0.4)),length(which(evalmetGBM$TSS<0.4)))
tabFU[4,] <-c(length(which(evalmetRF$TSS>=0.6)),length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6)),length(which(evalmetRF$TSS<=0.2)),length(which(evalmetRF$TSS>=0.4)),length(which(evalmetRF$TSS<0.4)))
tabFUperc[1,] <-c(length(which(evalmetGLM$TSS>=0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<=0.2))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.4))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<0.4))/nrow(evalmetGLM)*100)
tabFUperc[2,] <-c(length(which(evalmetGAM$TSS>=0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<=0.2))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.4))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<0.4))/nrow(evalmetGAM)*100)
tabFUperc[3,] <-c(length(which(evalmetGBM$TSS>=0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<=0.2))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.4))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<0.4))/nrow(evalmetGBM)*100)
tabFUperc[4,] <-c(length(which(evalmetRF$TSS>=0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<=0.2))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.4))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<0.4))/nrow(evalmetRF)*100)
save(tabFUperc,file="prop_ASV_well_modelled_fungi.Rda")

#best ASV modelled prot
GoodGLM<-which(evalmetGLM$TSS>=0.6)#10
GoodGLM<-which(evalmetGLM$TSS>=0.4)#617
GoodGLMGBM<-which(evalmetGLM$TSS>=0.4&evalmetGBM$TSS>=0.4)#458
GoodGLMGBMRF<-which(evalmetGLM$TSS>=0.4&evalmetGBM$TSS>=0.4&evalmetRF$TSS>=0.4)#179
GoodAll<-which(evalmetGLM$TSS>=quantile(evalmetGLM$TSS,0.975,na.rm=TRUE)&evalmetGBM$TSS>=quantile(evalmetGBM$TSS,0.975,na.rm=TRUE)&evalmetRF$TSS>=quantile(evalmetRF$TSS,0.975,na.rm=TRUE)&evalmetGAM$TSS>=quantile(evalmetGAM$TSS,0.975,na.rm=TRUE))#61
summary(evalmetGAM[GoodAll,]$TSS)
summary(evalmetGAM[GoodAll,]$TSS)
EZmodel_fung <- sample(GoodAll,10,replace=FALSE)
evalmetGLM[EZmodel_prot,]$TSS
evalmetGAM[EZmodel_prot,]$TSS
evalmetGBM[EZmodel_prot,]$TSS
evalmetRF[EZmodel_prot,]$TSS
save(EZmodel_fung,file="EZmodel_fungi.Rda")

BadGLM<-which(evalmetGLM$TSS<0.2)
BadAll<-which(evalmetGLM$TSS<0.4 & evalmetGLM$TSS>0 & evalmetGAM$TSS<0.4 & evalmetGAM$TSS>0 &
                evalmetGBM$TSS<0.4 & evalmetGBM$TSS>0 & evalmetRF$TSS<0.4 & evalmetRF$TSS>0)
BadAll<-which(evalmetGLM$TSS<0.2 & evalmetGLM$TSS>0.15 & evalmetGAM$TSS<0.2 & evalmetGAM$TSS>0.15 &
                evalmetGBM$TSS<0.2 & evalmetGBM$TSS>0.15 & evalmetRF$TSS<0.2 & evalmetRF$TSS>0.15)
summary(evalmetGLM[BadAll,]$TSS)
Hardmodel_fung <- sample(BadAll,10,replace=FALSE)
save(Hardmodel_fung,file="Hardmodel_fung.Rda")





library(ggplot2)
library(paletteer)
#color = prevalence
load(paste0("PA/FU/data/OTUdata.Rda"))
prevalence <-apply(OTUdata,2,function(X){sum(X)/length(X)})

TSStab<-data.frame(GLM=evalmetGLM$TSS,GAM=evalmetGAM$TSS,GBM=evalmetGBM$TSS,RF=evalmetRF$TSS)
GAMvsGLM_TSS_fungi <- ggplot(TSStab,aes(x=GAM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GAMvsGLM_TSS_fungi")
GBMvsGLM_TSS_fungi <- ggplot(TSStab,aes(x=GBM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGLM_TSS_fungi")
RFvsGLM_TSS_fungi <- ggplot(TSStab,aes(x=RF,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGLM_TSS_fungi")

GBMvsGAM_TSS_fungi <- ggplot(TSStab,aes(x=GBM,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGAM_TSS_fungi")
RFvsGAM_TSS_fungi <- ggplot(TSStab,aes(x=RF,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGAM_TSS_fungi")
RFvsGBM_TSS_fungi <- ggplot(TSStab,aes(x=GBM,y=RF,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsRF_TSS_fungi")
library(gridExtra)
library(grid)
PrevTestPlot<-grid.arrange(GAMvsGLM_TSS_fungi,GBMvsGLM_TSS_fungi,RFvsGLM_TSS_fungi,GBMvsGAM_TSS_fungi,RFvsGAM_TSS_fungi,RFvsGBM_TSS_fungi)
pdf(file="figures/testmodels/Fungi_prevalence_vs_algo_vs_predictability.pdf")


plot(PrevTestPlot)

dev.off()


tabBA <- tabBAperc <- as.data.frame(matrix(NA,nrow=4,ncol=5))
colnames(tabBA) <- colnames(tabBAperc) <- c("Good","UseBAl_Val", "Bad_Val","Good_Flav", "Bad_Flav")
rownames(tabBA) <- rownames(tabBAperc) <- c("GLM", "GAM", "GBM", "RF")

load("PA/BA/data/GLM/Eval.Rda")
evalmetGLM<-Eval
# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien BAevalence - model quality
load("PA/BA/data/GAM/Eval.Rda")
evalmetGAM<-Eval

load("PA/BA/data/GBM/Eval.Rda")
evalmetGBM<-Eval

load("PA/BA/data/RF/Eval.Rda")
evalmetRF <- Eval

tabBA[1,] <-c(length(which(evalmetGLM$TSS>=0.6)),length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6)),length(which(evalmetGLM$TSS<=0.2)),length(which(evalmetGLM$TSS>=0.4)),length(which(evalmetGLM$TSS<0.4)))
tabBA[2,] <-c(length(which(evalmetGAM$TSS>=0.6)),length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6)),length(which(evalmetGAM$TSS<=0.2)),length(which(evalmetGAM$TSS>=0.4)),length(which(evalmetGAM$TSS<0.4)))
tabBA[3,] <-c(length(which(evalmetGBM$TSS>=0.6)),length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6)),length(which(evalmetGBM$TSS<=0.2)),length(which(evalmetGBM$TSS>=0.4)),length(which(evalmetGBM$TSS<0.4)))
tabBA[4,] <-c(length(which(evalmetRF$TSS>=0.6)),length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6)),length(which(evalmetRF$TSS<=0.2)),length(which(evalmetRF$TSS>=0.4)),length(which(evalmetRF$TSS<0.4)))
tabBAperc[1,] <-c(length(which(evalmetGLM$TSS>=0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.2&evalmetGLM$TSS<0.6))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<=0.2))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS>=0.4))/nrow(evalmetGLM)*100,length(which(evalmetGLM$TSS<0.4))/nrow(evalmetGLM)*100)
tabBAperc[2,] <-c(length(which(evalmetGAM$TSS>=0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.2&evalmetGAM$TSS<0.6))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<=0.2))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS>=0.4))/nrow(evalmetGAM)*100,length(which(evalmetGAM$TSS<0.4))/nrow(evalmetGAM)*100)
tabBAperc[3,] <-c(length(which(evalmetGBM$TSS>=0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.2&evalmetGBM$TSS<0.6))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<=0.2))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS>=0.4))/nrow(evalmetGBM)*100,length(which(evalmetGBM$TSS<0.4))/nrow(evalmetGBM)*100)
tabBAperc[4,] <-c(length(which(evalmetRF$TSS>=0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.2&evalmetRF$TSS<0.6))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<=0.2))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS>=0.4))/nrow(evalmetRF)*100,length(which(evalmetRF$TSS<0.4))/nrow(evalmetRF)*100)
save(tabBAperc,file="prop_ASV_well_modelled_bactarchs.Rda")


#best ASV modelled prot
GoodGLM<-which(evalmetGLM$TSS>=0.6)#10
GoodGLM<-which(evalmetGLM$TSS>=0.4)#617
GoodGLMGBM<-which(evalmetGLM$TSS>=0.4&evalmetGBM$TSS>=0.4)#458
GoodGLMGBMRF<-which(evalmetGLM$TSS>=0.4&evalmetGBM$TSS>=0.4&evalmetRF$TSS>=0.4)#179
GoodAll<-which(evalmetGLM$TSS>=quantile(evalmetGLM$TSS,0.975,na.rm=TRUE)&evalmetGBM$TSS>=quantile(evalmetGBM$TSS,0.975,na.rm=TRUE)&evalmetRF$TSS>=quantile(evalmetRF$TSS,0.975,na.rm=TRUE)&evalmetGAM$TSS>=quantile(evalmetGAM$TSS,0.975,na.rm=TRUE))#61
summary(evalmetGAM[GoodAll,]$TSS)
summary(evalmetGAM[GoodAll,]$TSS)
EZmodel_bactarch <- sample(GoodAll,10,replace=FALSE)
evalmetGLM[EZmodel_prot,]$TSS
evalmetGAM[EZmodel_prot,]$TSS
evalmetGBM[EZmodel_prot,]$TSS
evalmetRF[EZmodel_prot,]$TSS
save(EZmodel_bactarch,file="EZmodel_bactarch.Rda")

BadGLM<-which(evalmetGLM$TSS<0.2)
BadAll<-which(evalmetGLM$TSS<0.4 & evalmetGLM$TSS>0 & evalmetGAM$TSS<0.4 & evalmetGAM$TSS>0 &
                evalmetGBM$TSS<0.4 & evalmetGBM$TSS>0 & evalmetRF$TSS<0.4 & evalmetRF$TSS>0)
BadAll<-which(evalmetGLM$TSS<0.2 & evalmetGLM$TSS>0.15 & evalmetGAM$TSS<0.2 & evalmetGAM$TSS>0.15 &
                evalmetGBM$TSS<0.2 & evalmetGBM$TSS>0.15 & evalmetRF$TSS<0.2 & evalmetRF$TSS>0.15)
summary(evalmetGLM[BadAll,]$TSS)
Hardmodel_bactarch <- sample(BadAll,10,replace=FALSE)
save(Hardmodel_fung,file="Hardmodel_bactarch.Rda")





library(ggplot2)
library(paletteer)
#color = prevalence
load(paste0("PA/BA/data/OTUdata.Rda"))
prevalence <-apply(OTUdata,2,function(X){sum(X)/length(X)})

TSStab<-data.frame(GLM=evalmetGLM$TSS,GAM=evalmetGAM$TSS,GBM=evalmetGBM$TSS,RF=evalmetRF$TSS)
GAMvsGLM_TSS_bactarch <- ggplot(TSStab,aes(x=GAM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GAMvsGLM_TSS_bactarch")
GBMvsGLM_TSS_bactarch <- ggplot(TSStab,aes(x=GBM,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGLM_TSS_bactarch")
RFvsGLM_TSS_bactarch <- ggplot(TSStab,aes(x=RF,y=GLM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGLM_TSS_bactarch")

GBMvsGAM_TSS_bactarch <- ggplot(TSStab,aes(x=GBM,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsGAM_TSS_bactarch")
RFvsGAM_TSS_bactarch <- ggplot(TSStab,aes(x=RF,y=GAM,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("RFvsGAM_TSS_bactarch")
RFvsGBM_TSS_bactarch <- ggplot(TSStab,aes(x=GBM,y=RF,colour=prevalence))+ 
  geom_point(size=1)+scale_color_gradient2(midpoint=0.5, low="lightgrey", mid="black",
                                           high="lightgrey") +
  ggtitle("GBMvsRF_TSS_bactarch")
library(gridExtra)
library(grid)
grid.arrange(GAMvsGLM_TSS_bactarch,GBMvsGLM_TSS_bactarch,RFvsGLM_TSS_bactarch,GBMvsGAM_TSS_bactarch,RFvsGAM_TSS_bactarch,RFvsGBM_TSS_bactarch)




#heatmap % of ASV being well modelled
# library(ggplot2)
tabPRperc$id<-rownames(tabPRperc)
ggplot(melt(tabPRperc,id.var="id"),aes(variable, id,fill=value))+geom_tile(colour = "white")+scale_fill_gradient(low="white",high="darkgrey")

tabFUperc$id<-rownames(tabFUperc)
ggplot(melt(tabFUperc,id.var="id"),aes(variable, id,fill=value))+geom_tile(colour = "white")+scale_fill_gradient(low="white",high="darkgrey")

tabBAperc$id<-rownames(tabBAperc)
ggplot(melt(tabBAperc,id.var="id"),aes(variable, id,fill=value))+geom_tile(colour = "white")+scale_fill_gradient(low="white",high="darkgrey")

