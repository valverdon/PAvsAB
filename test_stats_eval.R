library(tidyverse)
library(psych)
dataset=c(AR='#5CB800',BA='#FFD700',FU='#6B4C62',PR='#457EB0')

###data gathering###
#PR

load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.

#GLM
load("PA/PR/data/GLM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGLM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGLM_PR2<-evalmetGLM_PR[1:18]
evalmetGLM_PR2$group<-"PR"
# summary(Eval)
#GAM
load("PA/PR/data/GAM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGAM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGAM_PR2<-evalmetGAM_PR[1:18]
evalmetGAM_PR2$group<-"PR"

#GBM
load("PA/PR/data/GBM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGBM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGBM_PR2<-evalmetGBM_PR[1:18]
evalmetGBM_PR2$group<-"PR"
#RF
load("PA/PR/data/RF/Eval.Rda")

if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetRF_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetRF_PR2<-evalmetRF_PR[1:18]
evalmetRF_PR2$group<-"PR"

#FU
#GLM
load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.
load("PA/FU/data/GLM/Eval.Rda")
evalmetGLM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGLM_FU2<-evalmetGLM_FU[1:18]
evalmetGLM_FU2$group<-"FU"
#GAM
load("PA/FU/data/GAM/Eval.Rda")
evalmetGAM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGAM_FU2<-evalmetGAM_FU[1:18]
evalmetGAM_FU2$group<-"FU"
#GBM
load("PA/FU/data/GBM/Eval.Rda")
evalmetGBM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGBM_FU2<-evalmetGBM_FU[1:18]
evalmetGBM_FU2$group<-"FU"
#RF
load("PA/FU/data/RF/Eval.Rda")
evalmetRF_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetRF_FU2<-evalmetRF_FU[1:18]
evalmetRF_FU2$group<-"FU"

#BA
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum
#GLM
load("PA/BA/data/GLM/Eval.Rda")
evalmetGLM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA2<-evalmetGLM_BA[1:18]
evalmetGLM_BA2$group<-"BA"
evalmetGLM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR2<-evalmetGLM_AR[1:18]
evalmetGLM_AR2$group<-"AR"
#GAM
load("PA/BA/data/GAM/Eval.Rda")
evalmetGAM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA2<-evalmetGAM_BA[1:18]
evalmetGAM_BA2$group<-"BA"
evalmetGAM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR2<-evalmetGAM_AR[1:18]
evalmetGAM_AR2$group<-"AR"
#GBM
load("PA/BA/data/GBM/Eval.Rda")
evalmetGBM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA2<-evalmetGBM_BA[1:18]
evalmetGBM_BA2$group<-"BA"
evalmetGBM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR2<-evalmetGBM_AR[1:18]
evalmetGBM_AR2$group<-"AR"
#RF
load("PA/BA/data/RF/Eval.Rda")
evalmetRF_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA2<-evalmetRF_BA[1:18]
evalmetRF_BA2$group<-"BA"
evalmetRF_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR2<-evalmetRF_AR[1:18]
evalmetRF_AR2$group<-"AR"

evalmetGLM_all<-rbind(evalmetGLM_PR2,evalmetGLM_FU2,evalmetGLM_BA2,evalmetGLM_AR2)
obj_to_plot_GLM<-evalmetGLM_all[!(is.na(evalmetGLM_all$TSS_adj)),c("TSS_adj","group")]
evalmetGAM_all<-rbind(evalmetGAM_PR2,evalmetGAM_FU2,evalmetGAM_BA2,evalmetGAM_AR2)
obj_to_plot_GAM<-evalmetGAM_all[!(is.na(evalmetGAM_all$TSS_adj)),c("TSS_adj","group")]
evalmetGBM_all<-rbind(evalmetGBM_PR2,evalmetGBM_FU2,evalmetGBM_BA2,evalmetGBM_AR2)
obj_to_plot_GBM<-evalmetGBM_all[!(is.na(evalmetGBM_all$TSS_adj)),c("TSS_adj","group")]
evalmetRF_all<-rbind(evalmetRF_PR2,evalmetRF_FU2,evalmetRF_BA2,evalmetRF_AR2)
obj_to_plot_RF<-evalmetRF_all[!(is.na(evalmetRF_all$TSS_adj)),c("TSS_adj","group")]
sum(evalmetGLM_all$TSS_sign,na.rm=TRUE)/length(evalmetGLM_all$TSS_adj)
sum(evalmetGLM_all$TSS_adj>0.4,na.rm=TRUE)/length(evalmetGLM_all$TSS_adj)
Summary_TSSadjGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )


#taxo supplementary table
PRtaxo_supp<-PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
colnames(PRtaxo_supp)<-c("seq","level1","conf.1","level2","conf.2","level3","conf.3","level4","conf.4","level5","conf.5","level6","conf.6","level7","conf.7","level8","conf.8","level9","conf.9","Taxa")
# write.csv(PRtaxo_supp, file="figures/PAAB_selection/Figshare_tables/PRtaxo_supp.csv")
Taxo_BAARFU<-rbind(BAtaxo[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Archaea",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
unique(Taxo_BAARFU$Kingdom)
# write.csv(Taxo_BAARFU, file="figures/PAAB_selection/Figshare_tables/BAARFUtaxo_supp.csv")

#supp able model quality PA part
tobindPRGLM<-evalmetGLM_PR[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindPRGLM$marker<-"18S_protist"
tobindPRGLM<-tobindPRGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindFUGLM<-evalmetGLM_FU[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindFUGLM$marker<-"ITS_fungi"
tobindFUGLM<-tobindFUGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindBAGLM<-evalmetGLM_BA[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindBAGLM$marker<-"16S_bacteria-archaea"
tobindBAGLM<-tobindBAGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindARGLM<-evalmetGLM_AR[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindARGLM$marker<-"16S_bacteria-archaea"
tobindARGLM<-tobindARGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]

tobindPRGAM<-evalmetGAM_PR[c("auc","kap","TSS","TSS_adj")]
tobindFUGAM<-evalmetGAM_FU[c("auc","kap","TSS","TSS_adj")]
tobindBAGAM<-evalmetGAM_BA[c("auc","kap","TSS","TSS_adj")]
tobindARGAM<-evalmetGAM_AR[c("auc","kap","TSS","TSS_adj")]
tobindPRGBM<-evalmetGBM_PR[c("auc","kap","TSS","TSS_adj")]
tobindFUGBM<-evalmetGBM_FU[c("auc","kap","TSS","TSS_adj")]
tobindBAGBM<-evalmetGBM_BA[c("auc","kap","TSS","TSS_adj")]
tobindARGBM<-evalmetGBM_AR[c("auc","kap","TSS","TSS_adj")]
tobindPRRF<-evalmetRF_PR[c("auc","kap","TSS","TSS_adj")]
tobindFURF<-evalmetRF_FU[c("auc","kap","TSS","TSS_adj")]
tobindBARF<-evalmetRF_BA[c("auc","kap","TSS","TSS_adj")]
tobindARRF<-evalmetRF_AR[c("auc","kap","TSS","TSS_adj")]


evalmet_allPR<-cbind(tobindPRGLM,tobindPRGAM,tobindPRGBM,tobindPRRF)
colnames(evalmet_allPR)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allFU<-cbind(tobindFUGLM,tobindFUGAM,tobindFUGBM,tobindFURF)
colnames(evalmet_allFU)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allBA<-cbind(tobindBAGLM,tobindBAGAM,tobindBAGBM,tobindBARF)
colnames(evalmet_allBA)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allAR<-cbind(tobindARGLM,tobindARGAM,tobindARGBM,tobindARRF)
colnames(evalmet_allAR)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allall<-rbind(evalmet_allPR,evalmet_allFU,evalmet_allBA,evalmet_allAR)
# write.csv(evalmet_allall, file="figures/PAAB_selection/Figshare_tables/Model_quality_phylotype_PA.csv")             

###Statstics###
anovaGLM<-aov(TSS_adj ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM<-TukeyHSD(anovaGLM,ordered=TRUE)$group
Cohens_D_GLM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM)<-c("lower","effect","upper")
rownames(Cohens_D_GLM)<-rownames(tuk_GLM)
for (i in 1:nrow(tuk_GLM)){#i=1
  coupletolook <- word(rownames(tuk_GLM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGLM$mean[which(Summary_TSSadjGLM$group==coupletolook[1])] - Summary_TSSadjGLM$mean[which(Summary_TSSadjGLM$group==coupletolook[2])])/sd(evalmetGLM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GLM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGLM$count[which(Summary_TSSadjGLM$group==coupletolook[1])], n2=Summary_TSSadjGLM$count[which(Summary_TSSadjGLM$group==coupletolook[2])]) 
}

anovaGAM<-aov(TSS_adj ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM<-TukeyHSD(anovaGAM,ordered=TRUE)$group
Cohens_D_GAM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM)<-c("lower","effect","upper")
rownames(Cohens_D_GAM)<-rownames(tuk_GAM)
for (i in 1:nrow(tuk_GAM)){#i=1
  coupletolook <- word(rownames(tuk_GAM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGAM$mean[which(Summary_TSSadjGAM$group==coupletolook[1])] - Summary_TSSadjGAM$mean[which(Summary_TSSadjGAM$group==coupletolook[2])])/sd(evalmetGAM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GAM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGAM$count[which(Summary_TSSadjGAM$group==coupletolook[1])], n2=Summary_TSSadjGAM$count[which(Summary_TSSadjGAM$group==coupletolook[2])]) 
}

anovaGBM<-aov(TSS_adj ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM<-TukeyHSD(anovaGBM,ordered=TRUE)$group
Cohens_D_GBM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM)<-c("lower","effect","upper")
rownames(Cohens_D_GBM)<-rownames(tuk_GBM)
for (i in 1:nrow(tuk_GBM)){#i=1
  coupletolook <- word(rownames(tuk_GBM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGBM$mean[which(Summary_TSSadjGBM$group==coupletolook[1])] - Summary_TSSadjGBM$mean[which(Summary_TSSadjGBM$group==coupletolook[2])])/sd(evalmetGBM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GBM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGBM$count[which(Summary_TSSadjGBM$group==coupletolook[1])], n2=Summary_TSSadjGBM$count[which(Summary_TSSadjGBM$group==coupletolook[2])]) 
}

anovaRF<-aov(TSS_adj ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF<-TukeyHSD(anovaRF,ordered=TRUE)$group
Cohens_D_RF<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF)<-c("lower","effect","upper")
rownames(Cohens_D_RF)<-rownames(tuk_RF)
for (i in 1:nrow(tuk_RF)){#i=1
  coupletolook <- word(rownames(tuk_RF)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjRF$mean[which(Summary_TSSadjRF$group==coupletolook[1])] - Summary_TSSadjRF$mean[which(Summary_TSSadjRF$group==coupletolook[2])])/sd(evalmetRF_all$TSS_adj, na.rm = TRUE)
  Cohens_D_RF[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjRF$count[which(Summary_TSSadjRF$group==coupletolook[1])], n2=Summary_TSSadjRF$count[which(Summary_TSSadjRF$group==coupletolook[2])]) 
}


###plots###
# 
# Cohens_D[rownames(Cohens_D)=="AR-BA",2]
print("GLM")
print(Cohens_D_GLM)
print(tuk_GLM)
# c("a","ab","b","c")

stats_res_GLM<-c()
stats_res_GLM$tukey<-c("a","b","c","d")
stats_res_GLM$cohen<-c("a","ab","b","c")
stats_res_GLM<-as.data.frame(stats_res_GLM)
plotGLM<-ggplot(obj_to_plot_GLM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=20), 
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=20), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GLM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=10)+ ggtitle("GLM")


print("GAM")
print(Cohens_D_GAM)
print(tuk_GAM)
# c("a","ab","b","c")

stats_res_GAM<-c()
stats_res_GAM$tukey<-c("a","b","c","d")
stats_res_GAM$cohen<-c("a","ab","b","c")
stats_res_GAM<-as.data.frame(stats_res_GAM)
plotGAM<-ggplot(obj_to_plot_GAM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=20), 
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=20), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GAM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=10)+ ggtitle("GAM")

print("GBM")
print(Cohens_D_GBM)
print(tuk_GBM)
# c("a","ab","b","c")

stats_res_GBM<-c()
stats_res_GBM$tukey<-c("a","b","c","d")
stats_res_GBM$cohen<-c("a","ab","b","c")
stats_res_GBM<-as.data.frame(stats_res_GBM)
plotGBM<-ggplot(obj_to_plot_GBM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=20), 
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=20), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GBM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=10)+ ggtitle("GBM")

print("RF")
print(Cohens_D_RF)
print(tuk_RF)
# c("a","ab","b","c")

stats_res_RF<-c()
stats_res_RF$tukey<-c("a","b","c","d")
stats_res_RF$cohen<-c("a","ab","b","c")
stats_res_RF<-as.data.frame(stats_res_RF)
plotRF<-ggplot(obj_to_plot_RF, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=20), 
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=20), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_RF, aes(x=1:4,y=0.95,label=cohen),col="black",cex=10)+ ggtitle("RF")


pdf(file="figures/distribTSSadj/distribsTSS.pdf")
plotGLM
plotGAM
plotGBM
plotRF
dev.off()

png(file=paste0("figures/distribTSSadj/distribsTSS_GLM.png"),res=300,width=1961,height=1500)
plotGLM
dev.off()
png(file=paste0("figures/distribTSSadj/distribsTSS_GAM.png"),res=300,width=1961,height=1500)
plotGAM
dev.off()
png(file=paste0("figures/distribTSSadj/distribsTSS_GBM.png"),res=300,width=1961,height=1500)
plotGBM
dev.off()
png(file=paste0("figures/distribTSSadj/distribsTSS_RF.png"),res=300,width=1961,height=1500)
plotRF
dev.off()



#BA166
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum
#GLM
load("PA/BA_166/data/GLM/Eval.Rda")
evalmetGLM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA1662<-evalmetGLM_BA166[1:19]
evalmetGLM_BA1662$group<-"BA"
evalmetGLM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR1662<-evalmetGLM_AR166[1:19]
evalmetGLM_AR1662$group<-"AR"
#GAM
load("PA/BA_166/data/GAM/Eval.Rda")
evalmetGAM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA1662<-evalmetGAM_BA166[1:19]
evalmetGAM_BA1662$group<-"BA"
evalmetGAM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR1662<-evalmetGAM_AR166[1:19]
evalmetGAM_AR1662$group<-"AR"
#GBM
load("PA/BA_166/data/GBM/Eval.Rda")
evalmetGBM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA1662<-evalmetGBM_BA166[1:19]
evalmetGBM_BA1662$group<-"BA"
evalmetGBM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR1662<-evalmetGBM_AR166[1:19]
evalmetGBM_AR1662$group<-"AR"
#RF
load("PA/BA_166/data/RF/Eval.Rda")
evalmetRF_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA1662<-evalmetRF_BA166[1:19]
evalmetRF_BA1662$group<-"BA"
evalmetRF_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR1662<-evalmetRF_AR166[1:19]
evalmetRF_AR1662$group<-"AR"

boxplot(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj)
boxplot(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj)
t.test(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj,paired=TRUE)
t.test(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj,paired=TRUE)
cohenD<-(mean(evalmetGLM_BA$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_BA166$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_BA$TSS_adj))), n2=sum(!(is.na(evalmetGLM_BA166$TSS_adj))))
Cohens_D
cohenD<-(mean(evalmetGLM_AR$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_AR166$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_AR$TSS_adj))), n2=sum(!(is.na(evalmetGLM_AR166$TSS_adj))))
Cohens_D


t.test(evalmetGLM_BA166$TSS_adj,evalmetGLM_PR$TSS_adj)
t.test(evalmetGLM_AR166$TSS_adj,evalmetGLM_PR$TSS_adj)
cohenD<-(mean(evalmetGLM_BA166$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_PR$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_BA166$TSS_adj,evalmetGLM_PR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_BA166$TSS_adj))), n2=sum(!(is.na(evalmetGLM_PR$TSS_adj))))
Cohens_D
cohenD<-(mean(evalmetGLM_AR166$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_PR$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_AR166$TSS_adj,evalmetGLM_PR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_AR166$TSS_adj))), n2=sum(!(is.na(evalmetGLM_PR$TSS_adj))))
Cohens_D













###################################################
#rel. ab


#GLM
load("AB/PR/data/GLM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGLM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGLM_PR2<-evalmetGLM_PR[1:16]
evalmetGLM_PR2$group<-"PR"
# summary(Eval)
#GAM
load("AB/PR/data/GAM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGAM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGAM_PR2<-evalmetGAM_PR[1:16]
evalmetGAM_PR2$group<-"PR"

#GBM
load("AB/PR/data/GBM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGBM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGBM_PR2<-evalmetGBM_PR[1:16]
evalmetGBM_PR2$group<-"PR"


#RF
load("AB/PR/data/RF/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetRF_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetRF_PR2<-evalmetRF_PR[1:16]
evalmetRF_PR2$group<-"PR"


#FU
#GLM
load("AB/FU/data/GLM/Eval_Met.Rda")
evalmetGLM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGLM_FU2<-evalmetGLM_FU[1:16]
evalmetGLM_FU2$group<-"FU"
#GAM
load("AB/FU/data/GAM/Eval_Met.Rda")
evalmetGAM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGAM_FU2<-evalmetGAM_FU[1:16]
evalmetGAM_FU2$group<-"FU"
#GBM
load("AB/FU/data/GBM/Eval_Met.Rda")
evalmetGBM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGBM_FU2<-evalmetGBM_FU[1:16]
evalmetGBM_FU2$group<-"FU"
#RF
load("AB/FU/data/RF/Eval_Met.Rda")
evalmetRF_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetRF_FU2<-evalmetRF_FU[1:16]
evalmetRF_FU2$group<-"FU"

#BA
#GLM
load("AB/BA/data/GLM/Eval_Met.Rda")
evalmetGLM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA2<-evalmetGLM_BA[1:16]
evalmetGLM_BA2$group<-"BA"
evalmetGLM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR2<-evalmetGLM_AR[1:16]
evalmetGLM_AR2$group<-"AR"
#GAM
load("AB/BA/data/GAM/Eval_Met.Rda")
evalmetGAM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA2<-evalmetGAM_BA[1:16]
evalmetGAM_BA2$group<-"BA"
evalmetGAM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR2<-evalmetGAM_AR[1:16]
evalmetGAM_AR2$group<-"AR"
#GBM
load("AB/BA/data/GBM/Eval_Met.Rda")
evalmetGBM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA2<-evalmetGBM_BA[1:16]
evalmetGBM_BA2$group<-"BA"
evalmetGBM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR2<-evalmetGBM_AR[1:16]
evalmetGBM_AR2$group<-"AR"
#RF
load("AB/BA/data/RF/Eval_Met.Rda")
evalmetRF_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA2<-evalmetRF_BA[1:16]
evalmetRF_BA2$group<-"BA"
evalmetRF_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR2<-evalmetRF_AR[1:16]
evalmetRF_AR2$group<-"AR"

evalmetGLM_all<-rbind(evalmetGLM_PR2,evalmetGLM_FU2,evalmetGLM_BA2,evalmetGLM_AR2)
obj_to_plot_GLM<-evalmetGLM_all[!(is.na(evalmetGLM_all$TSS_adj)),c("Dspear","group")]
evalmetGAM_all<-rbind(evalmetGAM_PR2,evalmetGAM_FU2,evalmetGAM_BA2,evalmetGAM_AR2)
obj_to_plot_GAM<-evalmetGAM_all[!(is.na(evalmetGAM_all$TSS_adj)),c("Dspear","group")]
evalmetGBM_all<-rbind(evalmetGBM_PR2,evalmetGBM_FU2,evalmetGBM_BA2,evalmetGBM_AR2)
obj_to_plot_GBM<-evalmetGBM_all[!(is.na(evalmetGBM_all$TSS_adj)),c("Dspear","group")]
evalmetRF_all<-rbind(evalmetRF_PR2,evalmetRF_FU2,evalmetRF_BA2,evalmetRF_AR2)
obj_to_plot_RF<-evalmetRF_all[!(is.na(evalmetRF_all$TSS_adj)),c("Dspear","group")]
sum(evalmetGLM_all$Dspear,na.rm=TRUE)/length(evalmetGLM_all$Dspear)
sum(evalmetGLM_all$Dspear>0.4,na.rm=TRUE)/length(evalmetGLM_all$Dspear)
sum(evalmetGLM_all$Dspear>0.2,na.rm=TRUE)/length(evalmetGLM_all$Dspear)


#tentative correction of CVRMSE
# summary(evalmetGLM_BA$RMSE)
# load(paste0("AB/PR/data/OTUdata.Rda"))
# evalmetGLM_PR$CVRMSE<-NA
# for (i in 1:length(evalmetGLM_PR$RMSE)){
#   evalmetGLM_PR$CVRMSE[i]<-evalmetGLM_PR$RMSE[i]/mean(OTUdata[,i])
# }
# summary(evalmetGLM_PR$CVRMSE)
# sapply(evalmetGLM_PR$RMSE)
# #supp able model quality AB part
# summary(evalmetGLM_PR$RMSEs)
tobindPRGLM<-evalmetGLM_PR[c("Seq","Phylum","Dspear","RMSE")]
tobindPRGLM$marker<-"18S_protist"
tobindPRGLM<-tobindPRGLM[c("Seq","marker","Phylum","Dspear","RMSE")]
tobindFUGLM<-evalmetGLM_FU[c("Seq","Phylum","Dspear","RMSE")]
tobindFUGLM$marker<-"ITS_fungi"
tobindFUGLM<-tobindFUGLM[c("Seq","marker","Phylum","Dspear","RMSE")]
tobindBAGLM<-evalmetGLM_BA[c("Seq","Phylum","Dspear","RMSE")]
tobindBAGLM$marker<-"16S_bacteria-archaea"
tobindBAGLM<-tobindBAGLM[c("Seq","marker","Phylum","Dspear","RMSE")]
tobindARGLM<-evalmetGLM_AR[c("Seq","Phylum","Dspear","RMSE")]
tobindARGLM$marker<-"16S_bacteria-archaea"
tobindARGLM<-tobindARGLM[c("Seq","marker","Phylum","Dspear","RMSE")]

tobindPRGAM<-evalmetGAM_PR[c("Dspear","RMSE")]
tobindFUGAM<-evalmetGAM_FU[c("Dspear","RMSE")]
tobindBAGAM<-evalmetGAM_BA[c("Dspear","RMSE")]
tobindARGAM<-evalmetGAM_AR[c("Dspear","RMSE")]
tobindPRGBM<-evalmetGBM_PR[c("Dspear","RMSE")]
tobindFUGBM<-evalmetGBM_FU[c("Dspear","RMSE")]
tobindBAGBM<-evalmetGBM_BA[c("Dspear","RMSE")]
tobindARGBM<-evalmetGBM_AR[c("Dspear","RMSE")]
tobindPRRF<-evalmetRF_PR[c("Dspear","RMSE")]
tobindFURF<-evalmetRF_FU[c("Dspear","RMSE")]
tobindBARF<-evalmetRF_BA[c("Dspear","RMSE")]
tobindARRF<-evalmetRF_AR[c("Dspear","RMSE")]


evalmet_allPR<-cbind(tobindPRGLM,tobindPRGAM,tobindPRGBM,tobindPRRF)
colnames(evalmet_allPR)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-RMSE","GAM-Dspear","GAM-RMSE","GBM-Dspear","GBM-RMSE","RF-Dspear","RF-RMSE")
evalmet_allFU<-cbind(tobindFUGLM,tobindFUGAM,tobindFUGBM,tobindFURF)
colnames(evalmet_allFU)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-RMSE","GAM-Dspear","GAM-RMSE","GBM-Dspear","GBM-RMSE","RF-Dspear","RF-RMSE")
evalmet_allBA<-cbind(tobindBAGLM,tobindBAGAM,tobindBAGBM,tobindBARF)
colnames(evalmet_allBA)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-RMSE","GAM-Dspear","GAM-RMSE","GBM-Dspear","GBM-RMSE","RF-Dspear","RF-RMSE")
evalmet_allAR<-cbind(tobindARGLM,tobindARGAM,tobindARGBM,tobindARRF)
colnames(evalmet_allAR)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-RMSE","GAM-Dspear","GAM-RMSE","GBM-Dspear","GBM-RMSE","RF-Dspear","RF-RMSE")
evalmet_allall<-rbind(evalmet_allPR,evalmet_allFU,evalmet_allBA,evalmet_allAR)
# write.csv(evalmet_allall, file="figures/PAAB_selection/Figshare_tables/Model_quality_phylotype_AB.csv")             


