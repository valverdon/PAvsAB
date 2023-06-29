#creates violin plots of distribution of TSS values across kingdoms for 3 algos

library(tidyverse)
library(gridExtra)
library(reshape2)
library(stringi)
PAAB = "PA"
GtoM="BA"
algo="GLM"

PrepareToPlot <- function(datatab,taxa_level,PAAB){ #datatab=ggplotdata;taxa_level="Phylum";PAAB="PA"

  

  if(taxa_level=="Phylum"){
    if (PAAB=="PA"){
      datatab_nona <- datatab[!(is.na(datatab$TSS)),] 
      datatab_melt<-reshape2::melt(datatab_nona[c("kap","TSS","auc")])
      datatab2<-cbind(datatab_melt,datatab_nona)
      levels(datatab2$variable) <- c("Kappa","TSS","AUC")
    }
    if (PAAB=="AB"){
      datatab_nona <- datatab[!(is.na(datatab$Dpear)),] 
      # test<-reshape2::melt(FitEvalP[c("id","R2","D2","MAE","MAEs","MAE_scaled","RMSE","RMSEs", "RMSE_scaled","Dspear","Dpear","Pdispersion")])
      #
      datatab_melt<-reshape2::melt(datatab_nona[c("R2","MAE","RMSE","Dspear","Dpear")])
      datatab2<-cbind(datatab_melt,datatab_nona)
      levels(datatab2$variable) <- c("R2","MAE","RMSE","Dspear","Dpear")
      }
    datatab2$value <- sapply(datatab2$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
    return(datatab2)
  }
  if (taxa_level=="Kingdom"){
    if (PAAB=="PA"){
      datatab_nona <- datatab[!(is.na(datatab$TSS)),] 
      datatab_melt<-reshape2::melt(datatab_nona[c("kap","TSS","auc")])
      datatab2<-cbind(datatab_melt,datatab_nona)
      levels(datatab2$variable) <- c("Kappa","TSS","AUC")

    }
    if (PAAB=="AB"){
      datatab_nona <- datatab[!(is.na(datatab$Dpear)),] 
      datatab_melt<-reshape2::melt(datatab_nona[c("R2","MAE","RMSE","Dspear","Dpear")])
      datatab2<-cbind(datatab_melt,datatab_nona)
      levels(datatab2$variable) <- c("R2","MAE","RMSE","Dspear","Dpear")
    }
    datatab2$value <- sapply(datatab2$value,function(X){ifelse(is.na(X),NA,ifelse(X<0,0,X))}) #put negative fits to 0
    return(datatab2)
  }
  
  }

#figures poster ISME








############################# NEEED 166 BA !!!###################################################



#Figure 1 : PA
# for (GtoM in c("PR","BA","FU")){ 
Eval_met_mat_all <-c()
Taxoall<-c()
#Taxobact = 7  TaxoPR = 12  #TaxoFU =7
for (algo in c("GLM","GAM","RF","GBM")){

for (GtoM in c("PR","BA","FU","BA_166")){ #GtoM="BA_166"
  print(GtoM)
  load(paste0("../../ASV_data/ASV_taxo/",substr(GtoM,1,2),"taxo.Rda"))

    print(algo)
    load(paste0("PA/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    evalmet$GtoM<-GtoM
    evalmet$algo<-algo
    Eval_met_mat_all<-rbind(Eval_met_mat_all,evalmet[c("auc","TSS","kap","OTU","GtoM","algo")])
    if(GtoM=="BA_166"){
    tobind<-BAtaxo[c("Kingdom", "Phylum" , "Class"  , "Order"  , "Family"  ,"Genus"  , "Species")]
    tobind$Kingdom[which(!(is.na(tobind$Kingdom)))]<-paste(tobind$Kingdom[which(!(is.na(tobind$Kingdom)))],"_166",sep="")
    Taxoall<-rbind(Taxoall,tobind)
    }else
    Taxoall<-rbind(Taxoall,eval(parse(text=paste0(GtoM,"taxo")))[c("Kingdom", "Phylum" , "Class"  , "Order"  , "Family"  ,"Genus"  , "Species")])
}
}


ggplotdata<-cbind(Eval_met_mat_all,Taxoall)


#Paired stats tests between algos and between metrics
# cortest_aucTSS<-cor.test(Eval_met_mat_all$auc,Eval_met_mat_all$TSS)
# cortest_aucTSS$estimate^2
# cortest_kapTSS<-cor.test(Eval_met_mat_all$kap,Eval_met_mat_all$TSS)
# cortest_kapTSS$estimate^2
# cortest_auckap<-cor.test(Eval_met_mat_all$auc,Eval_met_mat_all$kap)
# cortest_auckap$estimate^2
# hist(Eval_met_mat_all$auc)
# hist(Eval_met_mat_all$TSS)
# hist(Eval_met_mat_all$kap)
# Eval_met_mat_all[Eval_met_mat_all$algo=="GLM",]
# Eval_met_mat_all[Eval_met_mat_all$algo=="GAM",]
# cortest_GLMGAM<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GLM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="GAM",]$TSS)
# cortest_GLMGAM$estimate^2
# cortest_GLMGBM<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GLM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="GBM",]$TSS)
# cortest_GLMGBM$estimate^2
# cortest_GLMRF<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GLM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="RF",]$TSS)
# cortest_GLMRF$estimate^2
# cortest_GBMGAM<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GBM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="GAM",]$TSS)
# cortest_GBMGAM$estimate^2
# cortest_GBMRF<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GBM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="RF",]$TSS)
# cortest_GBMRF$estimate^2
# cortest_GAMRF<-cor.test(Eval_met_mat_all[Eval_met_mat_all$algo=="GAM",]$TSS,Eval_met_mat_all[Eval_met_mat_all$algo=="RF",]$TSS)
# cortest_GAMRF$estimate^2
cortable<-matrix(0,nrow=12,ncol=12)
rownames(cortable)<-colnames(cortable)<-c("aucGLM","aucGAM","aucRF","aucGBM","TSSGLM","TSSGAM","TSSRF","TSSGBM","kapGLM","kapGAM","kapRF","kapGBM")
for (i in 1:ncol(cortable)){#i=1
  for (j in 1:ncol(cortable)){#j=1
    data1<-Eval_met_mat_all[Eval_met_mat_all$algo==substr(rownames(cortable)[i],4,6),which(colnames(Eval_met_mat_all)==substr(rownames(cortable)[i],1,3))]
    data2<-Eval_met_mat_all[Eval_met_mat_all$algo==substr(rownames(cortable)[j],4,6),which(colnames(Eval_met_mat_all)==substr(rownames(cortable)[j],1,3))]
    cortable[i,j]<-round(cor(data1,data2,use="complete.obs"),2)
    }
}
write.csv(cortable,file="figures/cortable_algo_metric.csv")

labs <- c(1:6)#i=1
for (i in 1:length(labs)){
  test<-ggplotdata[ggplotdata$algo=="GLM" & ggplotdata$Kingdom==na.omit(unique(ggplotdata$Kingdom))[i],]
  labs[i]<-nrow(test[!(is.na(test$TSS)),])
}
labs_m<-data.frame(Kingdom=na.omit(unique(ggplotdata$Kingdom)),labs=labs)

# for (algo in c("GLM","GAM","RF")){
      FitEvalP_Am<-PrepareToPlot(datatab=ggplotdata,taxa_level="Kingdom",PAAB="PA")
      p1 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("Kappa") + scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
        ylim(0,1)
      p2 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("TSS") + scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,-.5,0.2),"cm")) + ylim(0,1)

      p3 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("AUC") + scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,0,0.2),"cm")) + ylim(0.5,1)

      p4 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,-.5,0),"cm")) + ylim(0,1)

      p5 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,-.5,0),"cm")) + ylim(0,1)

      p6 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,0,0),"cm")) + ylim(0.5,1)

      p7 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

      p8 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

      p9 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,0,0),"cm")) + ylim(0.5,1)
      
      p10 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)
      
      p11 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",6)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)
      
      p12 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.5,0,0),"cm")) + ylim(0.5,1)
      
      
      # par(mfrow=c(4,3))
      # pdf(paste0("figures/poster/KappaTSS_",algo,".pdf"))
      library(grid)
      tGLM <- textGrob("GLM",hjust=-.2)
      tGAM <- textGrob("GAM",hjust=-.1)
      tRF <- textGrob(" RF",hjust=-.2)
      tGBM <- textGrob("GBM",hjust=-.1)
      pdf(paste0("figures/poster/figure1.pdf"))
      grid.arrange(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,tGLM,tGAM,tRF,tGBM), 
                   layout_matrix = rbind(c(13,14,15,16),
                                         c(1,4,7,10),
                                         c(2,5,8,11),
                                         c(3,6,9,12)),
                  heights=c(0.1,0.5,0.5,0.5))
      dev.off()
      # plot(p2,main="GLM")
      # 
      # 
      # plot(p1)

     # Tests stats 
      TSSGLM<-FitEvalP_Am[FitEvalP_Am$algo=="GLM"&FitEvalP_Am$variable=="TSS",]
      library(dplyr)
      Summary_TSSGLM<-group_by(TSSGLM, Kingdom) %>%
        summarise(
          count = n(),
          mean = mean(TSS, na.rm = TRUE),
          sd = sd(TSS, na.rm = TRUE)
        )
      TSSGLM.aov <- aov(TSS ~ Kingdom, data = TSSGLM)
summary(TSSGLM.aov)
TukeyHSD(TSSGLM.aov,ordered=TRUE)

TSSGLM_SD<-sd(TSSGLM$TSS)
couples_kingdom_diff<-TukeyHSD(TSSGLM.aov,ordered=TRUE)$Kingdom
print(TSSGLM_SD);print(Summary_TSSGLM)
library(stringr)

#cohen's d + pearson cor
Cohens_D<-c()
PearCor<-c()
for (i in 1:nrow(couples_kingdom_diff)){#i=1
  coupletolook <- word(rownames(couples_kingdom_diff)[i],1:2,sep = "-")
  Cohens_D[i] <-(Summary_TSSGLM$mean[which(Summary_TSSGLM$Kingdom==coupletolook[1])] - Summary_TSSGLM$mean[which(Summary_TSSGLM$Kingdom==coupletolook[2])])/TSSGLM_SD
}
#pearson cor

Stats_GLM_TSS<-cbind(couples_kingdom_diff,Cohens_D,PearCor)


met_algo.aov<-list()
met_algo_SD<-list()
meanvals<-list()

Stats_val2<-list()
for (metric in c("Kappa","TSS","AUC")){#metric="Kappa"
  Stats_vals<-list("GLM"=list(),"GAM"=list(),"RF"=list(),"GBM"=list())
  for(algo in c("GLM","GAM","RF","GBM")){#algo="GLM"
    metrics_vals<-FitEvalP_Am[FitEvalP_Am$algo==algo & FitEvalP_Am$variable==metric,]
    Summary_vals<-group_by(metrics_vals, Kingdom) %>%
      summarise(
        count = n(),
        meanTSS = mean(TSS, na.rm = TRUE),
        sdTSS = sd(TSS, na.rm = TRUE),
        meanKappa = mean(kap, na.rm = TRUE),
        sdKappa = sd(kap, na.rm = TRUE),
        meanAUC = mean(auc, na.rm = TRUE),
        sdAUC = sd(auc, na.rm = TRUE)
      )

    met_algo.aov[["Kappa"]] <- aov(kap ~ Kingdom, data = metrics_vals)
    met_algo.aov[["TSS"]] <- aov(TSS ~ Kingdom, data = metrics_vals)
    met_algo.aov[["AUC"]]  <- aov(auc ~ Kingdom, data = metrics_vals)
    # summary(met_algo.aov)[2]
    met_algo_SD[["Kappa"]]<-sd(metrics_vals$kap)
    met_algo_SD[["TSS"]]<-sd(metrics_vals$TSS)
    met_algo_SD[["AUC"]]<-sd(metrics_vals$auc)
    couples_kingdom_diff<-TukeyHSD(met_algo.aov[[metric]],ordered=TRUE)$Kingdom
    
    
    #cohen's d + pearson cor
    Cohens_D<-c()
    for (i in 1:nrow(couples_kingdom_diff)){#i=1
      coupletolook <- word(rownames(couples_kingdom_diff)[i],1:2,sep = "-")
      meanvals[["Kappa"]]<-c(Summary_vals$meanKappa[which(Summary_vals$Kingdom==coupletolook[1])],Summary_vals$meanKappa[which(Summary_vals$Kingdom==coupletolook[2])])
      meanvals[["TSS"]]<-c(Summary_vals$meanTSS[which(Summary_vals$Kingdom==coupletolook[1])],Summary_vals$meanTSS[which(Summary_vals$Kingdom==coupletolook[2])])
      meanvals[["AUC"]]<-c(Summary_vals$meanAUC[which(Summary_vals$Kingdom==coupletolook[1])],Summary_vals$meanAUC[which(Summary_vals$Kingdom==coupletolook[2])])
      Cohens_D[i] <-(meanvals[[metric]][1] - meanvals[[metric]][2])/met_algo_SD[[metric]]
    }
    Stats_vals[[algo]]<-cbind(couples_kingdom_diff,Cohens_D)
  }
  Stats_val2[[metric]]<-Stats_vals
}


Signif_aov_tukey_cohens<-Stats_val2
save(Signif_aov_tukey_cohens,file="figures/poster/figure1_signif.Rda")



#Figure 2 : AB
      Eval_met_mat_all <-c()
      Taxoall<-c()
      #Taxobact = 7  TaxoPR = 12  #TaxoFU =7
      for (algo in c("GLM","GAM","RF","GBM")){
        #algo="GLM"
        for (GtoM in c("PR","BA","FU")){ 
          # GtoM="PR"
          print(GtoM)
          load(paste0("../../ASV_data/ASV_taxo/",substr(GtoM,1,2),"taxo.Rda"))
          load(paste0("AB/",GtoM,"/data/dataTotSeqSum.Rda"))

          load(paste0("AB/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
          evalmet$GtoM<-GtoM
          evalmet$algo<-algo
          taxo<-eval(parse(text=paste0(GtoM,"taxo")))[c("Kingdom", "Phylum" , "Class"  , "Order"  , "Family"  ,"Genus"  , "Species")]
          if(nrow(evalmet)!=nrow(taxo)&algo!="GLM"){
          evalmet<-evalmet[complete.cases(evalmet),]
          }
          if(nrow(evalmet)==3160&GtoM=="PR"){
            taxo<-taxo[-3019,]#champi
          }
          print(nrow(evalmet))
          print(nrow(taxo))
          

          Eval_met_mat_all<-rbind(Eval_met_mat_all,evalmet[c("R2","MAE","RMSE","Dspear","Dpear","OTU","GtoM","algo")])
          Taxoall<-rbind(Taxoall,taxo)
        }
      }
      nrow(Eval_met_mat_all)
      nrow(Taxoall)
      ##########################
      
      ggplotdata<-cbind(Eval_met_mat_all,Taxoall)
      #number of ASV per group
      labs <- c(1:4)#i=1
      for (i in 1:4){
        test<-ggplotdata[ggplotdata$algo=="GLM" & ggplotdata$Kingdom==na.omit(unique(ggplotdata$Kingdom))[i],]
        labs[i]<-nrow(test[!(is.na(test$R2)),])
      }
      labs_m<-data.frame(Kingdom=na.omit(unique(ggplotdata$Kingdom)),labs=labs)
      
      # for (algo in c("GLM","GAM","RF")){
      FitEvalP_Am<-PrepareToPlot(datatab=ggplotdata,taxa_level="Kingdom",PAAB="AB")

      p1 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("") + scale_x_discrete(labels=rep("",4)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,0,0.2),"cm")) + 
        ylim(0,1)
      p2 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("") + scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0,0,0.2),"cm"))

      p3 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",4)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      
      p4 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+  scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      
      p5 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",4)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      
      p6 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      p7 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=rep("",4)) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      
      p8 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Kingdom))&FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Kingdom, y=value))+ geom_violin(scale="width") +
        ylab("")+ scale_x_discrete(labels=c("Eukaryota"="Protista")) +
        geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=0.5), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)
      
      
      
      # par(mfrow=c(4,3))

      library(grid)
      tGLM <- textGrob("GLM")
      tGAM <- textGrob("GAM")
      tRF <- textGrob(" RF")
      tGBM <- textGrob("GBM")
      tR2 <- textGrob(" R2")
      tDpear <- textGrob(" DPear",hjust=-.2)
      empty <- textGrob(" ",hjust=-.2)
      
      pdf(paste0("figures/poster/Figure2.pdf"))
      grid.arrange(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,tGLM,tGAM,tRF,tGBM,tR2,tDpear,empty), 
                   layout_matrix = rbind(c(15,13,14),
                                        c(9,1,2),
                                         c(10,3,4),
                                         c(11,5,6),
                                        c(12,7,8)),
                    heights=c(0.1,0.5,0.5,0.5,0.5),
                   widths=c(0.1,1,1))

      dev.off()





#figure 3
PAAB="PA"
Eval_met_mat_all <-c()
Taxoall<-c()
#Taxobact = 7  TaxoPR = 12  #TaxoFU =7
for (algo in c("GLM","GAM","RF","GBM")){
  
  for (GtoM in c("PR","BA","FU","BA_166")){ 
    print(GtoM)
    load(paste0("../../ASV_data/ASV_taxo/",substr(GtoM,1,2),"taxo.Rda"))
    
    print(algo)
    load(paste0("PA/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    evalmet$GtoM<-GtoM
    evalmet$algo<-algo
    Eval_met_mat_all<-rbind(Eval_met_mat_all,evalmet[c("auc","TSS","kap","OTU","GtoM","algo")])
    Taxoall<-rbind(Taxoall,eval(parse(text=paste0(substr(GtoM,1,2),"taxo")))[c("Kingdom", "Phylum" , "Class"  , "Order"  , "Family"  ,"Genus"  , "Species")])
  }
}


##########################
ggplotdata<-cbind(Eval_met_mat_all,Taxoall)

#number of ASV per group

labs <- c(1:length(na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Archaea",]$Phylum))))#i=1
for (i in 1:length(labs)){#i=1
  test<-ggplotdata[which(ggplotdata$algo=="GLM" & ggplotdata$GtoM==GtoM & ggplotdata$Kingdom=="Archaea" & ggplotdata$Phylum==na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Archaea",]$Phylum))[i]),]
  labs[i]<-nrow(test[!(is.na(test$TSS)),])
}
labs_m<-data.frame(Phylum=unique(ggplotdata[which(ggplotdata$Kingdom=="Archaea"),]$Phylum),labs=labs)
sum(BAtaxo$Phylum=="Asgardaeota",na.rm=TRUE)
BAtaxo[which(BAtaxo$Phylum=="Asgardaeota"),]

GtoM="BA"
labs <- c(1:length(na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Bacteria",]$Phylum))))#i=1
for (i in 1:length(labs)){#i=1
  test<-ggplotdata[which(ggplotdata$algo=="GLM" & ggplotdata$GtoM==GtoM & ggplotdata$Kingdom=="Bacteria" & ggplotdata$Phylum==unique(ggplotdata[which(ggplotdata$Kingdom=="Bacteria"),]$Phylum)[i]),]
  labs[i]<-nrow(test[!(is.na(test$TSS)),])
}
labs_m<-data.frame(Phylum=unique(ggplotdata[which(ggplotdata$Kingdom=="Bacteria"),]$Phylum),labs=labs)


# for (algo in c("GLM","GAM","RF")){
FitEvalP_Am<-PrepareToPlot(datatab=ggplotdata,taxa_level="Phylum",PAAB="PA")

p1 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GLM" & FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("Kappa") + scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), plot.margin = unit(c(0,0,-.5,0.2),"cm")) + 
  ylim(0,1)
p2 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("TSS") + scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), plot.margin = unit(c(0,0,-.5,0.2),"cm")) + ylim(0,1)

p3 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("AUC") + scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0,0,0.2),"cm")) + ylim(0.5,1)

p4 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0,-.5,0),"cm")) + ylim(0,1)

p5 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0,-.5,0),"cm")) + ylim(0,1)

p6 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0,0,0),"cm")) + ylim(0.5,1)

p7 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

p8 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

p9 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0.5,0,0),"cm")) + ylim(0.5,1)

p10 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("Kappa"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

p11 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("TSS"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.5,-.5,0),"cm")) + ylim(0,1)

p12 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" & FitEvalP_Am$algo=="GBM" &FitEvalP_Am$variable%in%c("AUC"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0.5,0,0),"cm")) + ylim(0.5,1)


# par(mfrow=c(4,3))

library(grid)
tGLM <- textGrob("GLM",hjust=-.2)
tGAM <- textGrob("GAM",hjust=-.1)
tRF <- textGrob(" RF",hjust=-.2)
tGBM <- textGrob("GBM",hjust=-.2)
pdf(paste0("figures/poster/figure2_Archaea.pdf"))
grid.arrange(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,tGLM,tGAM,tRF,tGBM), 
             layout_matrix = rbind(c(13,14,15,16),
                                   c(1,4,7,10),
                                   c(2,5,8,11),
                                   c(3,6,9,12)),
             heights=c(0.1,0.5,0.5,0.5))
dev.off()







#Figure 4 : AB
#RF cassé
# load(paste0("AB/PR/data/RF/Eval_Met.Rda"))
# evalmetRF<-evalmet
# load(paste0("AB/PR/data/GLM/Eval_Met.Rda"))
# evalmetGLM<-evalmet
PAAB
Eval_met_mat_all <-c()
Taxoall<-c()
#Taxobact = 7  TaxoPR = 12  #TaxoFU =7
for (algo in c("GLM","GAM","GBM")){#algo="RF"
  
  for (GtoM in c("PR","BA","FU")){ #GtoM="PR"
    print(GtoM)
    load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
    load(paste0("AB/",GtoM,"/data/dataTotSeqSum.Rda"))
    #length(TotSeqSum)
    load(paste0("AB/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    if(nrow(eval(parse(text=paste0(GtoM,"taxo")))) == 3161 ){
      # which(colnames(OTUdata_PR_AB)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
      #annoying  Champi = 3019
      PRtaxo <- PRtaxo[-3019,]
    }
    if(nrow(evalmet) == 3161 ){
      # which(colnames(OTUdata_PR_AB)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
      #annoying  Champi = 3019
      evalmet <- evalmet[-3019,]
    }
    evalmet$GtoM<-GtoM
    evalmet$algo<-algo
    print(nrow(eval(parse(text=paste0(GtoM,"taxo"))))==nrow(evalmet))
    Eval_met_mat_all<-rbind(Eval_met_mat_all,evalmet[c("R2","MAE","RMSE","Dspear","Dpear","OTU","GtoM","algo")])
    Taxoall<-rbind(Taxoall,eval(parse(text=paste0(GtoM,"taxo")))[c("Kingdom", "Phylum" , "Class"  , "Order"  , "Family"  ,"Genus"  , "Species")])
  }
}
# nrow(Eval_met_mat_all);nrow(Taxoall)
##########################
ggplotdata<-cbind(Eval_met_mat_all,Taxoall)

#number of ASV per group
labs <- c(1:length(na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Archaea",]$Phylum))))#i=1
for (i in 1:length(labs)){#i=1
  test<-ggplotdata[ggplotdata$algo=="GLM" & ggplotdata$Kingdom=="Archaea" & ggplotdata$Phylum==na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Archaea",]$Phylum))[i],]
  labs[i]<-nrow(test[!(is.na(test$R2)),])
}
labs_m<-data.frame(Phylum=na.omit(unique(ggplotdata[ggplotdata$Kingdom=="Archaea",]$Phylum)),labs=labs)

# for (algo in c("GLM","GAM","RF")){
FitEvalP_Am<-PrepareToPlot(datatab=ggplotdata,taxa_level="Phylum",PAAB="AB")

p1 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("R2") + scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), plot.margin = unit(c(0,0,0,0.2),"cm")) + 
  ylim(0,1)
p2 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="GLM" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("Spearman correlation") + scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(margin=margin(r=10)), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0,0,0.2),"cm")) + 
  ylim(0,1)

p3 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)

p4 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="GAM" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+  scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)

p5 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("R2"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=rep("",7)) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)

p6 <- ggplot(FitEvalP_Am[!(is.na(FitEvalP_Am$Phylum)) & FitEvalP_Am$Kingdom=="Archaea" &FitEvalP_Am$algo=="RF" &FitEvalP_Am$variable%in%c("Dpear"),], aes(x=Phylum, y=value))+ geom_violin(scale="width") +
  ylab("")+ scale_x_discrete(labels=c("uncultured"="unassigned")) +
  geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, vjust=1, hjust=1), plot.margin = unit(c(0,0.1,0,0),"cm")) + ylim(0,1)



# par(mfrow=c(4,3))
pdf(paste0("figures/poster/KappaTSS_",algo,".pdf"))
library(grid)
tGLM <- textGrob("GLM",hjust=-.2)
tGAM <- textGrob("GAM",hjust=-.1)
tRF <- textGrob(" RF",hjust=-.2)
pdf(paste0("figures/poster/figure1.pdf"))
grid.arrange(grobs=list(p1,p2,p3,p4,p5,p6,tGLM,tGAM,tRF), 
             layout_matrix = rbind(c(7,8,9),
                                   c(1,3,5),
                                   c(2,4,6)),
             heights=c(0.1,0.5,0.5))
dev.off()
plot(p3)








############################## WORK in prog ##################################

##############################$
#Test Carte

#Who with only spatial variable
load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
load("../../spatial_data/ENVstack.Rda")
ENVstack_df <- as.data.frame(ENVstack)
names(ENVstack)
names(ENVdata)
plot(ENVstack$pH)
  




#checking Eval metrics
# files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <- gsub("VarSel","",files)
VarSel_list<-c()
# PAAB = "PA"
# GtoM="BA"
# algo="GLM"
for(file in files){#file=files[8] #OTU23200
  
  load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/VarSel/VarSel",file)) #array 1
  temporary<-lapply(VarSel$ranking,function(X){
    if(length(X)==2){
      return(c("var1"=NA,"var2"=NA,"var3"=NA,"var4"=NA,"var5"=NA,"var6"=NA,"var7"=NA,"var8"=NA,"var9"=NA,"var10"=NA,"var11"=NA,"var12"=NA,"var13"=NA,"var14"=NA,"var15"=NA))
    }else{
      Xbody <- as.vector(X["var"])
      Xbody <- unlist(Xbody)
      length(Xbody) <- 15
      names(Xbody) <- c("var1","var2","var3","var4","var5","var6","var7","var8","var9","var10","var11","var12","var13","var14","var15")
      return(Xbody)
    }
  })
  varSel_tab <- do.call("rbind",temporary)
  VarSel_list <-rbind(VarSel_list,varSel_tab)

}

SpatOTU <- apply(VarSel_list,1,function(X){sum(X%in%names(ENVstack),na.rm=TRUE)/sum(!(is.na(X)))})
Spat_OTU <- VarSel_list[names(SpatOTU),][which(SpatOTU==1),]

load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo.Rda"))
BAtaxo ==length(SpatOTU)
BAtaxo$OTU<-1:nrow(BAtaxo)
BAtaxoSpat<-BAtaxo[BAtaxo$OTU%in%as.numeric(gsub("OTU","",rownames(Spat_OTU))),]
BAtaxoSpat[BAtaxoSpat$Kingdom=="Archaea",]$Phylum
BAtaxoSpat[BAtaxoSpat$Kingdom=="Archaea" & BAtaxoSpat$Phylum=="Asgardaeota",]
#9480 42656
VarSel_list[9480,]
VarSel_list[42656,]

library(h2o)
h2o.init()
OTUtoRun<- 42656
load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))

ENVdataSpat <- ENVdata[,colnames(ENVdata)%in%names(ENVstack)]

training_frame1 <- as.h2o( cbind(ENVdata,ASVdata=OTUdata[,OTUtoRun]) )
training_frame2 <- as.h2o( cbind(ENVdataSpat,ASVdata=OTUdata[,OTUtoRun]) )

Mglm1 <- h2o.glm(y = "ASVdata", training_frame = training_frame1,nfolds = 10)
glmPerf1 <- h2o.performance(Mglm1)
glmVarImp1 <- h2o.varimp(Mglm1)

metricsglm1 <-unlist(glmPerf1@metrics[c("RMSE","nobs", "r2","logloss", "AUC", "Gini", "null_deviance", "residual_deviance", "AIC")])

devs <- c(metricsglm1["null_deviance"], metricsglm1["residual_deviance"])
Eval_glm1 <- (devs[1]-devs[2])/devs[1]*100 #Deviance explained percentage (close to R2)
names(Eval_glm1) <- "expl.dev"



#table with data and prediction of the model
predF1_1<- as.data.frame(h2o.predict(Mglm1,as.h2o(ENVstack_df))[c(1,3)])

names(predF1_1) <- c("pred","prob")

confmatglm1 <- glmPerf1@metrics$cm$table
#maxTSS
maxTSSthresh1 <- glmPerf1@metrics$thresholds_and_metric_scores$threshold[which(glmPerf1@metrics$thresholds_and_metric_scores$recall+glmPerf1@metrics$thresholds_and_metric_scores$specificity==max(glmPerf1@metrics$thresholds_and_metric_scores$recall+glmPerf1@metrics$thresholds_and_metric_scores$specificity))]
predF1_1$predTSS<-predF1_1$prob>=maxTSSthresh1
coordinates(predF1_1)<-coordinates(ENVstack)
gridded(predF1_1)<-TRUE
predF1_1_rast<-stack(predF1_1)
values(predF1_1_rast$pred)<-values(predF1_1_rast$pred)-1


library(terra)
library(raster)
library(rgdal)
lakes<-readOGR("../../../../../common/50_data/GeoData_RECHALP_GRUYERE/hydro/lake1903+.shp")
leman <-readOGR("../../../../../common/50_data/GeoData_RECHALP_GRUYERE/LemanPart.shp")

rivers <- readOGR("../../../../../common/50_data/GeoData_RECHALP_GRUYERE/hydro/river1903+.shp")
plot(leman)


# load("../../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")
# OTUdataSpat_Pres<-as.data.frame(OTUdata[OTUdata[,OTUtoRun]==1,OTUtoRun])
# coordinates(OTUdataSpat_Pres) <- dataSoil[dataSoil$sampleNameBA %in% rownames(OTUdataSpat_Pres), ][c("x","y")]
# OTUdataSpat_Abs<-as.data.frame(OTUdata[OTUdata[,OTUtoRun]==0,OTUtoRun])
# coordinates(OTUdataSpat_Abs) <- dataSoil[dataSoil$sampleNameBA %in% rownames(OTUdataSpat_Abs), ][c("x","y")]

library(scales)
plot(mask(ENVstack$pH,leman,inverse=TRUE))
plot(ENVstack$pH)

plot(mask(mask(mask(predF1_1_rast$predTSS,lakes,inverse=TRUE),ENVstack$pH),leman,inverse=TRUE),col=c("white","orangered4"))
legend(x=2546000,y=1120000,legend=c("predicted presence","predicted absence"),fill=c("white","orangered4"),bty="n",)
# plot(OTUdataSpat_Pres,add=TRUE,col="black",pch=20)
# plot(OTUdataSpat_Abs,add=TRUE,col="lightgrey",pch=20,cex=1,alpha=50)


plot(mask(mask(mask(predF1_1_rast$prob,lakes,inverse=TRUE),ENVstack$pH),leman,inverse=TRUE))





Mglm2 <- h2o.glm(y = "ASVdata", training_frame = training_frame2,nfolds = 10)
glmPerf2 <- h2o.performance(Mglm2)
glmVarImp2 <- h2o.varimp(Mglm2)

metricsglm2 <-unlist(glmPerf2@metrics[c("RMSE","nobs", "r2","logloss", "AUC", "Gini", "null_deviance", "residual_deviance", "AIC")])

devs <- c(metricsglm2["null_deviance"], metricsglm2["residual_deviance"])
Eval_glm2 <- (devs[1]-devs[2])/devs[1]*100 #Deviance explained percentage (close to R2)
names(Eval_glm2) <- "expl.dev"

#table with data and prediction of the model
predF1_2<- as.data.frame(h2o.predict(Mglm2,as.h2o(ENVstack_df))[c(1,3)])

names(predF1_2) <- c("pred","prob")

confmatglm2 <- glmPerf2@metrics$cm$table
#maxTSS
maxTSSthresh2 <- glmPerf2@metrics$thresholds_and_metric_scores$threshold[which(glmPerf2@metrics$thresholds_and_metric_scores$recall+glmPerf2@metrics$thresholds_and_metric_scores$specificity==max(glmPerf2@metrics$thresholds_and_metric_scores$recall+glmPerf2@metrics$thresholds_and_metric_scores$specificity))]
predF1_2$predTSS<-predF1_2$prob>=maxTSSthresh2
coordinates(predF1_2)<-coordinates(ENVstack)
gridded(predF1_2)<-TRUE
predF1_2_rast<-stack(predF1_2)
values(predF1_2_rast$pred)<-values(predF1_2_rast$pred)-1

plot(mask(mask(mask(predF1_2_rast$predTSS,lakes,inverse=TRUE),ENVstack$pH),leman,inverse=TRUE),col=c("white","orangered4"))
legend(x=2546000,y=1120000,legend=c("predicted presence","predicted absence"),fill=c("white","orangered4"),bty="n",)
plot(mask(mask(mask(predF1_1_rast$prob,lakes,inverse=TRUE),ENVstack$pH),leman,inverse=TRUE))





VarSel_list[42656,]
#not god enough for abundance model :(


#Methods
# load("../../spatial_data/ENVstack.Rda")

plot(mask(mask(ENVstack$ndvi,lakes,inverse=TRUE),leman,inverse=TRUE))
plot(ENVstack$aspect,col=colorRampPalette(c("white","black"))(100))
hill <- raster(paste0("Z:\\common\\50_data\\GeoData_RECHALP_GRUYERE\\topo\\hillshade_mean2m1903+.tif"))%>%
  crop(.,extent(c(xmin=2552000,xmax=2587000,ymin=1114000,ymax=1157000)))%>%

  mask(.,biovars[[1]])
hill <- projectRaster(hill,ENVstack)
plot(mask(mask(mask(hill,ENVstack$bio1_t),lakes,inverse=TRUE),leman,inverse=TRUE),col=colorRampPalette(c("grey77","grey50"))(100))
load("../../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")

coordinates(dataSoil)<-dataSoil[c("x","y")]
crs(dataSoil)<-crs(ENVstack)
plot(dataSoil,add=TRUE,pch=20, col="blue")


hillCH<-readRDS("../../spatial_data/Valpar/ch_topo_alti3d2016_pixel_dem_mean2m_25.rds")
plot(hillCH,col=colorRampPalette(c("grey77","grey50"))(100))
