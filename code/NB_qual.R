PAAB="PA"
group="BA"
Mod="GLM"
load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval.Rda"))
load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
if(nrow(BAtaxo_grouped_phylum)!=nrow(Eval)){#securité
  message("problem nrow")
}
load(paste0("PA/",group,"/data/ENVdata.Rda"))
load(paste0("PA/",group,"/data/OTUdata.Rda"))
all(rownames(OTUdata)==rownames(ENVdata))
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
ENVdatacont<-apply(ENVdata[contvar],2,function(X){X/(max(X))})
ENVOTUdata<-cbind(ENVdatacont,OTUdata)
evalmet_BA<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",],Eval[BAtaxo_grouped_phylum$Kingdom=="Bacteria",])
evalmet_BA$Phylum<-factor(evalmet_BA$Phylum,levels=c(names(table(evalmet_BA$Phylum[!(evalmet_BA$Phylum%in%c("Others","unclassified_Bacteria"))])),"Others","unclassified_Bacteria"))
evalmet_AR<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",],Eval[BAtaxo_grouped_phylum$Kingdom=="Archaea",])
evalmet_AR$Phylum<-factor(evalmet_AR$Phylum,levels=c(names(table(evalmet_AR$Phylum[!(evalmet_AR$Phylum%in%c("Others","unclassified_Archaea"))])),"Others","unclassified_Archaea"))


library(ecospat)
# library(ade4)
# naxis<-9
# pca.env<-dudi.pca(df = apply(ENVdata, 2, as.numeric), scannf = FALSE, nf = naxis)
# pca.env$li
# ENVOTUdata<-cbind(pca.env$li,OTUdata)
# niche<-ecospat.nichePOSNB(df=ENVOTUdata,colvar=c(1:naxis),colfreq=c((naxis+1):ncol(ENVOTUdata)))
niche<-ecospat.nichePOSNB(df=ENVOTUdata,colvar=c(1:ncol(ENVdatacont)),colfreq=c((ncol(ENVdatacont)+1):ncol(ENVOTUdata)))
# ecospat.nichePOSNB
ENVOTUdata[1:10,1:ncol(ENVdata)]
load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
NBvsPQ2<-data.frame(niche[match(evalmet_BA$Seq,rownames(niche)),],TSSadj=evalmet_BA$TSS_adj,BAtaxo_grouped_phylum[match(evalmet_BA$Seq,BAtaxo_grouped_phylum$Seq),"Phylum"])
# NBvsPQ2[1:10,]

par(mfrow=c(4,4))
for (cov in contvar){#cov="ThinSand"
  NBvsPQ<-data.frame(NicheBreadth=NBvsPQ2[,paste0(cov,"_nb")],TSSadj=NBvsPQ2[,"TSSadj"])
  plot(NBvsPQ$NicheBreadth,NBvsPQ$TSSadj,main=cov, xlab="",ylab="")
  abline(glm(NBvsPQ$TSSadj~NBvsPQ$NicheBreadth))
  # NBPQ_loess<-loess(NBvsPQ$TSSadj~NBvsPQ$NicheBreadth,span=1.5)
  # xfit=seq(from=min(NBvsPQ$NicheBreadth),to=max(NBvsPQ$NicheBreadth),length.out=100)
  # points(xfit,yfit,type="l",lwd=2,col="red")
}
par(mfrow=c(1,1))
par(mfrow=c(3,3))
for (axis in names(pca.env$li)){

  NBvsPQ<-data.frame(NicheBreadth=NBvsPQ2[,paste0(axis,"_nb")],TSSadj=NBvsPQ2[,"TSSadj"])
  plot(NBvsPQ$NicheBreadth,NBvsPQ$TSSadj,main=axis)
  abline(glm(NBvsPQ$TSSadj~NBvsPQ$NicheBreadth))
  NBPQ_loess<-loess(NBvsPQ$TSSadj~NBvsPQ$NicheBreadth,span=1.5)
  xfit=seq(from=min(NBvsPQ$NicheBreadth),to=max(NBvsPQ$NicheBreadth),length.out=100)
  yfit=predict(NBPQ_loess,newdata=xfit)
  points(xfit,yfit,type="l",lwd=2,col="red")
}
par(mfrow=c(1,1))


library(hypervolume)
?hypervolume()



testGLM<-glm(NBvsPQ$TSSadj~NBvsPQ$NicheBreadth)
testGLM$deviance
testGLM$null.deviance
summary(testGLM)


nrow(niche[match(evalmet_BA$Seq,rownames(niche)),"pH_nb"])
BAtaxo_grouped_phylum$Seq
all(evalmet_BA$Seq==rownames(niche))
length(evalmet_BA$Seq);length(rownames(niche))
# PAAB="PA"
for(Mod in c("GBM","GAM","RF","GLM")){#Mod="GLM"
  # group="PR"
  for (group in c("PR","BA","FU")){#group="BA"
    load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval.Rda"))
    
    if(group=="BA"){
      # load(paste0("../../ASV_data/ASV_taxo/",substr(group,1,2),"taxo.Rda"))
      # Eval<- cbind(BAtaxo,Eval)
      load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
      if(nrow(BAtaxo_grouped_phylum)!=nrow(Eval)){#securité
        message("problem nrow")
      }
      evalmet_BA<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",],Eval[BAtaxo_grouped_phylum$Kingdom=="Bacteria",])
      evalmet_BA$Phylum<-factor(evalmet_BA$Phylum,levels=c(names(table(evalmet_BA$Phylum[!(evalmet_BA$Phylum%in%c("Others","unclassified_Bacteria"))])),"Others","unclassified_Bacteria"))
      evalmet_AR<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",],Eval[BAtaxo_grouped_phylum$Kingdom=="Archaea",])
      evalmet_AR$Phylum<-factor(evalmet_AR$Phylum,levels=c(names(table(evalmet_AR$Phylum[!(evalmet_AR$Phylum%in%c("Others","unclassified_Archaea"))])),"Others","unclassified_Archaea"))
      
      
      
      
      
      
      goodASV_phylum_OK_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_02_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.2&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.2&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_04_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.4&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.4&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_06_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.6&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.6&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_08_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.8&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.8&!(is.na(evalmet_BA$TSS))]))
      # number instead of proportions
      # summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])
      # Acidobacteriota      Actinobacteriota          Bacteroidota      Bdellovibrionota           Chloroflexi            Firmicutes 
      # 1171                   556                   368                    72                   528                    77 
      # Myxococcota       Patescibacteria       Planctomycetota        Proteobacteria     Verrucomicrobiota                Others 
      # 194                   189                   576                  1404                   285                   221 
      # unclassified_Bacteria 
      # 1334 
      # summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))]) #
      # Crenarchaeota        Euryarchaeota        Halobacterota        Iainarchaeota               Others unclassified_Archaea 
      # 20                    0                    0                    0                    2                    2 
      
      goodASV_phylum_OK_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_02_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.2&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.2&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_04_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.4&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.4&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_06_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.6&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.6&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_08_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.8&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.8&!(is.na(evalmet_AR$TSS))]))
      
      #proportion of each phylum in all modeled ASV pool
      AllASV_phylum_BA<-summary(factor(evalmet_BA$Phylum))/sum(summary(factor(evalmet_BA$Phylum)))
      AllASV_phylum_AR<-summary(factor(evalmet_AR$Phylum))/sum(summary(factor(evalmet_AR$Phylum)))
      
      comparison_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_OK_BA,goodASV_phylum_02_BA,goodASV_phylum_04_BA,goodASV_phylum_06_BA,goodASV_phylum_08_BA))
      comparison_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_OK_AR,goodASV_phylum_02_AR,goodASV_phylum_04_AR,goodASV_phylum_06_AR,goodASV_phylum_08_AR))
      
      #proportion of good modeled ASV in each phylum
      goodASV_OK_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/summary(factor(evalmet_BA$Phylum))
      goodASV_OK_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/summary(factor(evalmet_AR$Phylum))
      
      #proportion of verygood modeled ASV in each phylum
      goodASV_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/summary(factor(evalmet_BA$Phylum))
      if(length(goodASV_05_BA)==0){
        goodASV_05_BA<-goodASV_OK_BA*0
      }
      goodASV_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/summary(factor(evalmet_AR$Phylum))
      if(length(goodASV_05_AR)==0){
        goodASV_05_AR<-goodASV_OK_AR*0
      }
      
      molten_good_BA <- melt(data.frame(goodASV_OK=goodASV_OK_BA,phylum=names(goodASV_OK_BA),goodASV_05=goodASV_05_BA))
      molten_good_BA$phylum<-factor(molten_good_BA$phylum,levels=levels(evalmet_BA$Phylum))
      labels_BA <- paste(unique(molten_good_BA$phylum), " ( n =",c(table(evalmet_BA$Phylum)[-which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))],table(evalmet_BA$Phylum)[which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))]), ")")
      molten_good_AR <- melt(data.frame(goodASV_OK=goodASV_OK_AR,phylum=names(goodASV_OK_AR),goodASV_05=goodASV_05_AR))
      molten_good_AR$phylum<-factor(molten_good_AR$phylum,levels=levels(evalmet_AR$Phylum))
      labels_AR <- paste(unique(molten_good_AR$phylum), " ( n =",c(table(evalmet_AR$Phylum)[-which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))],table(evalmet_AR$Phylum)[which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))]), ")")
      
      p1_BA<-ggplot(comparison_BA, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_BA" = "Whole dataset", "goodASV_phylum_OK_BA" = "TSS >TSSnull", "goodASV_phylum_02_BA" = "TSSadj>0.2", "goodASV_phylum_04_BA" = "TSSadj>0.4", "goodASV_phylum_06_BA" = "TSSadj>0.6", "goodASV_phylum_08_BA" = "TSSadj>0.8"))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      p2_BA<-ggplot(molten_good_BA,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("Proportion of models")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA,
                          label=labels_BA)+ theme_bw() +
        theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS >TSSnull", "goodASV_05" = "TSSadj>0.5"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset[group]),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"),
              legend.background = element_blank())
      # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
      # p3_BA<-ggplot(comparison_05_BA, aes(x=Var1, y=value, fill=Var2)) +
      #   geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$BA)
      
      # p4_BA<-ggplot(data.frame(value=goodASV_05_BA,phylum=names(goodASV_05_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
      #   ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$BA)
      # grid.arrange(p3_BA,p4_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_05"))
      p1_AR<-ggplot(comparison_AR, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_AR" = "Whole dataset", "goodASV_phylum_OK_AR" = "TSS >TSSnull", "goodASV_phylum_02_AR" = "TSSadj>0.2", "goodASV_phylum_04_AR" = "TSSadj>0.4", "goodASV_phylum_06_AR" = "TSSadj>0.6", "goodASV_phylum_08_AR" = "TSSadj>0.8"))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      
      p2_AR<-ggplot(molten_good_AR,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR,
                          label=labels_AR)+ 
        theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS >TSSnull", "goodASV_05" = "TSSadj>0.5"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset["AR"]),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"))
      
      # nOTU_phylum_BA<- c(table(evalmet_BA$Phylum)[-which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))],table(evalmet_BA$Phylum)[which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))])
      # nOTU_phylum_AR<- c(table(evalmet_AR$Phylum)[-which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))],table(evalmet_AR$Phylum)[which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))])
      # ngoodOTU_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])
      # ngoodOTU_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])
      # table_phylum <-data.frame(BA_nOTU=nOTU_phylum_BA,AR_nOTU=nOTU_phylum_AR,BA_ngood=ngoodOTU_BA,AR_ngood=ngoodOTU_AR)
      # plot(nOTU_phylum_BA,goodASV_05_BA)
      # plot(nOTU_phylum_AR,goodASV_05_AR)
      
    }else{
      
      load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
      taxo_grouped_phylum_temp<-eval(parse(text=paste0(group,"taxo_grouped_phylum")))
      #remove not fungi
      table(taxo_grouped_phylum_temp$Phylum)
      
      if(nrow(taxo_grouped_phylum_temp)!=nrow(Eval)){#securité
        message("problem")
      }
      
      
      if(group=="PR"){
        evalmet<-cbind(Eval[taxo_grouped_phylum_temp$Phylum!="Not_Protist",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Protist",])
        evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others")))],"Others"))
      }
      if(group=="FU"){
        evalmet<-cbind(Eval[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",])
        evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others","Unidentified_Fungi")))],"Others","Unidentified_Fungi"))
      }
      table(evalmet$Phylum)
      evalmet$Phylum<-droplevels(evalmet$Phylum)
      #mean
      # mean(evalmet$TSS_adj,na.rm=TRUE) #FU 0.19  #PR  0.04
      # sd(evalmet$TSS_adj,na.rm=TRUE) #FU 0.20  #PR   0.14
      # levels(evalmet$Phylum)
      # evalmetGLM<-evalmet
      # fitmetGLM<-fitmet
      goodASV_phylum_OK<-summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))]))
      goodASV_phylum_02<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.2&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.2&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_04<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.4&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.4&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_05<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_06<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.6&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.6&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_08<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.8&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.8&!(is.na(evalmet$TSS_adj))]))
      # number instead of proportions
      # summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS))]) 
      #Fungi GLM
      # Ascomycota     Basidiomycota     Glomeromycota Mortierellomycota            Others      Unidentified 
      # 589                 43                 37                 89                  3                203 
      #Protists GLM
      # SAR      Amoebozoa Archaeplastida         Others           NA's 
      #        1              0              1              0             25 
      
      
      AllASV_phylum<-summary(factor(evalmet$Phylum))/sum(summary(factor(evalmet$Phylum)))
      
      comparison<-melt(rbind(AllASV_phylum,goodASV_phylum_OK,goodASV_phylum_02,goodASV_phylum_04,goodASV_phylum_06,goodASV_phylum_08))
      comparison$value[is.na(comparison$value)]<-0
      
      #proportion of good modeled ASV in each phylum
      goodASV_OK<-summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))])/summary(factor(evalmet$Phylum))
      
      #proportion of verygood modeled ASV in each phylum
      goodASV_05<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))])/summary(factor(evalmet$Phylum))
      if(length(goodASV_05)==0){
        goodASV_05<-goodASV_OK*0
      }
      
      
      molten_good <- melt(data.frame(goodASV_OK=goodASV_OK,phylum=names(goodASV_OK),goodASV_05=goodASV_05))
      molten_good$phylum<-factor(molten_good$phylum,levels=levels(evalmet$Phylum))
      # labels <- paste(unique(molten_good$phylum), " (n =",c(table(evalmet$Phylum)[-which(names(table(evalmet$Phylum))%in%c("Others","Fungi_unidentified"))],table(evalmet$Phylum)[which(names(table(evalmet$Phylum))%in%c("Others","Fungi_unidentified"))]), ")",sep="")
      labels <- paste(unique(molten_good$phylum), " (n =",table(evalmet$Phylum), ")",sep="")
      # sum(!(is.na(evalmet$TSS)))
      assign(paste0("p1_", group), ggplot(comparison, aes(x=Var1, y=value, fill=Var2)) +
               geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]])+ theme_bw() +
               scale_x_discrete(labels=c("AllASV_phylum" = "Whole dataset", "goodASV_phylum_OK" = "TSS > TSSnull", "goodASV_phylum_02" = "TSSadj>0.2", "goodASV_phylum_04" = "TSSadj>0.4", "goodASV_phylum_06" = "TSSadj>0.6", "goodASV_phylum_08" = "TSSadj>0.8"))+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     legend.title = element_blank(),
                     axis.title.x = element_blank())
      )#p1_FU
      if(group=="PR"){plotmarg<-unit(c(0,0.5,0,0),"cm")}else{plotmarg<-unit(c(0,0,0,0.5),"cm")}
      assign(paste0("p2_", group), ggplot(molten_good,aes(x=variable,y=value,fill=phylum,group=phylum))+
               geom_bar(position="dodge",stat="identity")+
               ylim(0,1) + xlab("") + ylab(ifelse(group=="PR","","Proportion of models"))+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]],
                                 label=labels)+ 
               theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS > TSSnull", "goodASV_05" = "TSSadj > 0.5"))+ 
               theme_classic()+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=10),
                     axis.text.y = element_text(size=6),
                     axis.text.x = element_text(size=10),
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     legend.title = element_text(size=7,colour = setcolors$dataset[group]),
                     legend.text = element_text(size=6),
                     legend.key.size = unit(0.3,"line"),
                     legend.position = c(.95,1),
                     legend.justification = c("right", "top"))
      )
    }#p2_FU
    # nOTU_phylum<- table(evalmet$Phylum)
    # ngoodOTU<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS))])
    # plot(as.numeric(nOTU_phylum[-length(nOTU_phylum)]),as.numeric(goodASV_05[-length(nOTU_phylum)]))
    
  }
  
  pdf(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Pylum_proportion",Mod,"_",PAAB,".pdf"),width=1961,height=1500)
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_proportion",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_relative_richness",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
}    
