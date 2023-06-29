#For presence- Absence
#Creates figures of proportion of each phylum in overall data and in models better than  null
#Creates figures of proportion of each phylum ASV that have better models than null
#Do the same for TSS > 0.5
#For abundances
#Creates figures of proportion of each phylum in overall data and in models with cor > 0.2
#Creates figures of proportion of each phylum ASV that have model with cor >0.2
#Do the same for cor > 0.5
#whogoodmodels
library(tidyverse)
library(ggplot2)
library(reshape2) 
library(grid)
library(gridExtra)


#For presence- Absence
PAAB="PA"



for(Mod in c("GLM","GBM","GAM","RF")){#Mod="GLM"
  # group="PR"
for (group in c("PR","FU","BA")){#group="FU"
# for (PAAB in c("PA","AB")){


# Mod="GLM"

  print(paste0(group,PAAB,Mod))
load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval.Rda"))

setcolors<-list(BA=c("#a6cee3", "#1f78b4" ,"#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99"),
                AR=c("#4053d3", "#ddb310", "#b51d14", "#00beff", "#fb49b0", "#00b25d", "#cacaca"),
                FU=c("#4053d3", "#ddb310", "#b51d14", "#00beff", "#fb49b0", "#00b25d"),
                PR=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"))
if(group=="BA"){
  load(paste0("../../ASV_data/ASV_taxo/",substr(group,1,2),"taxo.Rda"))
  Eval<- cbind(Eval,BAtaxo)
  evalmet_BA<-Eval[which(Eval$Kingdom=="Bacteria"),]
  evalmet_AR<-Eval[which(Eval$Kingdom=="Archaea"),]
  load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
  load(paste0("../../ASV_data/ASV_taxo/ARtaxo_grouped_phylum.Rda"))
# nrow(evalmet_BA)
# nrow(taxo_grouped_phylum_BA)
  taxo_grouped_phylum_BA<-BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",]
  taxo_grouped_phylum_AR<-ARtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",]
  
  # unique(taxo_grouped_phylum_BA$Species)
  #proportion of each phylum in good modeled ASV pool
  goodASV_phylum_02_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))]))
  goodASV_phylum_02_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))]))
  
  #proportion of each phylum in very good modeled ASV pool
  goodASV_phylum_05_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))]))
  # listofgoodsp<-unique(factor(taxo_grouped_phylum_BA$Species)[evalmet_BA$TSS>0.5&!(is.na(evalmet_BA$TSS))])
  # listofgoodsp<-listofgoodsp[substr(listofgoodsp,1,10)!="uncultured"]
  # listofgoodsp[grep("Pseudomonas",listofgoodsp)]
  # Pseudomonas<-c(taxo_grouped_phylum_BA[taxo_grouped_phylum_BA$Species=="Pseudomonas_syringae_pv._broussonetiae"&!(is.na(taxo_grouped_phylum_BA$Species)),],evalmet_BA)
  # 
  
  goodASV_phylum_05_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))]))

  #proportion of each phylum in all modeled ASV pool
  AllASV_phylum_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum))/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)))
  AllASV_phylum_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum))/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)))
  
  comparison_02_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_02_BA))
  comparison_05_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_05_BA))
  comparison_02_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_02_AR))
  comparison_05_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_05_AR))
  
  #proportion of good modeled ASV in each phylum
  goodASV_02_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/summary(factor(taxo_grouped_phylum_BA$Phylum))
  goodASV_02_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/summary(factor(taxo_grouped_phylum_AR$Phylum))
  
  #proportion of verygood modeled ASV in each phylum
  goodASV_05_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/summary(factor(taxo_grouped_phylum_BA$Phylum))
  if(length(goodASV_05_BA)==0){
    goodASV_05_BA<-goodASV_02_BA*0
  }
  goodASV_05_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/summary(factor(taxo_grouped_phylum_AR$Phylum))
  if(length(goodASV_05_AR)==0){
    goodASV_05_AR<-goodASV_02_AR*0
  }

  p1_BA<-ggplot(comparison_02_BA, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
    scale_fill_manual("Bacteria phylum",values = setcolors$BA)+ theme_bw() +
    scale_x_discrete(labels=c("AllASV_phylum_BA" = "Phyla relative \n gamma diversity", "goodASV_phylum_02_BA" = "Phyla relative diversity \n across good models"))+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=15),
          legend.title = element_blank())
  p2_BA<-ggplot(data.frame(value=goodASV_02_BA,phylum=names(goodASV_02_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
    ylim(0,1) + xlab("") + ylab("Proportion of all ASV")+
    scale_fill_manual("Bacteria phylum",values = setcolors$BA)+ theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.title.y = element_text(size=25),
          axis.text.y = element_text(size=15),
          axis.text.x = element_blank(),
          legend.title = element_text(size=25),
          legend.text = element_text(size=20),
          legend.key.size = unit(2,"line"))
  # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
  p3_BA<-ggplot(comparison_05_BA, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$BA)
  
  p4_BA<-ggplot(data.frame(value=goodASV_05_BA,phylum=names(goodASV_05_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
    ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$BA)
  # grid.arrange(p3_BA,p4_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_05"))
  p1_AR<-ggplot(comparison_02_AR, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$AR)+ scale_x_discrete(labels=c("AllASV_phylum_AR" = "Phyla relative \n gamma diversity", "goodASV_phylum" = "Phyla relative diversity \n across good models"))
  
  p2_AR<-ggplot(data.frame(value=goodASV_02_AR,phylum=names(goodASV_02_AR)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
    ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$AR)
  # grid.arrange(p1_AR,p2_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_02"))
  p3_AR<-ggplot(comparison_05_AR, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$AR)
  
  p4_AR<-ggplot(data.frame(value=goodASV_05_AR,phylum=names(goodASV_05_AR)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
    ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors$AR)
  # grid.arrange(p3_AR,p4_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_05"))

  
}else{
  
  load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
  taxo_grouped_phylum<-eval(parse(text=paste0(group,"taxo_grouped_phylum")))
  
  # evalmetGLM<-evalmet
  # fitmetGLM<-fitmet
  
  #how many good ASV per phylum
  # summary(factor(PRtaxo_grouped_phylum[evalmetGLM$TSS>0.2&!(is.na(evalmetGLM$TSS)),]$Phylum))
  
  #how many very good ASV per phylum
  # summary(factor(PRtaxo_grouped_phylum[evalmetGLM$TSS>0.5&!(is.na(evalmetGLM$TSS)),]$Phylum))
  
  #how many modeled ASV per phylum
  # summary(factor(PRtaxo_grouped_phylum$Phylum))
  
  #proportion of each phylum in good modeled ASV pool
  #goodmodeledASVpooltaxo_grouped_phylum[Eval$TSS_sign&!(is.na(Eval$TSS_sign)),]
  goodmodeledASVpool<-taxo_grouped_phylum[Eval$TSS_sign&!(is.na(Eval$TSS_sign)),]
  goodASV_phylum<-summary(factor(goodmodeledASVpool$Phylum))/sum(summary(factor(goodmodeledASVpool$Phylum)))
  #proportion of each phylum in very good modeled ASV pool
  goodASV_phylum_05<-summary(factor(taxo_grouped_phylum[Eval$TSS_adj>0.5&!(is.na(Eval$TSS)),]$Phylum))/sum(summary(factor(taxo_grouped_phylum[Eval$TSS_adj>0.5&!(is.na(Eval$TSS)),]$Phylum)))
  
  #proportion of each phylum in all modeled ASV pool
  AllASV_phylum<-summary(factor(taxo_grouped_phylum$Phylum))/sum(summary(factor(taxo_grouped_phylum$Phylum)))
  
  comparison<-melt(rbind(AllASV_phylum,goodASV_phylum))
  comparison_05<-melt(rbind(AllASV_phylum,goodASV_phylum_05))
  
  #proportion of good modeled ASV in each phylum
  goodASV<-summary(factor(goodmodeledASVpool$Phylum))/summary(factor(taxo_grouped_phylum$Phylum))
  
  #proportion of verygood modeled ASV in each phylum
  goodASV_05<-summary(factor(taxo_grouped_phylum[Eval$TSS_adj>0.5&!(is.na(Eval$TSS)),]$Phylum))/summary(factor(taxo_grouped_phylum$Phylum))
  if(length(goodASV_05)==0){
    goodASV_05<-goodASV*0
  }

  assign(paste0("p1_", group), ggplot(comparison, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("") + scale_x_discrete(labels=c("AllASV_phylum" = "Phyla relative \n gamma diversity", "goodASV_phylum" = "Phyla relative diversity \n across good models"))+
    scale_fill_manual(values = setcolors[[group]]))
  assign(paste0("p2_", group), ggplot(data.frame(value=goodASV,phylum=names(goodASV)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
           ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors[[group]]))
  # grid.arrange(p1,p2,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
  assign(paste0("p3_", group), ggplot(comparison_05, aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors[[group]]))
  assign(paste0("p4_", group), ggplot(data.frame(value=goodASV_05,phylum=names(goodASV_05)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
           ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
    scale_fill_manual(values = setcolors[[group]]))
  # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
  
}
}
  
  pdf(file=paste0("figures/test_phylum/goodmodels/",Mod,"_",PAAB,".pdf"))
  grid.arrange(p1_BA,p2_BA,p1_AR,p2_AR,p1_FU,p2_FU,p1_PR,p2_PR,nrow=4,top=paste0(Mod,"_",PAAB))
  grid.arrange(p3_BA,p4_BA,p3_AR,p4_AR,p3_FU,p4_FU,p3_PR,p4_PR,nrow=4,top=paste0(Mod,"_",PAAB,"_05"))
  dev.off()
}


PAAB="AB"
#For Abundances



for (group in c("PR","FU","BA")){#group="PR"
  # for (PAAB in c("PA","AB")){
  for(Mod in c("GLM","GBM","GAM","RF")){#Mod="GLM"
    # group="PR"
    
    # Mod="GLM"
    
    print(paste0(group,PAAB,Mod))
    load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval_Met.Rda"))
    
    
    if(group=="BA"){
      load("../../ASV_data/ASV_taxo/BAtaxo.Rda")
      evalmet_BA<-evalmet[BAtaxo$Kingdom=="Bacteria" & !(is.na(BAtaxo$Kingdom=="Bacteria")),]
      evalmet_AR<-evalmet[BAtaxo$Kingdom=="Archaea" & !(is.na(BAtaxo$Kingdom=="Bacteria")),]
      load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
      load(paste0("../../ASV_data/ASV_taxo/ARtaxo_grouped_phylum.Rda"))
      
      taxo_grouped_phylum_BA<-BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",]
      taxo_grouped_phylum_AR<-ARtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",]
      
      
      #proportion of each phylum in good modeled ASV pool
      goodASV_phylum_02_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.2&!(is.na(evalmet_BA$Dpear))])/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.2&!(is.na(evalmet_BA$Dpear))]))
      goodASV_phylum_02_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.2&!(is.na(evalmet_AR$Dpear))])/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.2&!(is.na(evalmet_AR$Dpear))]))
      
      #proportion of each phylum in very good modeled ASV pool
      goodASV_phylum_05_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.5&!(is.na(evalmet_BA$Dpear))])/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.5&!(is.na(evalmet_BA$Dpear))]))
      goodASV_phylum_05_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.5&!(is.na(evalmet_AR$Dpear))])/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.5&!(is.na(evalmet_AR$Dpear))]))
      
      #proportion of each phylum in all modeled ASV pool
      AllASV_phylum_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum))/sum(summary(factor(taxo_grouped_phylum_BA$Phylum)))
      AllASV_phylum_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum))/sum(summary(factor(taxo_grouped_phylum_AR$Phylum)))
      
      comparison_02_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_02_BA))
      comparison_05_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_05_BA))
      comparison_02_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_02_AR))
      comparison_05_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_05_AR))
      
      #proportion of good modeled ASV in each phylum
      goodASV_02_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.2&!(is.na(evalmet_BA$Dpear))])/summary(factor(taxo_grouped_phylum_BA$Phylum))
      goodASV_02_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.2&!(is.na(evalmet_AR$Dpear))])/summary(factor(taxo_grouped_phylum_AR$Phylum))
      
      #proportion of verygood modeled ASV in each phylum
      goodASV_05_BA<-summary(factor(taxo_grouped_phylum_BA$Phylum)[evalmet_BA$Dpear>0.5&!(is.na(evalmet_BA$Dpear))])/summary(factor(taxo_grouped_phylum_BA$Phylum))
      if(length(goodASV_05_BA)==0){
        goodASV_05_BA<-goodASV_02_BA*0
      }
      goodASV_05_AR<-summary(factor(taxo_grouped_phylum_AR$Phylum)[evalmet_AR$Dpear>0.5&!(is.na(evalmet_AR$Dpear))])/summary(factor(taxo_grouped_phylum_AR$Phylum))
      if(length(goodASV_05_AR)==0){
        goodASV_05_AR<-goodASV_02_AR*0
      }
      
      p1_BA<-ggplot(comparison_02_BA, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")
      p2_BA<-ggplot(data.frame(value=goodASV_02_BA,phylum=names(goodASV_02_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
        ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")
      # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
      p3_BA<-ggplot(comparison_05_BA, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")
      p4_BA<-ggplot(data.frame(value=goodASV_05_BA,phylum=names(goodASV_05_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
        ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")
      # grid.arrange(p3_BA,p4_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_05"))
      p1_AR<-ggplot(comparison_02_AR, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")
      p2_AR<-ggplot(data.frame(value=goodASV_02_AR,phylum=names(goodASV_02_AR)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
        ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")
      # grid.arrange(p1_AR,p2_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_02"))
      p3_AR<-ggplot(comparison_05_AR, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")
      p4_AR<-ggplot(data.frame(value=goodASV_05_AR,phylum=names(goodASV_05_AR)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
        ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")
      # grid.arrange(p3_AR,p4_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_05"))
      
      
    }
    else{
      
      load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
      taxo_grouped_phylum<-eval(parse(text=paste0(group,"taxo_grouped_phylum")))
      # evalmetGLM<-evalmet
      # fitmetGLM<-fitmet
      
      #how many good ASV per phylum
      # summary(factor(PRtaxo_grouped_phylum[evalmetGLM$TSS>0.2&!(is.na(evalmetGLM$TSS)),]$Phylum))
      
      #how many very good ASV per phylum
      # summary(factor(PRtaxo_grouped_phylum[evalmetGLM$TSS>0.5&!(is.na(evalmetGLM$TSS)),]$Phylum))
      
      #how many modeled ASV per phylum
      # summary(factor(PRtaxo_grouped_phylum$Phylum))
      
      #proportion of each phylum in good modeled ASV pool
      goodASV_phylum_02<-summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.2&!(is.na(evalmet$Dpear)),]$Phylum))/sum(summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.2&!(is.na(evalmet$Dpear)),]$Phylum)))
      #proportion of each phylum in very good modeled ASV pool
      goodASV_phylum_05<-summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.5&!(is.na(evalmet$Dpear)),]$Phylum))/sum(summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.5&!(is.na(evalmet$Dpear)),]$Phylum)))
      
      #proportion of each phylum in all modeled ASV pool
      AllASV_phylum<-summary(factor(taxo_grouped_phylum$Phylum))/sum(summary(factor(taxo_grouped_phylum$Phylum)))
      
      comparison_02<-melt(rbind(AllASV_phylum,goodASV_phylum_02))
      comparison_05<-melt(rbind(AllASV_phylum,goodASV_phylum_05))
      
      #proportion of good modeled ASV in each phylum
      goodASV_02<-summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.2&!(is.na(evalmet$Dpear)),]$Phylum))/summary(factor(taxo_grouped_phylum$Phylum))
      
      #proportion of verygood modeled ASV in each phylum
      goodASV_05<-summary(factor(taxo_grouped_phylum[evalmet$Dpear>0.5&!(is.na(evalmet$Dpear)),]$Phylum))/summary(factor(taxo_grouped_phylum$Phylum))
      if(length(goodASV_05)==0){
        goodASV_05<-goodASV_02*0
      }
      
      assign(paste0("p1_", group), ggplot(comparison_02, aes(x=Var1, y=value, fill=Var2)) +
               geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")) 
      assign(paste0("p2_", group), ggplot(data.frame(value=goodASV_02,phylum=names(goodASV_02)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
               ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab(""))
      # grid.arrange(p1,p2,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
      assign(paste0("p3_", group), ggplot(comparison_05, aes(x=Var1, y=value, fill=Var2)) +
               geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab(""))
      assign(paste0("p4_", group), ggplot(data.frame(value=goodASV_05,phylum=names(goodASV_05)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
               ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab(""))
      # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
      
    }
    pdf(file=paste0("figures/test_phylum/goodmodels/",Mod,"_",PAAB,".pdf"))
    grid.arrange(p1_BA,p2_BA,p1_AR,p2_AR,p1_FU,p2_FU,p1_PR,p2_PR,nrow=4,top=paste0(Mod,"_",PAAB,"_02"))
    grid.arrange(p3_BA,p4_BA,p3_AR,p4_AR,p3_FU,p4_FU,p3_PR,p4_PR,nrow=4,top=paste0(Mod,"_",PAAB,"_05"))
    dev.off()
  }
  
}
