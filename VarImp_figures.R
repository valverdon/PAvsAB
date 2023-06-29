#figures Varimp

library(ggplot2)
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)
library(paletteer)

PAAB="AB"
GtoM="BA"
algo="GLM"
# load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
#colnames(ENVdata)
climatic<-c("bio1_t", "bio10_twa", "bio11_tco", "bio12_p" ,"bio13_pwe","bio14_pdr","bio15_ps" ,"bio16_pwe","bio17_pdr","bio18_pwa","bio19_pco",
            "bio2_tdr" ,"bio3_tiso","bio4_ts" , "bio5_tmax","bio6_tmin","bio7_tar", "bio8_twet","bio9_tdry","GDD0" , "ETP" ,  "cumday_no","sRadY"  )
edaphic<-c("pH","pH.1","bulkSoilW","soilTemp","EC_1_5",    "TotalP",    "Nitrogen",  "Carbon",    "Hydrogen",  "Phyllosil", "Quartz",    "Feldspath","Plagiocla","MassiveLi",
           "Calcite",   "Indoses",   "SiO2",      "TiO2",     "Al2O3",     "Fe2O3" ,"MnO"  , "MgO",   "CaO",   "Na2O","Marlyshal","MarlShale",
           "K2O" ,  "P2O5"  , "OM"  ,  "Cr2O3",  "NiO",   "d15N" ,  "d13C","Silt_clay", "clay"  , "ThinSilt", "ThickSilt", "ThinSand", "ThickSand","Soil_aera","Soil_humu","Soil_mois","Soil_mois.1","Soil_nutr" )
topographic<-c("aspect","slope","Elevation","Altitude")
landcover<-c("forest_ag", "hydro_agg", "lowVeg_ag", "anthropos",   "deciduous")
col_table<-data.frame(categ=c("climatic","edaphic","landcover","ndmi","ndvi","noise","topographic"),color=c("blue","brown","orange","darkgreen","lightgreen","purple","black"))

for (PAAB in c("PA","AB")){
  print(PAAB)
  for(algo in c("GAM","RF","GBM","GLM")){#Pb empty Variable_ranks.Rda GAM+GBM PA --> Bug Figure Selected Var
    print(algo)
    for(GtoM in c("PR","FU","BA")){
      print(GtoM)
      load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
      colnames(ENVdata)<-substr(colnames(ENVdata),1,9)
      # plot(hclust(as.dist(1-abs(cor(apply(ENVdata,2,as.numeric), use = "pairwise.complete.obs")))))
      # abline(h=0.3,col="red")
      
      
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_ranks.Rda"))
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_preselected.Rda"))
      
      names(Variable_importance)<-substr(names(Variable_preselected),1,9)
      names(Variable_ranks)<-substr(names(Variable_ranks),1,9)
      names(Variable_preselected)<-substr(names(Variable_preselected),1,9)
      
      
      # summary(Variable_importance)
      # summary(Variable_ranks)
      # apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],1,function(X){sum(!is.na(X))})
      
      # summary(Variable_preselected)
      if(GtoM=="BA"){
        load(paste0("../../ASV_data/ASV_taxo/",substr(GtoM,1,2),"taxo.Rda"))
        # all(Variable_importance$seq==BAtaxo$Seq)
        # all(Variable_ranks$seq==BAtaxo$Seq)
        # all(Variable_preselected$seq==BAtaxo$Seq)
        
        Variable_importance_BA <- Variable_importance[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
        Variable_importance_AR <- Variable_importance[BAtaxo$Kingdom=="Archaea"&!is.na(BAtaxo$Kingdom),]
        Variable_ranks_BA <- Variable_ranks[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
        Variable_ranks_AR <- Variable_ranks[BAtaxo$Kingdom=="Archaea"&!is.na(BAtaxo$Kingdom),]
        Variable_preselected_BA <- Variable_preselected[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
        Variable_preselected_AR <- Variable_preselected[BAtaxo$Kingdom=="Archaea"&!is.na(BAtaxo$Kingdom),]

        
              #BA
        
        allpresel<-apply(Variable_preselected_BA[,-which(colnames(Variable_preselected_BA)=="seq")],2,function(X){sum(X)/length(X)})
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]
        allranks<-apply(Variable_ranks_BA[,-which(colnames(Variable_ranks_BA)=="seq")],2,function(X){na.exclude(X)})
        Var_ranks<-lapply(allranks,as.vector)
        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,variable=factor(names(allpresel2),levels=names(allpresel2)),categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                                                                                                         ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                                                                                                ifelse(names(allpresel2)%in%topographic,"topographic",
                                                                            ifelse(names(allpresel2)%in%landcover,"landcover",names(allpresel2))))))
        categinplot_allpresel3<-unique(allpresel3$categ)   
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        # ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
        #   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
        assign(paste0("pPres_BA_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
                 scale_y_continuous(limits=c(0,1)))


      #pPres_BA_GAM_PA
        #per variable, proportion of models in which they where selected
        allsel_BA<-apply(Variable_ranks_BA[,-which(colnames(Variable_ranks)=="seq")],2,function(X){ifelse(is.na(X),NA,1)})
        allsel2<-apply(allsel_BA,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        allsel3<-data.frame(value=allsel2,variable=factor(names(allsel2),levels=names(allsel2)),categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                                                                                             ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                                                                                    ifelse(names(allsel2)%in%topographic,"topographic",
                                                                                                                           ifelse(names(allsel2)%in%landcover,"landcover",names(allsel2))))))
        categinplot_allsel3<-unique(allsel3$categ)   
        allsel3_10<-allsel3[1:10,]
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_BA_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
          scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
            scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_BA_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
      #pSel_FU_RF_PA
        
        
        #per variable, it's mean rank whithin postselected variables
        reordering<-with(melt(Var_ranks),                       # Order boxes by median
                         reorder(L1,
                                 value,
                                 median))
        
        data_reordered<-melt(Var_ranks)
        data_reordered$L1<-factor(data_reordered$L1,
                                  levels = levels(reordering))
        data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
                                     ifelse(X%in%edaphic,"edaphic",
                                            ifelse(X%in%topographic,"topographic",
                                                   ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_pRank<-unique(data_reordered$categ)   
        assign(paste0("pRank_BA_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank])+
                 scale_y_continuous(limits=c(0,1)))
        
      #pRank_BA_GLM_PA
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance2<-t(apply(Variable_importance_BA[,-which(colnames(Variable_importance_BA)=="seq")],1,function(X){X/max(X)}))
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        allsel4<-apply(Variable_importance2*allsel_BA,2,function(X){na.exclude(X)})
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)
        reordering4<-with(melt(allsel4),                       # Order boxes by median
                          reorder(L1,
                                  value,
                                  median,decreasing=TRUE))
        
        data_reordered4<-melt(allsel4)
        data_reordered4$L1<-factor(data_reordered4$L1,
                                   levels = levels(reordering4))
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                            ifelse(X%in%edaphic,"edaphic",
                                                                                                   ifelse(X%in%topographic,"topographic",
                                                                                                          ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)
        # data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        assign(paste0("pVIMP_BA_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_BA_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
      #pVIMP_BA_GLM_PA
        

                #AR
        
        allpresel<-apply(Variable_preselected_AR[,-which(colnames(Variable_preselected_AR)=="seq")],2,function(X){sum(X)/length(X)})
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]
        allranks<-apply(Variable_ranks_AR[,-which(colnames(Variable_ranks_AR)=="seq")],2,function(X){na.exclude(X)})
        Var_ranks<-lapply(allranks,as.vector)
        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,variable=factor(names(allpresel2),levels=names(allpresel2)),categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                                                                                                         ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                                                                                                ifelse(names(allpresel2)%in%topographic,"topographic",
                                                                                                                                       ifelse(names(allpresel2)%in%landcover,"landcover",names(allpresel2))))))
        categinplot_allpresel3<-unique(allpresel3$categ)   
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        # ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
        #   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
        assign(paste0("pPres_AR_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
                 scale_y_continuous(limits=c(0,1)))
        
      #pPres_AR_GLM_PA
        #per variable, proportion of models in which they where selected
        allsel_AR<-apply(Variable_ranks_AR[,-which(colnames(Variable_ranks)=="seq")],2,function(X){ifelse(is.na(X),NA,1)})
        allsel2<-apply(allsel_AR,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        allsel3<-data.frame(value=allsel2,variable=factor(names(allsel2),levels=names(allsel2)),categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                                                                                             ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                                                                                    ifelse(names(allsel2)%in%topographic,"topographic",
                                                                                                                           ifelse(names(allsel2)%in%landcover,"landcover",names(allsel2))))))
        categinplot_allsel3<-unique(allsel3$categ)   
        allsel3_10<-allsel3[1:10,]
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_AR_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_AR_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
        
        #per variable, it's mean rank whithin postselected variables
        reordering<-with(melt(Var_ranks),                       # Order boxes by median
                         reorder(L1,
                                 value,
                                 median))
        
        data_reordered<-melt(Var_ranks)
        data_reordered$L1<-factor(data_reordered$L1,
                                  levels = levels(reordering))
        data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                            ifelse(X%in%edaphic,"edaphic",
                                                                                                   ifelse(X%in%topographic,"topographic",
                                                                                                          ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_pRank<-unique(data_reordered$categ)   
        assign(paste0("pRank_AR_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank])+
                 scale_y_continuous(limits=c(0,1)))
        
        #pRank_AR_GLM_PA
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance2<-t(apply(Variable_importance_AR[,-which(colnames(Variable_importance_AR)=="seq")],1,function(X){X/max(X)}))
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        allsel4<-apply(Variable_importance2*allsel_AR,2,function(X){na.exclude(X)})
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)
        reordering4<-with(melt(allsel4),                       # Order boxes by median
                          reorder(L1,
                                  value,
                                  median,decreasing=TRUE))
        
        data_reordered4<-melt(allsel4)
        data_reordered4$L1<-factor(data_reordered4$L1,
                                   levels = levels(reordering4))
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                              ifelse(X%in%edaphic,"edaphic",
                                                                                                     ifelse(X%in%topographic,"topographic",
                                                                                                            ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)
        # data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        assign(paste0("pVIMP_AR_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_AR_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
        
      #pVIMP_AR_GLM_AB
        
        #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
        if(PAAB=="PA"){
          #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
            load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval.Rda"))
            Eval_BA<-Eval[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
            Eval_AR<-Eval[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
            TSSadj<-ifelse(Eval_BA$TSS_adj>0,Eval_BA$TSS_adj,NA)
            varimp_percentage<-t(apply(Variable_importance_BA[,-which(colnames(Variable_ranks_BA)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
            varimp_percentage<-varimp_percentage*allsel_BA
            varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})

            # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
            reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
                              reorder(Var2,
                                      value,
                                      median,decreasing=TRUE))
            
            data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
            data_reordered5$Var2<-factor(data_reordered5$Var2,
                                         levels = levels(reordering5))
            assign(paste0("pVIMPw_BA_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                     theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt")))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
            #pVIMPw_BA_GLM_PA
            
          TSSadj<-ifelse(Eval_AR$TSS_adj>0,Eval_AR$TSS_adj,NA)
          varimp_percentage<-t(apply(Variable_importance_AR[,-which(colnames(Variable_ranks_AR)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
          varimp_percentage<-varimp_percentage*allsel_AR
          varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})
          
          # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
          reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
                            reorder(Var2,
                                    value,
                                    median,decreasing=TRUE))
          
          data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
          data_reordered5$Var2<-factor(data_reordered5$Var2,
                                       levels = levels(reordering5))
          assign(paste0("pVIMPw_AR_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt")))
        #pVIMPw_AR
        }
      }
      if(GtoM %in% c("PR","FU")){
        #FU &PR
        allpresel<-apply(Variable_preselected[,-which(colnames(Variable_preselected)=="seq")],2,function(X){sum(X,na.rm=TRUE)/length(X)})
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]
        allranks<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){na.exclude(X)})
        Var_ranks<-lapply(allranks,as.vector)
        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,variable=factor(names(allpresel2),levels=names(allpresel2)),categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                                                                                                         ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                                                                                                ifelse(names(allpresel2)%in%topographic,"topographic",
                                                                                                                                       ifelse(names(allpresel2)%in%landcover,"landcover",names(allpresel2))))))
        categinplot_allpresel3<-unique(allpresel3$categ)   
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        assign(paste0("pPres_",GtoM,"_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
            scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
              scale_y_continuous(limits=c(0,1)))
              #pPres_PR_GLM_AB
        
        #per variable, proportion of models in which they where selected
        allsel<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){ifelse(is.na(X),NA,1)})
        allsel2<-apply(allsel,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        allsel3<-data.frame(value=allsel2,variable=factor(names(allsel2),levels=names(allsel2)),categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                                                                                             ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                                                                                    ifelse(names(allsel2)%in%topographic,"topographic",
                                                                                                                           ifelse(names(allsel2)%in%landcover,"landcover",names(allsel2))))))
        categinplot_allsel3<-unique(allsel3$categ)   
        allsel3_10<-allsel3[1:10,]
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_",GtoM,"_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_",GtoM,"_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("")+ ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
              #pSel_PR_GBM_PA
        
        
        #per variable, it's mean rank whithin postselected variables
        reordering<-with(melt(Var_ranks),                       # Order boxes by median
             reorder(L1,
                     value,
                     median))
        
        data_reordered<-melt(Var_ranks)
        data_reordered$L1<-factor(data_reordered$L1,
               levels = levels(reordering))
        data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                            ifelse(X%in%edaphic,"edaphic",
                                                                                                   ifelse(X%in%topographic,"topographic",
                                                                                                          ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_pRank<-unique(data_reordered$categ)   
        assign(paste0("pRank_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank])+
                 scale_y_continuous(limits=c(0,1)))
        # boxplot(value~testing,melt(Var_ranks),main=paste0(PAAB,"_",GtoM,"_",algo,"_Rank_Distribution_when_selected"),las = 2)
        # boxplot(value~L1,melt(Var_ranks),main=paste0(PAAB,"_",GtoM,"_",algo,"_Rank_Distribution_when_selected"),las = 2)
        
              #pRank_PR_GBM_PA
        #careful, possibility to have a low rank in a model with few variables
        
        
        
        
        #per variable, raw variable importance (unit depends on algo), 
        #all model performances weighted the same (true ecological importance depends on model performance)
        # boxplot(Variable_importance[,-which(colnames(Variable_importance)=="seq")],main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_all"),las = 2)
        
        # allsel3<-apply(Variable_importance[,-which(colnames(Variable_ranks)=="seq")]*allsel,2,function(X){na.exclude(X)})
        # boxplot(allsel3,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_only_when_selected"),las = 2)
        
        
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance2<-t(apply(Variable_importance[,-which(colnames(Variable_importance)=="seq")],1,function(X){X/max(X)}))
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        allsel4<-apply(Variable_importance2*allsel,2,function(X){na.exclude(X)})
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)
        reordering4<-with(melt(allsel4),                       # Order boxes by median
                         reorder(L1,
                                 value,
                                 median,decreasing=TRUE))
        
        data_reordered4<-melt(allsel4)
        data_reordered4$L1<-factor(data_reordered4$L1,
                                  levels = levels(reordering4))
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                              ifelse(X%in%edaphic,"edaphic",
                                                                                                     ifelse(X%in%topographic,"topographic",
                                                                                                            ifelse(X%in%landcover,"landcover",X))))}))
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)
        # data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        assign(paste0("pVIMP_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
                #pVIMP_FU_GLM_PA
                #pVIMP_FU_GLM_AB
        #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
        if(PAAB=="PA"){
          load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval.Rda"))
          TSSadj<-ifelse(Eval$TSS_adj>0,Eval$TSS_adj,NA)
          varimp_percentage<-t(apply(Variable_importance[,-which(colnames(Variable_ranks)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
          varimp_percentage<-varimp_percentage*allsel
          varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})
          
          # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
          reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
                            reorder(Var2,
                                    value,
                                    median,decreasing=TRUE))
          
          data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
          data_reordered5$Var2<-factor(data_reordered5$Var2,
                                     levels = levels(reordering5))
          assign(paste0("pVIMPw_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt")))
          #pVIMPw_FU_GLM_PA
        }
      }
    }
  }
}




# tPA <- textGrob("Presence-Absence models",hjust=-.2)
# tAB <- textGrob("Rel. Abundance models",hjust=-.1)

# plot(Ppresel_GLM)
Ppresel_GLM<-grid.arrange(grobs=list(pPres_BA_GLM_PA,pPres_BA_GLM_AB,pPres_AR_GLM_PA,pPres_AR_GLM_AB,pPres_FU_GLM_PA,pPres_FU_GLM_AB,pPres_PR_GLM_PA,pPres_PR_GLM_AB),
                          nrow=4)
Psel_GLM<-grid.arrange(pSel_BA_GLM_PA,pSel_BA_GLM_AB,pSel_AR_GLM_PA,pSel_AR_GLM_AB,pSel_FU_GLM_PA,pSel_FU_GLM_AB,pSel_PR_GLM_PA,pSel_PR_GLM_AB,nrow=4)
Psel10_GLM<-grid.arrange(pSel10_BA_GLM_PA,pSel10_BA_GLM_AB,pSel10_AR_GLM_PA,pSel10_AR_GLM_AB,pSel10_FU_GLM_PA,pSel10_FU_GLM_AB,pSel10_PR_GLM_PA,pSel10_PR_GLM_AB,nrow=4)
PRank_GLM<-grid.arrange(pRank_BA_GLM_PA,pRank_BA_GLM_AB,pRank_AR_GLM_PA,pRank_AR_GLM_AB,pRank_FU_GLM_PA,pRank_FU_GLM_AB,pRank_PR_GLM_PA,pRank_PR_GLM_AB,nrow=4)
PVIMP_GLM<-grid.arrange(pVIMP_BA_GLM_PA,pVIMP_BA_GLM_AB,pVIMP_AR_GLM_PA,pVIMP_AR_GLM_AB,pVIMP_FU_GLM_PA,pVIMP_FU_GLM_AB,pVIMP_PR_GLM_PA,pVIMP_PR_GLM_AB,nrow=4)
PVIMP10_GLM<-grid.arrange(pVIMP10_BA_GLM_PA,pVIMP10_BA_GLM_AB,pVIMP10_AR_GLM_PA,pVIMP10_AR_GLM_AB,pVIMP10_FU_GLM_PA,pVIMP10_FU_GLM_AB,pVIMP10_PR_GLM_PA,pVIMP10_PR_GLM_AB,
                          nrow=4,heights=c(2,2,2,2))
PVIMPw_GLM<-grid.arrange(pVIMPw_BA_GLM_PA,pVIMPw_AR_GLM_PA,pVIMPw_FU_GLM_PA,pVIMPw_PR_GLM_PA,nrow=4)

Ppresel_GAM<-grid.arrange(pPres_BA_GAM_PA,pPres_BA_GAM_AB,pPres_AR_GAM_PA,pPres_AR_GAM_AB,pPres_FU_GAM_PA,pPres_FU_GAM_AB,pPres_PR_GAM_PA,pPres_PR_GAM_AB,nrow=4)
Psel_GAM<-grid.arrange(pSel_BA_GAM_PA,pSel_BA_GAM_AB,pSel_AR_GAM_PA,pSel_AR_GAM_AB,pSel_FU_GAM_PA,pSel_FU_GAM_AB,pSel_PR_GAM_PA,pSel_PR_GAM_AB,nrow=4)
Psel10_GAM<-grid.arrange(pSel10_BA_GAM_PA,pSel10_BA_GAM_AB,pSel10_AR_GAM_PA,pSel10_AR_GAM_AB,pSel10_FU_GAM_PA,pSel10_FU_GAM_AB,pSel10_PR_GAM_PA,pSel10_PR_GAM_AB,nrow=4)
PRank_GAM<-grid.arrange(pRank_BA_GAM_PA,pRank_BA_GAM_AB,pRank_AR_GAM_PA,pRank_AR_GAM_AB,pRank_FU_GAM_PA,pRank_FU_GAM_AB,pRank_PR_GAM_PA,pRank_PR_GAM_AB,nrow=4)
PVIMP_GAM<-grid.arrange(pVIMP_BA_GAM_PA,pVIMP_BA_GAM_AB,pVIMP_AR_GAM_PA,pVIMP_AR_GAM_AB,pVIMP_FU_GAM_PA,pVIMP_FU_GAM_AB,pVIMP_PR_GAM_PA,pVIMP_PR_GAM_AB,nrow=4)
PVIMP10_GAM<-grid.arrange(pVIMP10_BA_GAM_PA,pVIMP10_BA_GAM_AB,pVIMP10_AR_GAM_PA,pVIMP10_AR_GAM_AB,pVIMP10_FU_GAM_PA,pVIMP10_FU_GAM_AB,pVIMP10_PR_GAM_PA,pVIMP10_PR_GAM_AB,nrow=4)
PVIMPw_GAM<-grid.arrange(pVIMPw_BA_GAM_PA,pVIMPw_AR_GAM_PA,pVIMPw_FU_GAM_PA,pVIMPw_PR_GAM_PA,nrow=4)

Ppresel_GBM<-grid.arrange(pPres_BA_GBM_PA,pPres_BA_GBM_AB,pPres_AR_GBM_PA,pPres_AR_GBM_AB,pPres_FU_GBM_PA,pPres_FU_GBM_AB,pPres_PR_GBM_PA,pPres_PR_GBM_AB,nrow=4)
Psel_GBM<-grid.arrange(pSel_BA_GBM_PA,pSel_BA_GBM_AB,pSel_AR_GBM_PA,pSel_AR_GBM_AB,pSel_FU_GBM_PA,pSel_FU_GBM_AB,pSel_PR_GBM_PA,pSel_PR_GBM_AB,nrow=4)
Psel10_GBM<-grid.arrange(pSel10_BA_GBM_PA,pSel10_BA_GBM_AB,pSel10_AR_GBM_PA,pSel10_AR_GBM_AB,pSel10_FU_GBM_PA,pSel10_FU_GBM_AB,pSel10_PR_GBM_PA,pSel10_PR_GBM_AB,nrow=4)
PRank_GBM<-grid.arrange(pRank_BA_GBM_PA,pRank_BA_GBM_AB,pRank_AR_GBM_PA,pRank_AR_GBM_AB,pRank_FU_GBM_PA,pRank_FU_GBM_AB,pRank_PR_GBM_PA,pRank_PR_GBM_AB,nrow=4)
PVIMP_GBM<-grid.arrange(pVIMP_BA_GBM_PA,pVIMP_BA_GBM_AB,pVIMP_AR_GBM_PA,pVIMP_AR_GBM_AB,pVIMP_FU_GBM_PA,pVIMP_FU_GBM_AB,pVIMP_PR_GBM_PA,pVIMP_PR_GBM_AB,nrow=4)
PVIMP10_GBM<-grid.arrange(pVIMP10_BA_GBM_PA,pVIMP10_BA_GBM_AB,pVIMP10_AR_GBM_PA,pVIMP10_AR_GBM_AB,pVIMP10_FU_GBM_PA,pVIMP10_FU_GBM_AB,pVIMP10_PR_GBM_PA,pVIMP10_PR_GBM_AB,nrow=4)
PVIMPw_GBM<-grid.arrange(pVIMPw_BA_GBM_PA,pVIMPw_AR_GBM_PA,pVIMPw_FU_GBM_PA,pVIMPw_PR_GBM_PA,nrow=4)

Ppresel_RF<-grid.arrange(pPres_BA_RF_PA,pPres_BA_RF_AB,pPres_AR_RF_PA,pPres_AR_RF_AB,pPres_FU_RF_PA,pPres_FU_RF_AB,pPres_PR_RF_PA,pPres_PR_RF_AB,nrow=4)
Psel_RF<-grid.arrange(pSel_BA_RF_PA,pSel_BA_RF_AB,pSel_AR_RF_PA,pSel_AR_RF_AB,pSel_FU_RF_PA,pSel_FU_RF_AB,pSel_PR_RF_PA,pSel_PR_RF_AB,nrow=4)
Psel10_RF<-grid.arrange(pSel10_BA_RF_PA,pSel10_BA_RF_AB,pSel10_AR_RF_PA,pSel10_AR_RF_AB,pSel10_FU_RF_PA,pSel10_FU_RF_AB,pSel10_PR_RF_PA,pSel10_PR_RF_AB,nrow=4)
PRank_RF<-grid.arrange(pRank_BA_RF_PA,pRank_BA_RF_AB,pRank_AR_RF_PA,pRank_AR_RF_AB,pRank_FU_RF_PA,pRank_FU_RF_AB,pRank_PR_RF_PA,pRank_PR_RF_AB,nrow=4)
PVIMP_RF<-grid.arrange(pVIMP_BA_RF_PA,pVIMP_BA_RF_AB,pVIMP_AR_RF_PA,pVIMP_AR_RF_AB,pVIMP_FU_RF_PA,pVIMP_FU_RF_AB,pVIMP_PR_RF_PA,pVIMP_PR_RF_AB,nrow=4)
PVIMP10_RF<-grid.arrange(pVIMP10_BA_RF_PA,pVIMP10_BA_RF_AB,pVIMP10_AR_RF_PA,pVIMP10_AR_RF_AB,pVIMP10_FU_RF_PA,pVIMP10_FU_RF_AB,pVIMP10_PR_RF_PA,pVIMP10_PR_RF_AB,nrow=4)
PVIMPw_RF<-grid.arrange(pVIMPw_BA_RF_PA,pVIMPw_AR_RF_PA,pVIMPw_FU_RF_PA,pVIMPw_PR_RF_PA,nrow=4)



pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Presel_GLMGAMGBMRF.pdf"))
plot(Ppresel_GLM)
plot(Ppresel_GAM)
plot(Ppresel_GBM)
plot(Ppresel_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Sel_GLMGAMGBMRF.pdf"))
plot(Psel_GLM)
plot(Psel_GAM)
plot(Psel_GBM)
plot(Psel_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Sel10_GLMGAMGBMRF.pdf"))
plot(Psel10_GLM)
plot(Psel10_GAM)
plot(Psel10_GBM)
plot(Psel10_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_RF)
dev.off()



pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Rank_GLMGAMGBMRF.pdf"))
plot(PRank_GLM)
plot(PRank_GAM)
plot(PRank_GBM)
plot(PRank_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Rank.png"),res=300,width=1961,height=1500)
plot(PRank_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Rank.png"),res=300,width=1961,height=1500)
plot(PRank_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Rank.png"),res=300,width=1961,height=1500)
plot(PRank_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Rank.png"),res=300,width=1961,height=1500)
plot(PRank_RF)
dev.off()


pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMP_GLMGAMGBMRF.pdf"))
plot(PVIMP_GLM)
plot(PVIMP_GAM)
plot(PVIMP_GBM)
plot(PVIMP_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMP10_GLMGAMGBMRF.pdf"))
plot(PVIMP10_GLM)
plot(PVIMP10_GAM)
plot(PVIMP10_GBM)
plot(PVIMP10_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_RF)
dev.off()


pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMPw_GLMGAMGBMRF.pdf"))
plot(PVIMPw_GLM)
plot(PVIMPw_GAM)
plot(PVIMPw_GBM)
plot(PVIMPw_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMPw.png"),res=300,width=1961,height=1500)
plot(PVIMPw_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMPw.png"),res=300,width=1961,height=1500)
plot(PVIMPw_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMPw.png"),res=300,width=1961,height=1500)
plot(PVIMPw_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMPw.png"),res=300,width=1961,height=1500)
plot(PVIMPw_RF)
dev.off()
