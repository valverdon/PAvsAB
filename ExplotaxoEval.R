
load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
load("PA/PR/data/GLM/Eval.Rda")

if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGLM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(Eval)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGLM_PR_K<-evalmetGLM_PR[evalmetGLM_PR$Kingdom%in%names(summary(factor(evalmetGLM_PR$Kingdom)))[summary(factor(evalmetGLM_PR$Kingdom))>10],]
summary(factor(evalmetGLM_PR_K$Kingdom))
ggplot(evalmetGLM_PR_K, aes(x=Kingdom, y=TSS_adj, colour=Phylum))+ 
          geom_violin(aes(fill=Kingdom), scale="width",alpha=0.2)

evalmet_GLM_PR_Sup_P<-evalmetGLM_PR[evalmetGLM_PR$Sup_phylum%in%names(summary(factor(evalmetGLM_PR$Sup_phylum)))[summary(factor(evalmetGLM_PR$Sup_phylum))>10],]
ggplot(evalmet_GLM_PR_Sup_P, aes(x=Sup_phylum, y=TSS_adj, colour=Sup_phylum))+ 
  geom_violin(aes(fill=Phylum), scale="width",alpha=0.2)

evalmet_GLM_PR_P<-evalmetGLM_PR[evalmetGLM_PR$Phylum2%in%names(summary(factor(evalmetGLM_PR$Phylum2)))[summary(factor(evalmetGLM_PR$Phylum2))>10],]
ggplot(evalmet_GLM_PR_P, aes(x=Phylum2, y=TSS_adj, colour=Phylum2))+ 
  geom_violin(aes(fill=Phylum), scale="width",alpha=0.2)

evalmet_GLM_PR_C<-evalmetGLM_PR[evalmetGLM_PR$Class%in%names(summary(factor(evalmetGLM_PR$Class)))[summary(factor(evalmetGLM_PR$Class))>10],]
ggplot(evalmet_GLM_PR_C, aes(x=Class, y=TSS_adj, colour=Class))+ 
  geom_violin(aes(fill=Phylum), scale="width",alpha=0.2)

evalmetGLM_PR$Seq==evalmetGLM_PR$Seq
evalmet_GLM_PR_O<-evalmetGLM_PR[evalmetGLM_PR$Order%in%names(summary(factor(evalmetGLM_PR$Order)))[summary(factor(evalmetGLM_PR$Order))>10],]
summary(factor(evalmet_GLM_PR_O$Order))
ggplot(evalmet_GLM_PR_O, aes(x=Order, y=TSS_adj, colour=Order))+ 
  geom_violin(aes(fill=Phylum), scale="width",alpha=0.2) +
theme(axis.text.x = element_text(angle = 90))

ggplot(evalmetGLM_PR, aes(x=Phylum, y=TSS_adj, colour=Class))+ 
  geom_violin(aes(fill=Kingdom), scale="width",alpha=0.2)
