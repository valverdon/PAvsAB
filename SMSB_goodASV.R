#Get gud ASV to put in SMSB files
load("PA/BA/data/GLM/Eval.Rda")
load("PA/BA/data/GLM/Fit_Met.Rda")
evalmetGLM<-Eval
fitmetGLM<-fitmet
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo.Rda")
all(BAtaxo$Seq==colnames(OTUdata)) #check same sequences same order.
#number of prev >0.95
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Archaea"]],2,sum)/nrow(OTUdata))>0.95)
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Bacteria"]],2,sum)/nrow(OTUdata))>0.95)

evalmetGLM_BA<-evalmetGLM[BAtaxo$Kingdom=="Bacteria" & !(is.na(BAtaxo$Kingdom=="Bacteria")),]
evalmetGLM_BA_taxo<-cbind(evalmetGLM_BA,BAtaxo[BAtaxo$Kingdom=="Bacteria" & !(is.na(BAtaxo$Kingdom=="Bacteria")),])

listofgoodsp<-factor(evalmetGLM_BA_taxo$Species)[evalmetGLM_BA_taxo$TSS>0.5&!(is.na(evalmetGLM_BA_taxo$TSS))]
listofgoodsp<-listofgoodsp[!(is.na(listofgoodsp))]
listofgoodsp<-listofgoodsp[substr(listofgoodsp,1,10)!="uncultured"]
listofgoodsp<-listofgoodsp[-grep("sp.",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("Candida",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("metagenom",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("culture",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("bacterium",listofgoodsp)]
listofgoodsp<-listofgoodsp[-c(grep("unidentified",listofgoodsp),1)]
length(listofgoodsp)
listofgoodsp[grep("Wolbachia",listofgoodsp)]
unique(listofgoodsp)
table_goodsp<-evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species%in%listofgoodsp&!(is.na(evalmetGLM_BA_taxo$Species)),]
max(table_goodsp$TSS,na.rm=TRUE)
table_goodsp[table_goodsp$TSS==max(table_goodsp$TSS,na.rm=TRUE)& !(is.na(table_goodsp$TSS)),]
table_goodsp[table_goodsp$TSS==max(table_goodsp$TSS,na.rm=TRUE)& !(is.na(table_goodsp$TSS)),]$Species
evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Spirochaeta_aurantia"&!(is.na(evalmetGLM_BA_taxo$Species)),]
Spirochaeta_aurantia<-OTUdata[,evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Spirochaeta_aurantia"&!(is.na(evalmetGLM_BA_taxo$Species)),]$OTU]
Spirochaeta_aurantia<-apply(Spirochaeta_aurantia,1,any)
#Ecologie?
evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Pseudomonas_putida"&!(is.na(evalmetGLM_BA_taxo$Species)),]
Pseudomonas_putida<-OTUdata[,evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Pseudomonas_putida"&!(is.na(evalmetGLM_BA_taxo$Species)),]$OTU]
Pseudomonas_putida<-apply(Pseudomonas_putida,1,any)
sum(Pseudomonas_putida)

evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Pseudomonas_syringae_pv._broussonetiae"&!(is.na(evalmetGLM_BA_taxo$Species)),]
Pseudomonas_syringae<-OTUdata[,evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Pseudomonas_syringae_pv._broussonetiae"&!(is.na(evalmetGLM_BA_taxo$Species)),]$OTU]
sum(Pseudomonas_syringae)

evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Devosia_glacialis"&!(is.na(evalmetGLM_BA_taxo$Species)),]
Devosia_glacialis<-OTUdata[,evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Devosia_glacialis"&!(is.na(evalmetGLM_BA_taxo$Species)),]$OTU]
sum(Devosia_glacialis)

evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Nitrobacter_hamburgensis"&!(is.na(evalmetGLM_BA_taxo$Species)),]
Nitrobacter_hamburgensis<-OTUdata[,evalmetGLM_BA_taxo[evalmetGLM_BA_taxo$Species=="Nitrobacter_hamburgensis"&!(is.na(evalmetGLM_BA_taxo$Species)),]$OTU]
Nitrobacter_hamburgensis<-apply(Nitrobacter_hamburgensis,1,any)
sum(Nitrobacter_hamburgensis)

load("../../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")
dataSoil[dataSoil$sampleNameBA%in%names(Spirochaeta_aurantia),c("x","y")]
PAgoodBact<-cbind(Spirochaeta_aurantia,Pseudomonas_putida,Pseudomonas_syringae,Devosia_glacialis,Nitrobacter_hamburgensis,dataSoil[dataSoil$sampleNameBA%in%names(Spirochaeta_aurantia),c("x","y")])


load("PA/FU/data/GLM/Eval.Rda")
load("PA/FU/data/GLM/Fit_Met.Rda")
evalmetGLM<-Eval
fitmetGLM<-fitmet
load(paste0("PA/FU/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/FUtaxo.Rda")
all(FUtaxo$Seq==colnames(OTUdata)) #check same sequences same order.
#number of prev >0.95
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Archaea"]],2,sum)/nrow(OTUdata))>0.95)
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Bacteria"]],2,sum)/nrow(OTUdata))>0.95)

evalmetGLM_FU_taxo<-cbind(evalmetGLM,FUtaxo)
listofgoodsp[2837]
listofgoodsp<-factor(evalmetGLM_FU_taxo$Species)[evalmetGLM_FU_taxo$TSS>0.5&!(is.na(evalmetGLM_FU_taxo$TSS))]
listofgoodsp<-listofgoodsp[!(is.na(listofgoodsp))]

Cantharellus_appalachiensis_data<-evalmetGLM_FU_taxo[grep("Cantharellus_appalachiensis",evalmetGLM_FU_taxo$Species),]
Cantharellus_appalachiensis_seq<-as.vector(Cantharellus_appalachiensis_data$Seq)[which(Cantharellus_appalachiensis_data$TSS>0.65&!(is.na(Cantharellus_appalachiensis_data$TSS)))]
Cantharellus_appalachiensis<-OTUdata[,Cantharellus_appalachiensis_seq]
Cantharellus_appalachiensis<-apply(Cantharellus_appalachiensis,1,any)
sum(Cantharellus_appalachiensis)

evalmetGLM_FU_taxo[evalmetGLM_FU_taxo$Species%in%listofgoodsp[grep("Craterellus_tubaeformis",listofgoodsp)],]
Craterellus_tubaeformis_data<-evalmetGLM_FU_taxo[grep("Craterellus_tubaeformis",evalmetGLM_FU_taxo$Species),]
Craterellus_tubaeformis_seq<-as.vector(Craterellus_tubaeformis_data$Seq)
Craterellus_tubaeformis<-OTUdata[,Craterellus_tubaeformis_seq]
Craterellus_tubaeformis<-apply(Craterellus_tubaeformis,1,any)
sum(Craterellus_tubaeformis)

evalmetGLM_FU_taxo[evalmetGLM_FU_taxo$Species%in%listofgoodsp[grep("Lactarius_ilicis",listofgoodsp)],]
Lactarius_ilicis_data<-evalmetGLM_FU_taxo[grep("Lactarius_ilicis",evalmetGLM_FU_taxo$Species),]
Lactarius_ilicis_seq<-as.vector(Lactarius_ilicis_data$Seq)
Lactarius_ilicis<-OTUdata[,Lactarius_ilicis_seq]
sum(Lactarius_ilicis)

evalmetGLM_FU_taxo[evalmetGLM_FU_taxo$Species%in%listofgoodsp[grep("Morchella_populiphila|SH218667.06FU",listofgoodsp)],]
Morchella_populiphila_data<-evalmetGLM_FU_taxo[grep("Morchella_populiphila|SH218667.06FU",evalmetGLM_FU_taxo$Species),]
Morchella_populiphila_seq<-as.vector(Morchella_populiphila_data$Seq)
Morchella_populiphila<-OTUdata[,Morchella_populiphila_seq]
Morchella_populiphila<-apply(Morchella_populiphila,1,any)
sum(Morchella_populiphila)

dataSoil[dataSoil$sampleNameFU%in%names(Morchella_populiphila),c("x","y")]
PAgoodFung<-cbind(Cantharellus_appalachiensis,Craterellus_tubaeformis,Lactarius_ilicis,Morchella_populiphila,dataSoil[dataSoil$sampleNameFU%in%names(Morchella_populiphila),c("x","y")])



load("PA/PR/data/GLM/Eval.Rda")
load("PA/PR/data/GLM/Fit_Met.Rda")
evalmetGLM<-Eval
fitmetGLM<-fitmet
load(paste0("PA/PR/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/PRtaxo.Rda")
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.
which(PRtaxo$Seq=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
PRtaxo<-PRtaxo[PRtaxo$Seq!="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT",]
length(PRtaxo$Seq)
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.
#number of prev >0.95
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Archaea"]],2,sum)/nrow(OTUdata))>0.95)
# sum((apply(OTUdata[,colnames(OTUdata)%in%BAtaxo$Seq[BAtaxo$Kingdom=="Bacteria"]],2,sum)/nrow(OTUdata))>0.95)
PRspecieslist<-as.vector(apply(PRtaxo,1,function(X){X[sum(!(is.na(X)))]}))

evalmetGLM_PR_taxo<-cbind(evalmetGLM,PRtaxo)

listofgoodsp<-factor(PRspecieslist)[evalmetGLM_PR_taxo$TSS>0.5&!(is.na(evalmetGLM_PR_taxo$TSS))]
listofgoodsp<-listofgoodsp[!(is.na(listofgoodsp))]
listofgoodsp<-listofgoodsp[-grep("metagenom",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("culture",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("sp.",listofgoodsp)]
listofgoodsp<-listofgoodsp[-grep("environmental_sample",listofgoodsp)]


#oomycete
evalmetGLM_PR_taxo[PRspecieslist=="Achlya_hypogyna",]
Achlya_hypogyna_data<-evalmetGLM_PR_taxo[grep("Achlya_hypogyna",PRspecieslist),]
Achlya_hypogyna_seq<-as.vector(Achlya_hypogyna_data$Seq)
Achlya_hypogyna<-OTUdata[,Achlya_hypogyna_seq]
sum(Achlya_hypogyna)

#tested amoebae
evalmetGLM_PR_taxo[PRspecieslist=="Assulina_muscorum",]
Assulina_muscorum_data<-evalmetGLM_PR_taxo[grep("Assulina_muscorum",PRspecieslist),]
Assulina_muscorum_seq<-as.vector(Assulina_muscorum_data$Seq)
Assulina_muscorum<-OTUdata[,Assulina_muscorum_seq]
sum(Assulina_muscorum)

#Ciliophora
evalmetGLM_PR_taxo[PRspecieslist=="Colpoda_inflata",]
Colpoda_inflata_data<-evalmetGLM_PR_taxo[grep("Colpoda_inflata",PRspecieslist),]
Colpoda_inflata_seq<-as.vector(Colpoda_inflata_data$Seq)
Colpoda_inflata<-OTUdata[,Colpoda_inflata_seq]
Colpoda_inflata<-apply(Colpoda_inflata,1,any)
sum(Colpoda_inflata)

dataSoil[dataSoil$sampleNamePR%in%names(Colpoda_inflata),c("x","y")]
PAgoodProt<-cbind(Achlya_hypogyna,Assulina_muscorum,Colpoda_inflata,dataSoil[dataSoil$sampleNamePR%in%names(Achlya_hypogyna),c("x","y")])

save(PAgoodProt,file="PAgoodProt.Rda")
save(PAgoodFung,file="PAgoodFung.Rda")
save(PAgoodBact,file="PAgoodBact.Rda")
