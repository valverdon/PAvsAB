###############################################################
# ENVdata
#Get topoclimatic data from different sources                            #
# extract values at sampling points                           #
# merge everything in Envdata_full                            #
#     #
###############################################################
#rm(list=ls())
library(terra)
library(tidyverse)
##################################### climatic and topographic variables ####################################################

load("../../spatial_data/ENVstack.Rda")
# names(ENVstack)
# plot(ENVstack[["Altitude"]])

#problem : 5 points out of Altitude layer

########################################## edaphic variables #####################################################################
#point data
# load OTU data from 265 sites in addition to soil data to find matching sites
load("PA\\BA\\data\\BA_OTU_265samples.Rdata")
# colnames(OTU_data_bac_265samples)
load("PA\\FU\\data\\FU_OTU_232samples.Rdata")
load("PA\\PR\\data\\PR_OTU_179samples.Rdata")

load("../../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")
edaph_data_BA <- dataSoil[dataSoil$sampleNameBA %in% colnames(OTU_data_bac_265samples), ] # pick only sites for which there is BA OTUdata
edaph_data_FU <- dataSoil[dataSoil$sampleNameFU %in% colnames(OTU_data_fun_232samples), ] # pick only sites for which there is FU OTUdata
edaph_data_PR <- dataSoil[dataSoil$sampleNamePR %in% colnames(OTU_data_pro_179samples), ] # pick only sites for which there is PR OTUdata
# edaph_data_BA$altitude
from_raster_BA <- extract(ENVstack, edaph_data_BA[,c("x","y")])
from_raster_FU <- extract(ENVstack, edaph_data_FU[,c("x","y")])
from_raster_PR <- extract(ENVstack, edaph_data_PR[,c("x","y")])

# plot(ENVstack$bio1_t)
# points(edaph_data_BA[which(edaph_data_BA$sampleNameBA%in%rownames(OTUdata_BA_PA)),c("x","y")],col="red")
# points(edaph_data_FU[which(edaph_data_FU$sampleNameFU%in%rownames(OTUdata_BA_FU)),c("x","y")],col="red")
# points(edaph_data_PR[which(edaph_data_PR$sampleNamePR%in%rownames(OTUdata_BA_PR)),c("x","y")],col="red")

# ENVdatatemp_BA <- data.frame(from_raster_BA[,colnames(from_raster_BA)!="d13C"], #for stack with d13C that double the var
#                              edaph_data_BA[,c(13:ncol(edaph_data_BA))],
#                              AllSand=edaph_data_BA$ThickSand+edaph_data_BA$ThinSand, 
#                              AllSilt=edaph_data_BA$ThickSilt+edaph_data_BA$ThinSilt)
# ENVdatatemp_FU <- data.frame(from_raster_FU[,colnames(from_raster_FU)!="d13C"], 
#                              edaph_data_FU[,c(13:ncol(edaph_data_FU))],
#                              AllSand=edaph_data_FU$ThickSand+edaph_data_FU$ThinSand, 
#                              AllSilt=edaph_data_FU$ThickSilt+edaph_data_FU$ThinSilt)
# ENVdatatemp_PR <- data.frame(from_raster_PR[,colnames(from_raster_PR)!="d13C"], 
#                              edaph_data_PR[,c(13:ncol(edaph_data_PR))],
#                              AllSand=edaph_data_PR$ThickSand+edaph_data_PR$ThinSand, 
                             # AllSilt=edaph_data_PR$ThickSilt+edaph_data_PR$ThinSilt)
ENVdatatemp_BA <- data.frame(from_raster_BA[,-which(colnames(from_raster_BA)%in%c("d13C"))],
                             edaph_data_BA[,c(13:ncol(edaph_data_BA))],
                             AllSand=edaph_data_BA$ThickSand+edaph_data_BA$ThinSand,
                             AllSilt=edaph_data_BA$ThickSilt+edaph_data_BA$ThinSilt)



# summary(ENVdatatemp_BA) ; nrow(ENVdatatemp_BA)
#remove some var 
#c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt")
#dry and wet because they are just here to have water content (multicol)
#EC because we have 2 ways of measuring the same thing
#AllSand and Allsilt because they are sum of more precise ones already in the table
#removed Ankerite, because "04_Clust_UnivMod.R" cant fit univariate gam model for OTU15660 for that variable. And its not an important variable.
#removed Dolomite, because not enough values != from 0 to perform variable selection with cubic regression.
#(later classified in "trash")
# varcor <- cor(ENVdatatemp_BA,use = "pairwise.complete.obs")
# plot(hclust(as.dist(1-abs(var_cor)), "average"))
# abline(h=0.3,lty=2, col="red", lwd=2)
# abline(h=0.2,lty=2, col="black", lwd=2)
# abline(h=0.25,lty=2, col="blue", lwd=2)

#Tmax removed, not interpretable and missing data
#C.N removed : Obvious multicolinearity (literally Carbon / Nitrogen) and missing data
#TOC removed : missing data and highly correlated with Nitrogen/Carbo /BSWC group
#MINC removed : missing data and highly correlated with CaO/pH group
#Phyllosilicates + Quartz Feldspath_K Plagioclase_Na    Calcite Goethite   Indoses : not very informative (except calcite but correlated), win more plots.
ENVdatatemp_BA <- ENVdatatemp_BA[,-which(names(ENVdatatemp_BA)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]

#remove HI OI cause too many empty data
ENVdatatemp_BA <- ENVdatatemp_BA[,-which(names(ENVdatatemp_BA)%in%c("HI","OI"))]

ENVdatatemp_BA <- ENVdatatemp_BA[,-which(names(ENVdatatemp_BA)%in%c("AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]

rownames(ENVdatatemp_BA) <- edaph_data_BA$sampleNameBA
#Remove sites with NA
#nrow(ENVdatatemp_BA)
# tail(ENVdatatemp_BA[!(complete.cases(ENVdatatemp_BA)),])
# summary(ENVdatatemp_BA[!(complete.cases(ENVdatatemp_BA)),])
ENVdatatemp_BA <- ENVdatatemp_BA[complete.cases(ENVdatatemp_BA),] #264 -> 250
# summary(ENVdatatemp_BA) ; nrow(ENVdatatemp_BA)

#remove binary var with not enough points with the var

#factorization of remaining binary var
binvar <- c()
for (i in names(ENVdatatemp_BA)){#i =names(ENVdata) [31]
  if(length(unique(ENVdatatemp_BA[,i]))==2){
    ENVdatatemp_BA[,i]<-as.factor(ENVdatatemp_BA[,i])
    binvar <- c(binvar,i)
  }
}
# summary(ENVdatatemp_BA)
# nrow(ENVdatatemp_BA) 250
# ncol(ENVdatatemp_BA) 77 dont 8 bin
ENVdata_BA <-ENVdatatemp_BA



ENVdatatemp_FU <- data.frame(from_raster_FU[,-which(colnames(from_raster_FU)%in%c("d13C"))], 
                             edaph_data_FU[,13:ncol(edaph_data_FU)], 
                             AllSand=edaph_data_FU$ThickSand+edaph_data_FU$ThinSand, 
                             AllSilt=edaph_data_FU$ThickSilt+edaph_data_FU$ThinSilt)

ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]

#remove HI OI cause too many empty data
ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("HI","OI"))]
rownames(ENVdatatemp_FU) <- edaph_data_FU$sampleNameFU
#Remove sites with NA

#nrow(ENVdata_PR)
ENVdatatemp_FU <- ENVdatatemp_FU[complete.cases(ENVdatatemp_FU),] #232 -> 217
# summary(ENVdatatemp_FU)
#nrow(ENVdatatemp_FU)

#remove binary var with not enough points with the var
ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]

#factorization of remaining binary var
binvar <- c()
for (i in names(ENVdatatemp_FU)){#i =names(ENVdata) [31]
  if(length(unique(ENVdatatemp_FU[,i]))==2){
    ENVdatatemp_FU[,i]<-as.factor(ENVdatatemp_FU[,i])
    binvar <- c(binvar,i)
  }
}
# summary(ENVdatatemp_FU)
# nrow(ENVdatatemp_FU) 217
# ncol(ENVdatatemp_FU) 77 dont 8 bin
ENVdata_FU <-ENVdatatemp_FU


ENVdatatemp_PR <- data.frame(from_raster_PR[,-which(colnames(from_raster_PR)%in%c("d13C"))], 
                             edaph_data_PR[,13:ncol(edaph_data_PR)], 
                             AllSand=edaph_data_PR$ThickSand+edaph_data_PR$ThinSand, 
                             AllSilt=edaph_data_PR$ThickSilt+edaph_data_PR$ThinSilt)

ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]

#remove HI OI cause too many empty data
ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("HI","OI"))]

rownames(ENVdatatemp_PR) <- edaph_data_PR$sampleNamePR
#Remove sites with NA

#nrow(ENVdata_PR)
ENVdatatemp_PR <- ENVdatatemp_PR[complete.cases(ENVdatatemp_PR),] #174 -> 166
# summary(ENVdatatemp_PR)
#nrow(ENVdatatemp_PR)

#remove binary var with not enough points with the var
ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]

#factorization of remaining binary var
binvar <- c()
for (i in names(ENVdatatemp_PR)){#i =names(ENVdata) [31]
  if(length(unique(ENVdatatemp_PR[,i]))==2){
    ENVdatatemp_PR[,i]<-as.factor(ENVdatatemp_PR[,i])
    binvar <- c(binvar,i)
  }
}
# summary(ENVdatatemp_PR)
# nrow(ENVdatatemp_PR) 166
# ncol(ENVdatatemp_PR) 77 dont 8 bin
ENVdata_PR <-ENVdatatemp_PR






### Only to have an idea about taxonomy in data filtering steps
load("../../../30_data/ASVtables_taxo/BAtaxo2019_full.Rda")
load("../../../30_data/ASVtables_taxo/PRtaxo2023_full.Rda")
load("../../../30_data/ASVtables_taxo/FUtaxo2023_full.Rda")

EUKtaxofull<-euk_silva.tax
PRtaxofull<-read.csv("../../../30_data/ASVtables_taxo/Tax_table_Protists.csv",sep=",")



#ASV data preparation
#Bacteria
OTUdata_BA_AB <- OTU_data_bac_265samples[,which(colnames(OTU_data_bac_265samples)%in%dataSoil$sampleNameBA)]

# load(paste0("../../ASV_data/ASV_taxo/BAtaxo.Rda"))
OTUdata_BA_AB_taxo<-cbind(OTUdata_BA_AB ,BAtaxo2019[match(rownames(OTUdata_BA_AB),BAtaxo2019[,"Seq"]),])
all(OTUdata_BA_AB_taxo$Seq==rownames(OTUdata_BA_AB_taxo))

OTUdata_BA_AB <- t(OTUdata_BA_AB)
# nrow(OTUdata_BA_AB) #264
# ncol(OTUdata_BA_AB) #60567 ASV
# sum(OTUdata_BA_AB_taxo$Kingdom=="Bacteria",na.rm = TRUE) 52177 Bacteria
# sum(OTUdata_BA_AB_taxo$Kingdom=="Archaea",na.rm = TRUE)  204 Archaea
# sum(OTUdata_BA_AB_taxo$Kingdom=="unclassified_Root",na.rm = TRUE) 8186

#length(complete.cases(ENVdatatemp_BA))
#remove plots that dont have measures for all variables




OTUdata_BA_AB <- OTUdata_BA_AB[rownames(ENVdata_BA), ] 
# nrow(OTUdata_BA_AB) #250

OTUdata_BA_AB <- OTUdata_BA_AB[,colSums(OTUdata_BA_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
#ncol(OTUdata_BA_AB) #58462
OTUdata_BA_AB_taxo2<-BAtaxo2019[match(colnames(OTUdata_BA_AB),BAtaxo2019[,"Seq"]),]
all(OTUdata_BA_AB_taxo2$Seq==colnames(OTUdata_BA_AB))
# sum(OTUdata_BA_AB_taxo2$Kingdom=="Bacteria",na.rm = TRUE) 50353 Bacteria
# sum(OTUdata_BA_AB_taxo2$Kingdom=="Archaea",na.rm = TRUE)  187 Archaea
# sum(OTUdata_BA_AB_taxo2$Kingdom=="unclassified_Root",na.rm = TRUE)  7922 ??

#easier datset for future (remove unassigned taxa directly instead of post-modelling):
OTUdata_BA_AB <- OTUdata_BA_AB[,OTUdata_BA_AB_taxo$Kingdom!="unclassified_Root"]
OTUdata_BA_AB_taxo2<-BAtaxo2019[match(colnames(OTUdata_BA_AB),BAtaxo2019[,"Seq"]),]

#Need treshold to remove very high and very low prevalence (algos wont be able to adjust to data)
#5 or 10% ?
numberpresence_BA <- apply(OTUdata_BA_AB,2,function(X){sum(X!=0)})
# plot(numberpresence_BA,pch=20,cex=.5, main="number of plots in which the ASV \n is present")
# # abline(h=nrow(OTUdata_BA_AB)/10,col="blue")#10% cut
# abline(h=nrow(OTUdata_BA_AB)/20,col="red") #5% cut
# # abline(h=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/10),col="blue")#10% cut
# abline(h=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),col="red") #5% cut
# plot(density(numberpresence_BA))

# # abline(v=nrow(OTUdata_BA_AB)/10,col="blue")#10% cut
# abline(v=nrow(OTUdata_BA_AB)/20,col="red")#5% cut
# # abline(v=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/10),col="blue")#10% cut
# abline(v=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),col="red") #5% cut
# densityx<-lapply(density(numberpresence_BA)$x,function(X){ifelse(X>nrow(OTUdata_BA_AB)/20,nrow(OTUdata_BA_AB)/20,X)})
# polygon(densityx,density(numberpresence_BA)$y,col="grey")
# densityx2<-lapply(density(numberpresence_BA)$x,function(X){ifelse(X<nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),X)})
# polygon(densityx2,density(numberpresence_BA)$y,col="grey")

removed_toomuchrare <- sum(numberpresence_BA<nrow(OTUdata_BA_AB)/20)
removed_toogeneral <- sum(numberpresence_BA>nrow(OTUdata_BA_AB)-nrow(OTUdata_BA_AB)/20)
removed_toorare_taxo<-OTUdata_BA_AB_taxo2[numberpresence_BA<nrow(OTUdata_BA_AB)/20,]
removed_toogeneral_taxo<-OTUdata_BA_AB_taxo2[numberpresence_BA>nrow(OTUdata_BA_AB)-nrow(OTUdata_BA_AB)/20,]

# removed_toomuchrare 2813
# removed_toogeneral 817   #they are only removed in PA models
# removed_toomuchrare+removed_toogeneral #3513
sum(removed_toorare_taxo$Kingdom=="Bacteria",na.rm = TRUE)#2037
sum(removed_toorare_taxo$Kingdom=="Archaea",na.rm = TRUE)#24
sum(removed_toorare_taxo$Kingdom=="unclassified_Root",na.rm = TRUE)# 752    if 0 -> removal ok
sum(removed_toogeneral_taxo$Kingdom=="Bacteria",na.rm = TRUE)#796
sum(removed_toogeneral_taxo$Kingdom=="Archaea",na.rm = TRUE)#0
sum(removed_toogeneral_taxo$Kingdom=="unclassified_Root",na.rm = TRUE)#21

#remove less prevalent ones for AB
OTUdata_BA_AB <- OTUdata_BA_AB[,numberpresence_BA>=nrow(OTUdata_BA_AB)/20]
OTUdata_BA_PA <- OTUdata_BA_AB[,(numberpresence_BA>=nrow(OTUdata_BA_AB)/20)*(numberpresence_BA<nrow(OTUdata_BA_AB)-nrow(OTUdata_BA_AB)/20)]

#double check
# load("PA/BA/data/OTUdata.Rda")
# all(colnames(OTUdata_BA_AB)==colnames(OTUdata))
# ncol(OTUdata_BA_AB)   #Bact+Arch +unassigned 55649
# ncol(OTUdata_BA_PA)# Bact+Arch +unassigned54832
#transformation AB --> PA
OTUdata_BA_PA <- apply(OTUdata_BA_AB,c(1,2),FUN=function(X){X>0}) 


seqvec <- colnames(OTUdata_BA_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_BA_AB))
# save(seqvec,file="seqvecBA.Rda")

load("../../ASV_data/ASV_taxo/BAtaxo.Rda")
all(BAtaxo$Seq==seqvec)
sum(BAtaxo$Kingdom=="Bacteria",na.rm=TRUE)
sum(BAtaxo$Kingdom=="Archaea",na.rm=TRUE)

# ncol(OTUdata_BA_PA) #55262 ASV (817 removed if prevalence <5%) 55262-817
# ncol(OTUdata_BA_AB) #55262 ASV
#253arrays of 220 ASV
TotSeqSum_BA <- rowSums(OTUdata_BA_AB)




#################################Fungi

OTUdata_FU_AB <- OTU_data_fun_232samples[,which(colnames(OTU_data_fun_232samples)%in%dataSoil$sampleNameFU)]
OTUdata_FU_AB_taxo<-cbind(OTUdata_FU_AB ,FUtaxo2021[match(rownames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),])
all(OTUdata_FU_AB_taxo$Seq==rownames(OTUdata_FU_AB_taxo))

OTUdata_FU_AB <- t(OTUdata_FU_AB)
# nrow(OTUdata_FU_AB) #232
#ncol(OTUdata_FU_AB) #95387
# sum(OTUdata_FU_AB_taxo$Kingdom=="Fungi",na.rm = TRUE) 31112 Fungi
# sum(OTUdata_FU_AB_taxo$Kingdom!="Fungi",na.rm = TRUE) 64261 not Fungi (or unclassified)
# sum(is.na(OTUdata_FU_AB_taxo$Kingdom)) 0

OTUdata_FU_AB <- OTUdata_FU_AB[rownames(ENVdata_FU), ] # #remove plots that dont have measures for all variables
# nrow(OTUdata_FU_AB) #217
OTUdata_FU_AB <- OTUdata_FU_AB[,colSums(OTUdata_FU_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them

# ncol(OTUdata_FU_AB) #92753
OTUdata_FU_AB_taxo2<-FUtaxo2021[match(colnames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),]
all(OTUdata_FU_AB_taxo2$Seq==colnames(OTUdata_FU_AB))

# sum(OTUdata_FU_AB_taxo2$Kingdom=="Fungi",na.rm = TRUE)# 30366 Fungi
# sum(OTUdata_FU_AB_taxo2$Kingdom!="Fungi",na.rm = TRUE)# 62387 not Fungi (or unclassified)
# sum(is.na(OTUdata_FU_AB_taxo2$Kingdom))# 0

#easier datset for future (remove unassigned taxa directly instead of post-modelling): #ignored for publication 2023
OTUdata_FU_AB <- OTUdata_FU_AB[,OTUdata_FU_AB_taxo2$Kingdom=="Fungi"]
OTUdata_FU_AB_taxo2<-FUtaxo2021[match(colnames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),]
#

numberpresence_FU <- apply(OTUdata_FU_AB,2,function(X){sum(X!=0)})
# plot(numberpresence_FU,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
# abline(h=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut 21
# abline(h=nrow(OTUdata_FU_AB)/20,col="red") #5% cut 10
# # abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# # abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20),col="red") #5% cut
# plot(density(numberpresence_FU))
# abline(v=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut
# abline(v=nrow(OTUdata_FU_AB)/20,col="red")#5% cut
# #abline(v=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# #abline(v=nrow(OTUdata_FU_AB)-(nrow(OTUdata_FU_AB)/20),col="red") #5% cut



removed_toorare <- sum(numberpresence_FU<nrow(OTUdata_FU_AB)/20)
removed_toogeneral <- sum(numberpresence_FU>nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20)
removed_toorare_taxo<-OTUdata_FU_AB_taxo2[numberpresence_FU<nrow(OTUdata_FU_AB)/20,]
removed_toogeneral_taxo<-OTUdata_FU_AB_taxo2[numberpresence_FU>nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20,]

# removed_toorare #40890
# removed_toogeneral #34  #they are only removed in PA models
# removed_toorare+removed_toogeneral #40924
sum(removed_toorare_taxo$Kingdom=="Fungi",na.rm = TRUE)#13021
sum(removed_toorare_taxo$Kingdom!="Fungi",na.rm = TRUE)#27869
sum(is.na(removed_toorare_taxo$Kingdom!="Fungi"))

sum(removed_toogeneral_taxo$Kingdom=="Fungi",na.rm = TRUE)#27
sum(removed_toogeneral_taxo$Kingdom!="Fungi",na.rm = TRUE)#7

# OTUdata_FU_PA <- OTUdata_FU_PA[,numberpresence_FU_PA>=nrow(OTUdata_FU_PA)/20]

#remove less prevalent ones for AB
OTUdata_FU_AB <- OTUdata_FU_AB[,numberpresence_FU>=nrow(OTUdata_FU_AB)/20]
OTUdata_FU_PA <- OTUdata_FU_AB[,(numberpresence_FU>=nrow(OTUdata_FU_AB)/20)*(numberpresence_FU<nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20)]

# ncol(OTUdata_FU_AB)# 51863
# ncol(OTUdata_FU_PA)# 51829

OTUdata_FU_PA <- apply(OTUdata_FU_AB,c(1,2),FUN=function(X){X>0})


seqvec <- colnames(OTUdata_FU_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_FU_AB))
# save(seqvec,file="seqvecFU.Rda")

# ncol(OTUdata_FU_AB) #51863
# ncol(OTUdata_FU_PA) #51863

TotSeqSum_FU <- rowSums(OTUdata_FU_AB)




#Protists
OTUdata_PR_AB <- OTU_data_pro_179samples[,which(colnames(OTU_data_pro_179samples)%in%dataSoil$sampleNamePR)]


OTUdata_PR_AB_taxo<-cbind(OTUdata_PR_AB ,PRtaxo2023[match(rownames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),])
all(OTUdata_PR_AB_taxo$Seq==rownames(OTUdata_PR_AB_taxo))
OTUdata_PR_AB <- t(OTUdata_PR_AB)

# nrow(OTUdata_PR_AB) #174
# ncol(OTUdata_PR_AB) #11239
# sum(OTUdata_PR_AB_taxo$Phylum!="Not_Protist",na.rm = TRUE) 2331 Protists
# sum(OTUdata_PR_AB_taxo$Phylum=="Not_Protist",na.rm = TRUE) 8908 not Fungi (or unclassified)
# sum(is.na(OTUdata_PR_AB_taxo$Phylum)) 0

# TotSeqSum_PR_AB <- rowSums(OTUdata_PR) #sequencing depth for "precleaned" data
OTUdata_PR_AB <- OTUdata_PR_AB[rownames(ENVdata_PR), ] # #remove plots that dont have measures for all variables

# nrow(OTUdata_PR_AB) #166
OTUdata_PR_AB <- OTUdata_PR_AB[,colSums(OTUdata_PR_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
# ncol(OTUdata_PR_AB) #10860
OTUdata_PR_AB_taxo2<-PRtaxo2023[match(colnames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),]
all(OTUdata_PR_AB_taxo2$Seq==colnames(OTUdata_PR_AB))

# sum(OTUdata_PR_AB_taxo2$Phylum!="Not_Protist",na.rm = TRUE) 2270 Protists
# sum(OTUdata_PR_AB_taxo2$Phylum=="Not_Protist",na.rm = TRUE) 8590 not Fungi (or unclassified)
# sum(is.na(OTUdata_PR_AB_taxo$Phylum)) 0

# GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT
#Remove non protists
# troph_euk<-read.csv("../../ASV_data/Troph_Prot.csv",sep=";")
# troph_euk<-read.csv("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/Tax_table_Protists.csv",sep=",")

# OTUdata_PR_AB <- OTUdata_PR_AB[,which(colnames(OTUdata_PR_AB)%in%troph_euk$X)]


#easier datset for future (remove unassigned taxa directly instead of post-modelling): #ignored for publication 2023
OTUdata_PR_AB <- OTUdata_PR_AB[,OTUdata_PR_AB_taxo2$Phylum!="Not_Protist"]
OTUdata_PR_AB_taxo2<-PRtaxo2023[match(colnames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),]
#
# ncol(OTUdata_PR_AB ) #2270

numberpresence_PR <- apply(OTUdata_PR_AB,2,function(X){sum(X!=0)})
# plot(numberpresence_PR,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
# abline(h=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
# abline(h=nrow(OTUdata_PR_AB)/20,col="red") #5% cut
# # abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# # abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut
# plot(density(numberpresence_PR))
# abline(v=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
# abline(v=nrow(OTUdata_PR_AB)/20,col="red")#5% cut
# # abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# # abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut

removed_toorare <- sum(numberpresence_PR<nrow(OTUdata_PR_AB)/20)
removed_toogeneral <- sum(numberpresence_PR>nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20)
removed_toorare_taxo<-OTUdata_PR_AB_taxo2[numberpresence_PR<nrow(OTUdata_PR_AB)/20,]
removed_toogeneral_taxo<-OTUdata_PR_AB_taxo2[numberpresence_PR>nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20,]
# removed_toorare 108
# removed_toogeneral 8  #they are only removed in PA models
# removed_toorare+removed_toogeneral #116
sum(removed_toorare_taxo$Phylum!="Not_Protist",na.rm = TRUE)#108
sum(removed_toorare_taxo$Phylum=="Not_Protist",na.rm = TRUE)#0
sum(is.na(removed_toorare_taxo$Phylum))

sum(removed_toogeneral_taxo$Phylum!="Not_Protist",na.rm = TRUE)#27
sum(removed_toogeneral_taxo$Phylum=="Not_Protist",na.rm = TRUE)#7


#remove less prevalent ones for AB
OTUdata_PR_AB <- OTUdata_PR_AB[,numberpresence_PR>=nrow(OTUdata_PR_AB)/20]
OTUdata_PR_PA <- OTUdata_PR_AB[,(numberpresence_PR>=nrow(OTUdata_PR_AB)/20)*(numberpresence_PR<nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20)]

# ncol(OTUdata_PR_AB)# 2162
# ncol(OTUdata_PR_PA)# 2154

OTUdata_PR_PA <- apply(OTUdata_PR_AB,c(1,2),FUN=function(X){X>0})
# ncol(OTUdata_PR_AB) #2162
# which(colnames(OTUdata_PR_AB)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
#annoying  Champi = 3019




# ncol(OTUdata_PR_PA) #3160

seqvec <- colnames(OTUdata_PR_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_PR_AB))
# save(seqvec,file="seqvecPR.Rda")

#cut at 10% : BA = 26, FU = 23, PR = 17 
#cut at 5% : BA = 13 FU = 11 PR = 8
#Where to cut here? 
#10% removes most of fungi (most of them are present in a very few plots)
#5% leaves some protists with only 9 plots to fit the models
TotSeqSum_PR <- rowSums(OTUdata_PR_AB)

# datapoints<-runif(100,10,40)
# dporder<-datapoints[order(datapoints)]
# 
# data<-matrix(data=c(datapoints,))
# plot(dporder,dnorm(dporder,25,3)*7.5,type="l",lwd=5,ylim=c(0,1))
# abline(h=0)
# abline(h=1)
# points(dporder,rbinom(100,1,dnorm(dporder,25,3)*7.5),pch=20,col="blue",cex=2)

colnames(ENVdata_BA)
#save all as "OTUdata"
OTUdata <- OTUdata_BA_PA
ENVdata <- ENVdata_BA
# save(OTUdata,file="PA/BA/Outputs/OTUdata.Rda")

# save(OTUdata,file="PA/BA/data/OTUdata.Rda")
# save(ENVdata,file="PA/BA/Outputs/ENVdata.Rda")


OTUdata <- OTUdata_BA_AB
TotSeqSum <- TotSeqSum_BA
# save(OTUdata,file="AB/BA/Outputs/OTUdata.Rda")
# save(ENVdata,file="AB/BA/Outputs/ENVdata.Rda")
# save(TotSeqSum,file="AB/BA/Outputs/dataTotSeqSum.Rda")


OTUdata <- OTUdata_FU_PA
ENVdata <- ENVdata_FU
# save(OTUdata,file="PA/FU/Outputs/OTUdata.Rda")
# save(ENVdata,file="PA/FU/Outputs/ENVdata.Rda")
# save(ENVdata,file="PA/FU/data/ENVdata.Rda")

OTUdata <- OTUdata_FU_AB
TotSeqSum <- TotSeqSum_FU
# save(OTUdata,file="AB//FU/Outputs/OTUdata.Rda")
# save(ENVdata,file="AB//FU/Outputs/ENVdata.Rda")
# save(TotSeqSum,file="AB//FU/Outputs/dataTotSeqSum.Rda")

OTUdata <- OTUdata_PR_PA
ENVdata <- ENVdata_PR
# save(OTUdata,file="PA/PR/Outputs/OTUdata.Rda")
# save(ENVdata,file="PA/PR/Outputs/ENVdata.Rda")
# save(ENVdata,file="PA/PR/data/ENVdata.Rda")
OTUdata <- OTUdata_PR_AB
TotSeqSum <- TotSeqSum_PR
# save(OTUdata,file="AB//PR/Outputs/OTUdata.Rda")
# save(ENVdata,file="AB//PR/Outputs/ENVdata.Rda")
# save(TotSeqSum,file="AB//PR/Outputs/dataTotSeqSum.Rda")

# write.csv(ENVdata,file="ENVdata_FU")

#plot points
hillshade <- crop(readRDS("../../spatial_data/Valpar/ch_topo_alti3d2016_pixel_hillshade_mean.rds"),extent(c(xmin=2552000,xmax=2587000,ymin=1114000,ymax=1157000)))

plot(ENVstack$aspect)
hillshade<-mask(hillshade,ENVstack$aspect)
plot(hillshade,col=gray.colors(1000))

points(edaph_data_BA[,c("x","y")],pch=16,col="#aec800",cex=1.2)
points(edaph_data_FU[,c("x","y")],pch=17,col=alpha("#6B4C62",0.8))
points(edaph_data_PR[,c("x","y")],pch=18,col=alpha("#457EB0",0.7))
legend(2550000,1130000,c("Bacteria\n & Archaea (16S)","Fungi (ITS)", "Protista (18S)"),pch=c(16,17,18), col=c(alpha("#aec800",0.8),alpha("#6B4C62",0.8),"#457EB0"),bty="n")


pdf(file="figures/Datasets_sampling.pdf")
plot(hillshade,col=gray.colors(1000))

points(edaph_data_BA[,c("x","y")],pch=16,col="#aec800",cex=1.2)
points(edaph_data_FU[,c("x","y")],pch=17,col=alpha("#6B4C62",0.8))
points(edaph_data_PR[,c("x","y")],pch=18,col=alpha("#457EB0",0.7))
legend(2540000,1130000,c("Bacteria & \nArchaea (16S)","Fungi (ITS)", "Protista (18S)"),pch=c(16,17,18), col=c(alpha("#aec800",0.8),alpha("#6B4C62",0.8),"#457EB0"),bty="n")
dev.off()
png(file=paste0("figures/Datasets_sampling.png"),res=300,width=1961,height=1500)
plot(hillshade,col=gray.colors(1000))
points(edaph_data_BA[,c("x","y")],pch=16,col="#aec800",cex=1.2)
points(edaph_data_FU[,c("x","y")],pch=17,col=alpha("#6B4C62",0.8))
points(edaph_data_PR[,c("x","y")],pch=18,col=alpha("#457EB0",0.7))
legend(2540000,1130000,c("Bacteria & \nArchaea (16S)","Fungi (ITS)", "Protista (18S)"),pch=c(16,17,18), col=c(alpha("#aec800",0.8),alpha("#6B4C62",0.8),"#457EB0"),bty="n")
dev.off()
