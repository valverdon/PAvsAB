library(dada2)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = '3.16')

library(DECIPHER)
library(RSQLite)

# load("../../ASV_data/ASV_taxo/UNITE_v2021_May2021.RData")
# load("../../ASV_data/ASV_taxo/FUtaxo2021.Rda")
# FUtaxo2021[FUtaxo2021$Kingdom=="Eukaryota_kgd_Incertae_sedis",]$Seq[1:10]
# 
# fastaFung<-paste0(">Seq",rownames(FUtaxo2021[FUtaxo2021$Kingdom=="Eukaryota_kgd_Incertae_sedis",][1:100,]),"\n",FUtaxo2021[FUtaxo2021$Kingdom=="Eukaryota_kgd_Incertae_sedis",][1:100,]$Seq)
# # fastaFung<-paste0(">Seq",1:nrow(OTU_data_fun_232samples),"\n",rownames(OTU_data_fun_232sam
# write(fastaFung,file="FastaFung_test10_unknown.fasta")

fas <- "FastaFung_test10_unknown.fasta"
seqs <- readDNAStringSet(fas) 
seqs <- RemoveGaps(seqs)

# ids <- IdTaxa(seqs,
#               trainingSet,
#               strand="both", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
#               threshold=30, # 60 (cautious) or 50 (sensible)
#               processors=NULL) # use all available processors
# 
# ids2 <- IdTaxa(seqs,
#               trainingSet,
#               strand="top", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
#               threshold=5, # 60 (cautious) or 50 (sensible)
#               processors=NULL) # use all available processors
# 
# plot(ids)
ids3<-assignTaxonomy(seqs,"Y:/FAC/FBM/DEE/aguisan/sometalp/D2c/Valentin_storage/unite/sh_general_release_dynamic_all_29.11.2022.fasta",outputBootstraps = TRUE,multithread =TRUE)
ids3[1:5]
unique(unname(ids3)[,1])
rownames(ids3)
