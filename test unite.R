library(dada2)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
library(DECIPHER)
library(RSQLite)
# 
# refpath<-"Y:/FAC/FBM/DEE/aguisan/sometalp/D2c/SometAlp1.0/sh_general_release_dynamic_29.11.2022.fasta"
# unitedb<-read.table(refpath)

# load("../../ASV_data/ASV_taxo/FUtaxo.Rda")
load("FU_OTU_232samples.Rdata")

# fastaFung<-paste0(">Seq",rownames(FUtaxo),"\n",FUtaxo$Seq)
fastaFung<-paste0(">Seq",1:nrow(OTU_data_fun_232samples),"\n",rownames(OTU_data_fun_232samples))
# write(fastaFung,file="../../ASV_data/ASV_taxo/FastaFung.fasta")
write(fastaFung,file="../../ASV_data/ASV_taxo/FastaFung.fasta")


fas <- "FastaFung.fasta"
seqs <- readDNAStringSet(fas) 
seqs <- RemoveGaps(seqs)
load("../../ASV_data/ASV_taxo/UNITE_v2021_May2021.RData")
ids <- IdTaxa(seqs,
              trainingSet,
              strand="top", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
              threshold=60, # 60 (cautious) or 50 (sensible)
              processors=NULL) # use all available processors


save(ids,file="../../ASV_data/ASV_taxo/FU_assign_unite2021.Rda")
load("../../ASV_data/ASV_taxo/FU_assign_unite2021.Rda")

dbConn <- dbConnect(SQLite(),
                    ":memory:")

# print(ids[1:10])
# FUtaxo[FUtaxo$Phylum=="Zygomycota",-1]
# print(ids[FUtaxo$Phylum=="Zygomycota"])
# plot(ids,trainingSet) #super long
# plot(ids)
# plot(ids[FUtaxo$Phylum=="Zygomycota"],main="Zygomycota")
# plot(ids[FUtaxo$Phylum=="Basidiomycota"],main="Basidiomycota")
# plot(ids[FUtaxo$Phylum=="Ascomycota"],main="Ascomycota")
# plot(ids[FUtaxo$Phylum=="Fungi_unidentified"],main="Fungi_unidentified")
# plot(ids[FUtaxo$Phylum=="Glomeromycota"],main="Glomeromycota")
# plot(ids[FUtaxo$Phylum=="Chytridiomycota"],main="Chytridiomycota")
# plot(ids[FUtaxo$Phylum=="Neocallimastigomycota"],main="Neocallimastigomycota")
# print(ids[FUtaxo$Phylum=="Neocallimastigomycota"])

# ids[,c("rootrank","phylum")]#to look at certain taxo level
# ids[threshold=70] #to change threshold to higher confidence

assignment <- sapply(ids,
                     function(x)
                       paste(x$taxon,
                             collapse=";"))

kingdom <- sapply(ids,
                  function(x) {
                    w <- which(x$rank=="kingdom")
                    if (length(w) != 1) {
                      "unknown"
                    } else {
                      x$taxon[w]
                    }
                  })
# table(kingdom)

phylum <- sapply(ids,
                 function(x) {
                   w <- which(x$rank=="phylum")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
# table(phylum[kingdom =="Fungi"])

class <- sapply(ids,
                  function(x) {
                    w <- which(x$rank=="class")
                    if (length(w) != 1) {
                      "unknown"
                    } else {
                      x$taxon[w]
                    }
                  })
# table(class)

order <- sapply(ids,
                function(x) {
                  w <- which(x$rank=="order")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(order)

family <- sapply(ids,
                function(x) {
                  w <- which(x$rank=="family")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(family)

genus <- sapply(ids,
                 function(x) {
                   w <- which(x$rank=="genus")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
# table(genus)

species <- sapply(ids,
                function(x) {
                  w <- which(x$rank=="species")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(species)
# colnames(FUtaxo)
#check sequences in same order
# all(gsub("ASV","",names(kingdom)) == rownames(FUtaxo))
FUtaxo2021<-data.frame(Seq=FUtaxo$Seq,Kingdom=kingdom,Phylum=phylum,Class=class,Order=order,Family=family,Genus=genus,Species=species)
# save(FUtaxo2021,file="../../ASV_data/ASV_taxo/FUtaxo2021.Rda")
save(FUtaxo2021,file="FUtaxo2021.Rda")

q("no")