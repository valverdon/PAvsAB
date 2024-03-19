library(DHARMa)
library(glmnet)
####Dharma test over underdispersion zeroinf.



GtoM<-"BA"
Mod<-"GAM"
load(file=paste0("AB/",GtoM,"/Outputs/",Mod,"/Fit_data/Fit_data_temp41.Rda"))
load(file=paste0("AB/",GtoM,"/Outputs/",Mod,"/Models/Models_temp41.Rda"))
Modtotest<-Model_list
FittedVals<-Fit_data_list
# Modtotest<-readRDS(file=paste0("AB/",GtoM,"/Outputs/",Mod,"/Models/Models_temp41.Rds"))
# FittedVals<-readRDS(file=paste0("AB/",GtoM,"/Outputs/",Mod,"/Fit_data/Fit_data_temp41.Rds"))


simul_signif <- matrix(NA, nrow=length(Modtotest), ncol=4, dimnames = list(names(Modtotest),
                                                                          c("Uniformity", "OverDispersion", "UnderDispersion", "ZeroInflation")))
simul_values <- matrix(NA, nrow=length(Modtotest), ncol=4, dimnames = list(names(Modtotest),
                                                                          c("Uniformity", "OverDispersion", "UnderDispersion", "ZeroInflation")))

for (i in 1:length(names(Modtotest))){#i=names(Modtotest[1])
  if(!(is.na(unlist(FittedVals[[i]][2])[1]))){
    if(!(any(unlist(FittedVals[[i]][2])==Inf))){
      sims<-replicate(100,rpois(nrow(FittedVals[[i]][2]),unlist(FittedVals[[i]][2])))
      testDHARMA<-createDHARMa(simulatedResponse = sims, 
                               observedResponse = unlist(FittedVals[[i]][1]),
                               fittedPredictedResponse = unlist(FittedVals[[i]][2]), 
                               integer = TRUE)
      # plot(testDHARMA, quantreg = FALSE)
      
      
      simul_signif[names(Modtotest)[i], "Uniformity"] <- testUniformity(testDHARMA, plot=FALSE)$p.value # to test over/underdispersion #qqplot #Kolmogorov Smirnov test
      # plotQQunif(simulationOutput) # left plot in plot.DHARMa()
      # plotResiduals(simulationOutput) # right plot in plot.DHARMa()
      # testDispersion(simulationOutput)
      # hist(OTUdata[,otu], xlab = "Response", main = "")
      simul_values[names(Modtotest)[i], "Uniformity"] <- testUniformity(testDHARMA, plot=FALSE)$statistic # to test over/underdispersion
      simul_signif[names(Modtotest)[i], "OverDispersion"] <- testDispersion(testDHARMA, plot=FALSE, alternative = "greater")$p.value # to test overdispersion
      simul_values[names(Modtotest)[i], "OverDispersion"] <- testDispersion(testDHARMA, plot=FALSE, alternative = "greater")$statistic # to test overdispersion
      simul_signif[names(Modtotest)[i], "UnderDispersion"] <- testDispersion(testDHARMA, plot=FALSE, alternative = "less")$p.value # to test underdispersion
      simul_values[names(Modtotest)[i], "UnderDispersion"] <- testDispersion(testDHARMA, plot=FALSE, alternative = "less")$statistic # to test underdispersion
      simul_signif[names(Modtotest)[i], "ZeroInflation"] <- testZeroInflation(testDHARMA, plot=FALSE)$p.value
      simul_values[names(Modtotest)[i], "ZeroInflation"] <- testZeroInflation(testDHARMA, plot=FALSE)$statistic
    }
  }
}
