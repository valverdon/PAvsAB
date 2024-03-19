
Model<-readRDS(file=paste0(savepath,"/",GtoM,"/Outputs/GLM/Models/Models_temp", arrayID, ".Rds"))
ENVstack<-readRDS(file="../../spatial_data/ENVstack.rds")
Presels<-readRDS(file=paste0(savepath,"/",GtoM,"/Outputs/GLM/VarSel/Varsel_temp", arrayID, ".Rds"))

Mod_to_proj<-Model[[1]]
predict.glmnet(Mod_to_proj)

##############################
########TODO##################18/03
##############################

form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar])),collapse=" + ")))
#glmnet needs to transform the function into a matrix format :
ModMat <- model.matrix(form,data = ENVdata)
