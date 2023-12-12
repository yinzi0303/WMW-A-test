

rm(list= ls())
gc()

library(monocle)
require(edgeR)
source("zingR.R")
require(DEsingle)

library(aggregateBioVar)
library(SummarizedExperiment, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(DESeq2, quietly = TRUE)

library(muscat)
source('Real_data/LP/Comp_code.R')

# ******************************************************************************************* #

# *** Read in the raw data:

dataset = readRDS('Real_data/LP/GSE54456.rds') 

data = dataset[["expr"]]
label = dataset$label



num = 100
n1=c(5,5,5,5,10,15)
n2=c(5,10,15,20,10,15)

Res = list()
sig_score = 0.01
result = list()


casesamples = which(label==1)
contsamples = which(label==2)

for (t in 1:num){
  
  
  # result$Res$ROTS[[t]] = list()
  # result$Res$DESeq[[t]] = list()
  result$Res$Limma[[t]] = list()
  # result$Res$MAST[[t]] = list()
  result$Res$monocle[[t]] = list()
  # result$Res$edgeR[[t]] = list()
  result$Res$wilx[[t]] = list()
  # result$Res$scde[[t]] = list()
  result$Res$ttest[[t]] = list()
  result$Res$ttestR[[t]] = list()
  result$Res$welch[[t]] = list()
  result$Res$zingeR[[t]] = list()
  result$Res$DEsingle[[t]] = list()
  # result$Res$Agg[[t]] = list()
  result$Res$muscat[[t]] = list()
  
  
  #length(n1)
  for (j in 1:length(n1)){
    
    sel_1 = sample(casesamples,n1[j])
    sel_2 = sample(contsamples,n2[j])
    sel = c(sel_1,sel_2)
    ex = data[,sel]
    
    sml <- c(rep(0,n1[j]),rep(1,n2[j]))
    sml = as.character(sml)
    gs <- factor(sml)
    
    
    ####################################################################################
    # Res_ROTS = Run_ROTS(rawData = ex, gr1= n1[j], gr2=n2[j], B=1000, K=3000 , FDR = sig_score )
    # Res_DESeq = Run_DESeq(rawData = ex, gs)
    Res_Limma = Run_Limma(rawData = ex, gs )
    # Res_MAST = Run_MAST(rawData = ex, gs)
    Res_monocle = Run_monocle(rawData = ex, gr1= n1[j], gr2=n2[j], gs )
    # Res_edgeR = Run_edgeR(rawData = ex, gr1= n1[j], gr2=n2[j], gs )
    Res_wilx = Run_wilx(rawData = ex, gr1= n1[j], gr2=n2[j])
    # Res_scde = Run_scde(rawData = ex, gs)
    Res_ttest = Run_ttest(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_ttestR = Run_ttestR(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_welch = Run_welch(rawData = ex, gr1= n1[j], gr2=n2[j])
    
    individual = colnames(ex)
    cluster = c(rep("spike",n1[j]+n2[j]))
    diagnosis = c(rep("healthy",n1[j]),rep("disease",n2[j]))
    metas = data.frame(individual=individual, diagnosis=diagnosis, cluster=cluster)
    Data <- data_process(Data = ex, group = sml, is.normalized = T)
    
    Res_zingeR <- Execute_zingeR.edgeR(Data,  maxit.EM = 200)
    Res_DEsingle <- Execute_DEsingle(object = Data, DEsingle.parallel = F)
    # Res_Agg = tryCatch(Run_Agg(rawdata = ex,metas = metas), error = function(e) {NA})
    Res_muscat = tryCatch(Run_muscat(rawdata = ex, metas = metas), error = function(e) {NA})
    
    ####################################################################################
    
    
    
    # result$Res$ROTS[[t]][[j]] = NULL
    # result$Res$DESeq[[t]][[j]] = Res_DESeq
    result$Res$Limma[[t]][[j]] = Res_Limma
    # result$Res$MAST[[t]][[j]] = Res_MAST
    result$Res$monocle[[t]][[j]] = Res_monocle
    # result$Res$edgeR[[t]][[j]] = Res_edgeR
    result$Res$wilx[[t]][[j]] = Res_wilx
    # result$Res$scde[[t]][[j]] = NULL
    result$Res$ttest[[t]][[j]] = Res_ttest
    result$Res$ttestR[[t]][[j]] = Res_ttestR
    result$Res$welch[[t]][[j]] = Res_welch
    result$Res$zingeR[[t]][[j]] = Res_zingeR
    result$Res$DEsingle[[t]][[j]] = Res_DEsingle
    result$Res$muscat[[t]][[j]] = Res_muscat
    

  }

  saveRDS(result,"Real_data/LP/compresult.rds")
  
}



