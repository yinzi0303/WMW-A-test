
source("Real_data/MS/WMWA.R")
source("Real_data/MS/z_generate.R")

ex_new = readRDS('Real_data/MS/ex_new.rds') #count data


sig_score = 0.01
result = list()
n1=c(5,5,5,5,10,15,20,25)
n2=c(5,10,15,20,10,15,20,25)


num= 100

for (t in 1:num){

  
  result$Res$WMWAN[[t]] = list()
  result$Res$WMWALN[[t]] = list()
  result$Res$WMWANM[[t]] = list()
  result$Res$WMWALNM[[t]] = list()
  result$Res$WMWANB[[t]] = list()
  result$Res$WMWANBM[[t]] = list()
  result$Res$WMWAreal[[t]] = list()
  
  
  

  #length(n1)
  for (j in 1:length(n1)){
    sel_1 = sample(1:623,n1[j],replace=FALSE)#sample(1:dim(casesamples)[2],n1[j])
    sel_2 = 623+sample(1:691,n2[j],replace=FALSE)#dim(casesamples)[2]+sample(1:dim(contsamples)[2],n2[j])
    sel = c(sel_1,sel_2)
    ex = ex_new[,sel]
    ex_xy = apply(ex, 2, as.numeric)
    sml <- c(rep(0,n1[j]),rep(1,n2[j]))
    sml = as.character(sml)
    gs <- factor(sml)
    
    n_xy = length(sel)
    n_z = 5 * n_xy
    
    ####################################################################################
    
    ex_z1 = apply(ex_xy, 1, z_generateNL, n=n_z)
    ex_z1[which(is.na(ex_z1))]=0
    ex = cbind(ex_xy, t(ex_z1))
    N = dim(ex)[2]
    
    Res_WMWAN = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWAN=as.vector(Res_WMWAN)
    result$Res$WMWAN[[t]][[j]] = Res_WMWAN
    
    ####################################################################################
    
    ex_xy_tmp = log10(ex_xy+1)
    ex_z2 = apply(ex_xy_tmp, 1, z_generateLN, n=n_z)
    ex_z2[which(is.na(ex_z2))]=0
    ex = cbind(ex_xy_tmp, t(ex_z2))
    # ex_z2 = apply(ex_xy, 1, z_generateLN, n=n_z)
    # ex_z2[which(is.na(ex_z2))]=0
    # ex = cbind(ex_xy, t(ex_z2))
    N = dim(ex)[2]
    
    Res_WMWALN = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWALN=as.vector(Res_WMWALN)
    result$Res$WMWALN[[t]][[j]] = Res_WMWALN
    
    ####################################################################################
    
    ex_z3 = apply(ex_xy, 1, z_generateNLM, n=n_z)
    ex_z3[which(is.na(ex_z3))]=0
    ex = cbind(ex_xy, t(ex_z3))
    N = dim(ex)[2]
    
    Res_WMWANM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANM=as.vector(Res_WMWANM)
    result$Res$WMWANM[[t]][[j]] = Res_WMWANM
    
    ####################################################################################
    
    ex_z4 = apply(ex_xy, 1, z_generateLNM, n=n_z)
    ex_z4[which(is.na(ex_z4))]=0
    ex = cbind(ex_xy, t(ex_z4))
    N = dim(ex)[2]
    
    Res_WMWALNM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWALNM=as.vector(Res_WMWALNM)
    result$Res$WMWALNM[[t]][[j]] = Res_WMWALNM
    
    
    
    ####################################################################################
    
    ex_z5 = apply(ex_xy, 1, z_generateNB, n=n_z)
    ex_z5[which(is.na(ex_z5))]=0
    ex = cbind(ex_xy, t(ex_z5))
    N = dim(ex)[2]
    
    Res_WMWANB = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANB = as.vector(Res_WMWANB)
    result$Res$WMWANB[[t]][[j]] = Res_WMWANB
    
    ####################################################################################
    
    
    ex_z6 = apply(ex_xy, 1, z_generateNBM, n=n_z)
    ex_z6[which(is.na(ex_z6))]=0
    ex = cbind(ex_xy, t(ex_z6))
    N = dim(ex)[2]
    
    Res_WMWANBM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANBM=as.vector(Res_WMWANBM)
    result$Res$WMWANBM[[t]][[j]] = Res_WMWANBM
    
    
    
    ####################################################################################
    left = setdiff(1:dim(ex_new)[2], sel)
    sel_z = sample(left, n_z, replace = FALSE)
    ex_z7 = ex_new[,sel_z]
    ex_z7[which(is.na(ex_z7))]=0
    ex = cbind(ex_xy, ex_z7)
    N = dim(ex)[2]
    
    Res_WMWAreal = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWAreal = as.vector(Res_WMWAreal)
    result$Res$WMWAreal[[t]][[j]] = Res_WMWAreal
    
    ####################################################################################
   
    
  
  }
  
  # DE=list(DE_gene=DE_gene,UNDE_gene=UNDE_gene)
  # saveRDS(DE,'degene9.rds')
  
  # result = list(rate=rate, auc=auc)
  saveRDS(result,"WMWAresult.rds")
  
}



