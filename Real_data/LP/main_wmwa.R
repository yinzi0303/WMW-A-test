
# install.packages("SIBER", repos="http://R-Forge.R-project.org")

source("Real_data/LP/WMWA.R")
source("Real_data/LP/z_generate.R")


dataset = readRDS('Real_data/LP/GSE54456.rds') 

data = dataset[["expr"]]
label = dataset$label


num = 100
Res = list()


sig_score = 0.01
result = list()

n1=c(5,5,5,5,10,15)
n2=c(5,10,15,20,10,15)


casesamples = which(label==1)
contsamples = which(label==2)


for (t in 1:num){
  
  
  result$Res$WMWANLM[[t]] = list()
  result$Res$WMWALNM[[t]] = list()
  result$Res$WMWANL[[t]] = list()
  result$Res$WMWALN[[t]] = list()
  result$Res$WMWAtrue[[t]] = list()
  
  for (j in 1:length(n1)){
    sel_1 = sample(casesamples,n1[j])
    sel_2 = sample(contsamples,n2[j])
    sel = c(sel_1,sel_2)
    ex_xy = data[,sel]
    ex_xy = apply(ex_xy, 2, as.numeric)
    ####################################################################################
    
    ex_z1 = apply(ex_xy, 1, z_generateNLM, n=(n1[j]+n2[j])*5)
    ex_z1[which(is.na(ex_z1))]=0
    ex = cbind(ex_xy, t(ex_z1))
    N = dim(ex)[2]
    
    Res_WMWANLM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANLM=as.vector(Res_WMWANLM)
    result$Res$WMWANLM[[t]][[j]] = Res_WMWANLM
    
    ####################################################################################
    
    ex_z2 = apply(ex_xy, 1, z_generateLNM, n=(n1[j]+n2[j])*5)
    ex_z2[which(is.na(ex_z2))]=0
    ex = cbind(ex_xy, t(ex_z2))
    N = dim(ex)[2]
    
    Res_WMWALNM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWALNM=as.vector(Res_WMWALNM)
    result$Res$WMWALNM[[t]][[j]] = Res_WMWALNM
    
    ####################################################################################
    
    ex_z3 = apply(ex_xy, 1, z_generateNL, n=(n1[j]+n2[j])*5)
    ex_z3[which(is.na(ex_z3))]=0
    ex = cbind(ex_xy, t(ex_z3))
    N = dim(ex)[2]
    
    Res_WMWANL = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANL=as.vector(Res_WMWANL)
    result$Res$WMWANL[[t]][[j]] = Res_WMWANL
    
    ####################################################################################
    
    ex_z4 = apply(ex_xy, 1, z_generateLN, n=(n1[j]+n2[j])*5)
    ex_z4[which(is.na(ex_z4))]=0
    ex = cbind(ex_xy, t(ex_z4))
    N = dim(ex)[2]
    
    Res_WMWALN = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWALN=as.vector(Res_WMWALN)
    result$Res$WMWALN[[t]][[j]] = Res_WMWALN
    
    ###################################################################################
    
    if (j<6){
      sel_left = setdiff(1:length(label), sel)
      sel_z = sample(sel_left, (n1[j]+n2[j])*5)
      ex_z5 = data[,sel_z]
      ex = cbind(ex_xy, ex_z5)
      N = dim(ex)[2]
      
      Res_WMWAtrue = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
      Res_WMWAtrue=as.vector(Res_WMWAtrue)
      result$Res$WMWAtrue[[t]][[j]] = Res_WMWAtrue
    }
    
    
    ####################################################################################
    
  }
  
  saveRDS(result,"Real_data/LP/WMWAresult.rds")
  
}

