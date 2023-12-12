
source("Real_data/MESC/WMWA.R")
source("Real_data/MESC/z_generate.R")


dataset = readRDS("Real_data/MESC/data_count_MESC.rds")
ex_new = dataset[['expr']]
label = dataset$label

casesamples = which(label==1)
contsamples = which(label==2)

num = 100
Res = list()
sig_score = 0.01
result = list()

n1=c(5,5,5,5,10,15,20,25)
n2=c(5,10,15,20,10,15,20,25)



for (t in 1:num){
  result$Res$WMWANB[[t]] = list()
  result$Res$WMWANNM[[t]] = list()
  result$Res$WMWA_count[[t]] = list()
  
  for (j in 1:length(n1)){
    sel_1 = sample(1:length(casesamples),n1[j])
    sel_2 = length(casesamples)+sample(1:length(contsamples),n2[j])
    sel = c(sel_1,sel_2)
    ex_xy = ex_new[,sel]
    sml <- c(rep(0,n1[j]),rep(1,n2[j]))
    sml = as.character(sml)
    gs <- factor(sml)
    
    ####################################################################################
    
    ex_z1 = apply(ex_xy, 1, z_generateNB, n=(n1[j]+n2[j])*3)
    ex_z1[which(is.na(ex_z1))]=0
    ex = cbind(ex_xy, t(ex_z1))
    N = dim(ex)[2]
    
    Res_WMWANB = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANB = as.vector(Res_WMWANB)
    result$Res$WMWANB[[t]][[j]] = Res_WMWANB
    
    ####################################################################################
    
    ex_z = apply(ex_xy, 1, z_generate, n=(n1[j]+n2[j])*3)
    ex_z[which(is.na(ex_z))]=0
    ex = cbind(ex_xy, t(ex_z))
    N = dim(ex)[2]
    
    Res_WMWANNM = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWANNM=as.vector(Res_WMWANNM)
    result$Res$WMWANNM[[t]][[j]] = Res_WMWANNM
    
    ####################################################################################
    n_z = (n1[j]+n2[j])*3
    left = setdiff(1:dim(ex_new)[2], sel)
    sel_z = sample(left, n_z, replace = FALSE)
    ex_z7 = ex_new[,sel_z]
    ex_z7[which(is.na(ex_z7))]=0
    ex = cbind(ex_xy, ex_z7)
    N = dim(ex)[2]
    
    Res_WMWA_count = apply(ex, 1, wmwa, m=n1[j], n=n2[j], N=N, nperm=1000)
    Res_WMWA_count=as.vector(Res_WMWA_count)
    result$Res$WMWA_count[[t]][[j]] = Res_WMWA_count
    
  }
  
  saveRDS(result,"WMWAresult_count.rds")
  
}

