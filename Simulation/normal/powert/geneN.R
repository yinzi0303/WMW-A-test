setwd("C:/Users/dell/Desktop/wmwa/WMWA-R/normal")
source("generatexy.R")

geneN = function(para){
  
  M = para[["M"]]
  N = para[["N"]]
  alphax = para[["alphax"]]
  Alphay = para[["Alphay"]]
  beta = para[["beta"]]
  Times = para[["Times"]]
  
  DX = list()
  DY = list()
  DZ = list()
  
  for (d in 1:length(Alphay)){
    DX[[d]] = list()
    DY[[d]] = list()
    DZ[[d]] = list()
    
    alphay = Alphay[d]
    
    for (i in 1:length(N)){
      m = M[i]
      n = N[i]
      data = generatexy(Times = Times, m = m, n = n, alphax = alphax, alphay=alphay, 
                        beta = beta, distribution = "Normal")
      DX[[d]][[i]] = data$X
      DY[[d]][[i]] = data$Y
      DZ[[d]][[i]] = data$Z
      
    }
  }
  
  Data = list(X=DX, Y=DY, Z=DZ, M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
  
  return(Data)
}



##########################################################################################

setwd("C:/Users/dell/Desktop/wmwa/WMWA-R/normal")

##balanced case
M = c(2,4,5,8,10,12,15,18,20,25)
N = c(2,4,5,8,10,12,15,18,20,25)
alphax = 0
Alphay = c(0, 0.5, 1, 1.5, 2)
beta = 1
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneN(input_para)
saveRDS(Data, "data_normal_eq.rds")


M = c(rep(5,10))
N = c(2,4,8,10,15,20,25,30,40,50)
alphax = 0
Alphay = c(0, 0.5, 1, 1.5, 2)
beta = 1
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneN(input_para)
saveRDS(Data, "data_normal_ueq.rds")


##########################################################################################
