
source("generatexy.R")

geneNB1 = function(para){
  
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
      data = generatexy(Times = Times, m = m, n = n, alphax = alphax, alphay=alphay, beta = beta, distribution = "NB1")
      DX[[d]][[i]] = data$X
      DY[[d]][[i]] = data$Y
      DZ[[d]][[i]] = data$Z
      
    }
  }
  
  Data = list(X=DX, Y=DY, Z=DZ, M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
  
  return(Data)
}


geneNB2 = function(para){
  
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
      data = generatexy(Times = Times, m = m, n = n, alphax = alphax, alphay=alphay, beta = beta, distribution = "NB2")
      DX[[d]][[i]] = data$X
      DY[[d]][[i]] = data$Y
      DZ[[d]][[i]] = data$Z
    }
  }
  
  Data = list(X=DX, Y=DY, Z=DZ, M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
  
  return(Data)
}

##########################################################################################




##balanced case
M = c(2,4,5,8,10,12,15,18,20,25)
N = c(2,4,5,8,10,12,15,18,20,25)
alphax = 0.5
Alphay = c(0.4,0.42,0.45,0.48)
beta = 100
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB1(input_para)
saveRDS(Data, "data_nb1_eq.rds")


M = c(rep(5,10))
N = c(2,4,8,10,15,20,25,30,40,50)
alphax = 0.5
Alphay = c(0.4,0.42,0.45,0.48)
beta = 100
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB1(input_para)
saveRDS(Data, "data_nb1_ueq.rds")


#################################

##balanced case
M = c(2,4,5,8,10,12,15,18,20,25)
N = c(2,4,5,8,10,12,15,18,20,25)
alphax = 100
Alphay = c(75,80,90,95)
beta = 0.5
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB2(input_para)
saveRDS(Data, "data_nb2_eq.rds")


M = c(rep(5,10))
N = c(2,4,8,10,15,20,25,30,40,50)
alphax = 100
Alphay = c(75,80,90,95)
beta = 0.5
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB2(input_para)
saveRDS(Data, "data_nb2_ueq.rds")




##########################################################################################



##balanced case
M = c(2,4,5,8,10,12,15,18,20,25)
N = c(2,4,5,8,10,12,15,18,20,25)
alphax = 0.5
Alphay = c(0.5)
beta = 100
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB1(input_para)
saveRDS(Data, "data_nb1_eqt.rds")


M = c(rep(5,10))
N = c(2,4,8,10,15,20,25,30,40,50)
alphax = 0.5
Alphay = c(0.5)
beta = 100
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB1(input_para)
saveRDS(Data, "data_nb1_ueqt.rds")


#################################

##balanced case
M = c(2,4,5,8,10,12,15,18,20,25)
N = c(2,4,5,8,10,12,15,18,20,25)
alphax = 100
Alphay = c(100)
beta = 0.5
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB2(input_para)
saveRDS(Data, "data_nb2_eqt.rds")


M = c(rep(5,10))
N = c(2,4,8,10,15,20,25,30,40,50)
alphax = 100
Alphay = c(100)
beta = 0.5
Times = 1000
input_para = list(M=M, N=N, alphax=alphax, Alphay=Alphay, beta=beta, Times=Times)
Data = geneNB2(input_para)
saveRDS(Data, "data_nb2_ueqt.rds")


##########################################################################################
