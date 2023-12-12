
source("wmwa.R")

wresult = function(data, k, nperm){
  
  DX = data$X
  DY = data$Y
  DZ = data$Z
  Times = data$Times
  
  wp = list()
  
  emp_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  for (d in 1:length(DX)){
    
    wp[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    
    for (i in 1:length(DX[[1]])){
      X = DX[[d]][[i]]
      Y = DY[[d]][[i]]
      Z = DZ[[d]][[i]]
      m = dim(X)[2]
      n = dim(Y)[2]
      l = k*(m+n)
      N = m+n+l
      
      ALL <- cbind(X, cbind(Y, Z))
      wp[[d]][,i] <- apply(ALL, 1, wmwa, m, n, N, nperm)
      
      emp_power[d,i] <- length(which(wp[[d]][,i] <= 0.05)) / Times
    }
  }
  
  re1 = list(pval=wp, power=emp_power)
  
  return(re1)
  
}


data = readRDS("data_nb1_eq.rds")
re_NB1_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_eq, "result_NB1_eq.rds")


data = readRDS("data_nb1_ueq.rds")
re_NB1_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_ueq, "result_NB1_ueq.rds")


data = readRDS("data_nb2_eq.rds")
re_NB2_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_eq, "result_NB2_eq.rds")


data = readRDS("data_nb2_ueq.rds")
re_NB2_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_ueq, "result_NB2_ueq.rds")
# 
# ##############################################################


data = readRDS("data_nb1_eqt.rds")
re_NB1_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_eq, "result_NB1_eqt.rds")


data = readRDS("data_nb1_ueqt.rds")
re_NB1_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_ueq, "result_NB1_ueqt.rds")


data = readRDS("data_nb2_eqt.rds")
re_NB2_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_eq, "result_NB2_eqt.rds")


data = readRDS("data_nb2_ueqt.rds")
re_NB2_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_ueq, "result_NB2_ueqt.rds")
