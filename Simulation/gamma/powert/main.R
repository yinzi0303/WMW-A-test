
source("Simulation/gamma/powert/wmwa.R")


require("SIBER")

# mixture normal distribution
z_generateNM <- function(comb_xy, num, k){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitNL(comb_xy)
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  if (length(f)==7){
    samp=f[5]*rnorm(num, f[1], f[3]) + (1-f[5])*rnorm(n, f[2], f[4])
  }else{
    samp=f
  }
  
  return(samp)
}



library("MASS")
# normal distribution
z_generateN <- function(comb_xy, num, k){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'normal')
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  
  if (length(f)==num){
    samp=f
  }else{
    samp=rnorm(num, f[["estimate"]][["mean"]], f[["estimate"]][["sd"]]) 
  }
  
  return(samp)
}







wresult = function(data, k, nperm){
  
  DX = data$X
  DY = data$Y
  DZ = data$Z
  Times = data$Times
  
  wp = list()
  wpn = list()
  wpnm = list()
  
  emp_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  empn_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  empnm_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  for (d in 1:length(DX)){
    
    wp[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    wpn[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    wpnm[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    
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
      
      
      Z1 = apply(cbind(X, Y), 1, z_generateN, num=l, k=k)
      Z1[which(is.na(Z1))]=0
      ALL1 <- cbind(cbind(X,Y), t(Z1))
      wpn[[d]][,i] = apply(ALL1, 1, wmwa, m, n, N, nperm)
      empn_power[d,i] <- length(which(wpn[[d]][,i] <= 0.05)) / Times
      
      
      Z2 = apply(cbind(X, Y), 1, z_generateNM, num=l, k=k)
      Z2[which(is.na(Z2))]=0
      ALL2 <- cbind(cbind(X,Y), t(Z2))
      wpnm[[d]][,i] = apply(ALL2, 1, wmwa, m, n, N, nperm)
      empnm_power[d,i] <- length(which(wpnm[[d]][,i] <= 0.05)) / Times
      
    }
  }
  
  re1 = list(pval=list(wp, wpn, wpnm), power=list(emp_power, empn_power, empnm_power))
  
  return(re1)
  
}



data = readRDS("Simulation/gamma/powert/data_gamma_eq.rds")
re_gamma_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_gamma_eq, "Simulation/gamma/powert/result_gamma_eq.rds")


data = readRDS("Simulation/gamma/powert/data_gamma_ueq.rds")
re_gamma_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_gamma_ueq, "Simulation/gamma/powert/result_gamma_ueq.rds")

