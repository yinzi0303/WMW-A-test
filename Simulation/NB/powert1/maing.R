
source("wmwa.R")


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
    samp=f[5]*rnorm(num, f[1], f[3]) + (1-f[5])*rnorm(num, f[2], f[4])
  }else{
    samp=f
  }
  
  return(samp)
}

z_generateNBM <- function(comb_xy, num, k){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitNB(comb_xy)
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  if (length(f)==7){
    samp=f[5]*rnbinom(num, size=1/(f[3]), mu=f[1]) + (1-f[5])*rnbinom(num, size=1/(f[4]), mu=f[2])
  }else{
    samp=f
  }
  
  return(samp)
}

z_generateNB <- function(comb_xy, num, k){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'negative binomial')
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  
  if (length(f)==num){
    samp=f
  }else{
    samp = rnbinom(num,size=f[["estimate"]][["size"]],mu=f[["estimate"]][["mu"]])
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
  
  wpnb = list()
  wpnbm = list()
  wpn = list()
  wpnm = list()
  
  empnb_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  empnbm_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  empn_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  empnm_power = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  for (d in 1:length(DX)){
    
    wpnb[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    wpnbm[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
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
      
      Z3 = apply(cbind(X, Y), 1, z_generateNB, num=l, k=k)
      Z3[which(is.na(Z3))]=0
      ALL3 <- cbind(cbind(X,Y), t(Z3))
      wpnb[[d]][,i] = apply(ALL3, 1, wmwa, m, n, N, nperm)
      empnb_power[d,i] <- length(which(wpnb[[d]][,i] <= 0.05)) / Times
      
      
      Z4 = apply(cbind(X, Y), 1, z_generateNBM, num=l, k=k)
      Z4[which(is.na(Z4))]=0
      ALL4 <- cbind(cbind(X,Y), t(Z4))
      wpnbm[[d]][,i] = apply(ALL4, 1, wmwa, m, n, N, nperm)
      empnbm_power[d,i] <- length(which(wpnbm[[d]][,i] <= 0.05)) / Times
      
    }
  }
  
  re1 = list(pval=list(wpnb, wpnbm, wpn, wpnm), power=list(empnb_power=empnb_power, 
                          empnbm_power=empnbm_power, empn_power=empn_power, empnm_power=empnm_power))
  
  return(re1)
  
}



# ##############################################################

data = readRDS("data_nb1_eqt.rds")
re_NB1_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_eq, "result_NB1_eqtg.rds")


data = readRDS("data_nb1_ueqt.rds")
re_NB1_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_ueq, "result_NB1_ueqtg.rds")


data = readRDS("data_nb2_eqt.rds")
re_NB2_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_eq, "result_NB2_eqtg.rds")


data = readRDS("data_nb2_ueqt.rds")
re_NB2_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_ueq, "result_NB2_ueqtg.rds")

# ##############################################################

data = readRDS("data_nb1_eq.rds")
re_NB1_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_eq, "result_NB1_eqg.rds")


data = readRDS("data_nb1_ueq.rds")
re_NB1_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB1_ueq, "result_NB1_ueqg.rds")


data = readRDS("data_nb2_eq.rds")
re_NB2_eq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_eq, "result_NB2_eqg.rds")


data = readRDS("data_nb2_ueq.rds")
re_NB2_ueq = wresult(data, k = 5, nperm = 1000)
saveRDS(re_NB2_ueq, "result_NB2_ueqg.rds")
