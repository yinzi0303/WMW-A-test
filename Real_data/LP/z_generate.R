
library("SIBER")
library("MASS")

z_generateNBM <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitNB(comb_xy)
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rnbinom(n, size=1/(f[3]), mu=f[1]) + (1-f[5])*rnbinom(n, size=1/(f[4]), mu=f[2])
  }else{
    samp=rep(comb_xy, 5)
  }
  
  return(samp)
}


# mixture normal distribution
z_generateNLM <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitNL(comb_xy)
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rnorm(n, f[1], f[3]) + (1-f[5])*rnorm(n, f[2], f[4])
  }else{
    samp=rep(comb_xy, 5)
  }
  
  return(samp)
}


# mixture log normal distribution
z_generateLNM <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitLN(comb_xy)
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rlnorm(n, f[1], f[3]) + (1-f[5])*rlnorm(n, f[2], f[4])
  }else{
    samp=rep(comb_xy, 5)
  }
  
  return(samp)
}



# normal distribution
z_generateNL <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'normal')
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  
  if (length(f)==n){
    samp=f
  }else{
    samp=rnorm(n, f[["estimate"]][["mean"]], f[["estimate"]][["sd"]]) 
  }
  
  return(samp)
}


# lognormal distribution
z_generateLN <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'lognormal')
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  
  if (length(f)==n){
    samp=f
  }else{
    samp=rlnorm(n, f[["estimate"]][["meanlog"]], f[["estimate"]][["sdlog"]]) 
  }
  
  return(samp)
}


z_generateNB <- function(comb_xy, n){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, 5)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'negative binomial')
      
    },error = function(e){
      rep(comb_xy, 5)
    })
  }
  
  
  if (length(f)==n){
    samp=f
  }else{
    samp = rnbinom(n,size=f[["estimate"]][["size"]],mu=f[["estimate"]][["mu"]])
  }
  
  return(samp)
}

