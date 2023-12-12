
library("SIBER")
library("MASS")

z_generateNBM <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitNB(comb_xy)
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rnbinom(num_l, size=1/(f[3]), mu=f[1]) + (1-f[5])*rnbinom(num_l, size=1/(f[4]), mu=f[2])
  }else{
    samp=rep(comb_xy, k)
  }
  
  return(samp)
}


# mixture normal distribution
z_generateNLM <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitNL(comb_xy)
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rnorm(num_l, f[1], f[3]) + (1-f[5])*rnorm(num_l, f[2], f[4])
  }else{
    samp=rep(comb_xy, k)
  }
  
  return(samp)
}


# mixture log normal distribution
z_generateLNM <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitLN(comb_xy)
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  if (length(f)==7 & all(as.data.frame(is.na(f))==FALSE)){
    samp=f[5]*rlnorm(num_l, f[1], f[3]) + (1-f[5])*rlnorm(num_l, f[2], f[4])
  }else{
    samp=rep(comb_xy, k)
  }
  
  return(samp)
}



# normal distribution
z_generateNL <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'normal')
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  
  if (length(f)==num_l){
    samp=f
  }else{
    samp=rnorm(num_l, f[["estimate"]][["mean"]], f[["estimate"]][["sd"]]) 
  }
  
  return(samp)
}


# lognormal distribution
z_generateLN <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'lognormal')
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  
  if (length(f)==num_l){
    samp=f
  }else{
    samp=rlnorm(num_l, f[["estimate"]][["meanlog"]], f[["estimate"]][["sdlog"]]) 
  }
  
  return(samp)
}


z_generateNB <- function(comb_xy, k, num_l){
  
  if (all(comb_xy==comb_xy[1])){
    f = rep(comb_xy, k)
  }else{
    f = tryCatch({
      fitdistr(comb_xy,'negative binomial')
      
    },error = function(e){
      rep(comb_xy, k)
    })
  }
  
  
  if (length(f)==num_l){
    samp=f
  }else{
    samp = rnbinom(num_l,size=f[["estimate"]][["size"]],mu=f[["estimate"]][["mu"]])
  }
  
  return(samp)
}

