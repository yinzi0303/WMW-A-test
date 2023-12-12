wmwa <- function(data, m, n, N, nperm){
  
  x <- data[1:m]
  y <- data[(m+1):(m+n)]
  z <- data[(m+n+1):N]
  comb_xyz <- c(x, y, z)
  
  if (all(comb_xyz==comb_xyz[1])){
    pVal = 1
  }else{
    allmin <- min(c(x, y))
    allmax <- max(c(x, y))
    idx <- intersect(which(z >= allmin), which(z <= allmax))
    z = z[idx]
    # z = z[sample(length(idx))[1:2*(m+n)]]
    
    l <- length(z)
    N <- m + n + l
    
    wmwa_obs_x <- sum(rank(c(x, y, z))[1:m]) 
    wmwa_obs_y <- sum(rank(c(y, x, z))[1:n]) 
    wmwa_obs <- wmwa_obs_x/m - wmwa_obs_y/n
    
    comb_xyz <- c(x, y, z)
    label <- NULL
    label[1:m] <- 1
    label[(m+1):(m+n)] = 2
    label[(m+n+1):N] = 3
    
    wmwa_null <- NULL
    for (j in 1:nperm){
      plabel <- sample(label, N)
      px <- comb_xyz[plabel==1]
      py <- comb_xyz[plabel==2]
      pz <- comb_xyz[plabel==3]
      Wx <- sum(rank(c(px, py, pz))[1:m]) 
      Wy <- sum(rank(c(py, px, pz))[1:n]) 
      wmwa_null[j] <- Wx/m - Wy/n
    }
    
    
    p_left <- sum(wmwa_null <= wmwa_obs) / nperm
    p_right <- sum(wmwa_null >= wmwa_obs) / nperm
    p_tail = min(p_left, p_right)
    pVal = min(2 * p_tail, 1)
  }
  
  return(pVal)
  
}
