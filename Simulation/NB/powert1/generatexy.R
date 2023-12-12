generatexy <- function (Times = 1000, m = 5, n = 5, k = 5, delta = 1, alphax = 0, 
                        alphay=1,beta = 1, distribution = "Normal"){
  X <- matrix(data = NA, nrow = Times, ncol = m)
  Y <- matrix(data = NA, nrow = Times, ncol = n)
  Z <- matrix(data = NA, nrow = Times, ncol = k * (m + n))
  l1 = k * m
  l2 = k * n
  if (distribution=="Normal"){
    for (i in 1:Times){
      X[i,] <- rnorm(m, alphax, beta)
      Y[i,] <- rnorm(n, alphay, beta)
      Z[i,] <- c(rnorm(l1, alphax, beta), rnorm(l2, alphay, beta))
    }
  }
  else if (distribution=="Gamma"){
    for (i in 1:Times){
      X[i,] <- rgamma(m, alphax, beta)
      Y[i,] <- rgamma(n, alphay, beta)
      Z[i,] <- c(rgamma(l1, alphax, beta), rgamma(l2, alphay, beta))
    }
  }
  else if (distribution=="NB1"){
    for (i in 1:Times){
      X[i,] <- rnbinom(m, prob=alphax, size=beta) 
      Y[i,] <- rnbinom(n, prob=alphay, size=beta)
      Z[i,] <- c(rnbinom(l1, prob=alphax, size=beta), rnbinom(l2, prob=alphay, size=beta))
    }
  }
  else if (distribution=="NB2"){
    for (i in 1:Times){
      X[i,] <- rnbinom(n=m, size=alphax, prob=beta) 
      Y[i,] <- rnbinom(n=n, size=alphay, prob=beta)
      Z[i,] <- c(rnbinom(n=l1, size=alphax, prob=beta), rnbinom(n=l2, size=alphay, prob=beta))
    }
  }
  
  data <- list()
  data <- list(X = X, Y = Y, Z = Z)
  return(data)
}

