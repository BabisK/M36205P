recinvert <- function(X) {
  if (!is.matrix(X)) {
    return(X)
  }
  
  if (ncol(X) == 2) {
    a <- nrow(X)
    b <- sum(X[,2])
    d <- sum(X[,2]^2)
    XTX <- matrix(c(a, b, b, d), nrow = 2, ncol = 2)
    XTXi <- solve(XTX)
    return(XTXi)
  }
  else{
    Ai <- recinvert(X[,-ncol(X)])
    B <- crossprod(X[,ncol(X)], X[,-ncol(X)])
    B <- t(B)
    d <- sum(X[,ncol(X)]*X[,ncol(X)])
    k <- d-crossprod(B,Ai%*%B)
    
    final <- rbind(cbind(Ai + 1/k[1,1]*(Ai%*%B%*%crossprod(B,Ai)), -1/k[1,1]*(Ai%*%B)),cbind(-1/k[1,1]*(crossprod(B,Ai)), 1/k[1,1]))
    
    return(final)
  }
}

rand <- replicate(100, rnorm(10000))
rand <- cbind(replicate(10000, 1), rand)

i1 <- system.time(solve(crossprod(rand)))

i2 <- system.time(recinvert(rand))
