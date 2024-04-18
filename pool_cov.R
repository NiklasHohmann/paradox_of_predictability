# function to compute pooled covariance matrix from a list, and flexibly return
# a new list with some or all elements replaced by the pooled S
# all elements based on fewer than min.N observations (from nn vector) will
# be replaced by the pooled matrix.  
# S is a list of covariance matrices, nn is corresponding vector of sample sizes
# if return.list = TRUE, the modified list of covariance matrices is returned
# if return.list = FALSE, the pooled covariance matrix alone is returned
pool_cov <- function(S, nn, return.list = TRUE, min.N = 50){
  nS <- length(S)
  nN <- length(nn)
  if(nS != nN) stop("Lengths of S and nn do not match.") else nsamp <- nS
  
  Ntot <- sum(nn)
  nvar <- ncol(S[[1]])
  Sp <- matrix(0, nrow = nvar, ncol = nvar)
  for(i in 1:nsamp) 
    Sp <- Sp + (nn[i]-1) * S[[i]]
  Sp <- Sp/(Ntot - nsamp)
  
  if(!return.list) return(Sp) else{
    Sret <- S
    lowN <- nn < min.N
    for(i in 1:nsamp)
      if(lowN[i]) Sret[[i]] <- Sp else Sret[[i]] <- S[[i]]
    
    return(Sret) 
  }
}


