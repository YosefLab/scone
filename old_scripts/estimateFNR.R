Lik_Fun = function(Z,X,Beta){
  return((1-Z) + (2*Z - 1)/(1+exp(-X %*% Beta)))
}

Post_Z = function(Y,W,Alpha,X,Beta){
  return(Lik_Fun(Y,W,Alpha)*Lik_Fun(1,X,Beta)/(Lik_Fun(Y,W,Alpha)*Lik_Fun(1,X,Beta) + (Y == 0)*Lik_Fun(0,X,Beta)))
}

estimateFNR = function(x, bulk_model = FALSE, is.control = TRUE, K = 0, X = NULL, Alpha = NULL,
                       is.expressed = NULL, niter = 100, ntresh = 0){
  
  require(matrixStats)
  Y = t(matrix(as.numeric(x > ntresh ),nrow = dim(x)[1]))
  
  if(bulk_model){
    x2 = x
    x2[x2 <= ntresh] = NA
    Alpha = t(matrix(log(rowMedians(x2, na.rm=TRUE))))
  }
  
  # Prepare Wanted Variation and Gene-Intrinsic Intercept
  X = cbind(X,rep(1,dim(x)[2]))
  is.intercept = dim(X)[2] == 1
  Beta = matrix(0,dim(X)[2],dim(x)[1])
  
  # Prepare Unwanted Variation (K factors) and Sample-Intrinsic Intercept
  if( K == 0 || !is.null(Alpha) ){
    Alpha = rbind(Alpha,rep(1,dim(x)[1]))
    K = dim(Alpha)[1] - 1
  }else{
    M = 2*Y[is.control,] - 1
    M = t(t(M - rowMeans(M)) - colMeans(M)) + mean(M)
    Alpha = rbind(t(svd(M,nv = K)$v[,1:K]),rep(1,dim(x)[1]))
  }
  W = matrix(0,dim(x)[2],K+1)
  
  # Initialize Z
  Z = Post_Z(Y,W,Alpha,X,Beta)
  
  # If control genes are supplied, W fit once, according to control genes
  if(!is.null(is.expressed)){
    cY = Y[,is.expressed]
    cAlpha = Alpha[,is.expressed]
    for (i in 1:dim(Z)[1]){
      W[i,] = glm(cbind(cY[i,],1-cY[i,]) ~ 0 + t(cAlpha),quasibinomial)$coefficients
    }   
  }
  
  ## DR: it seems strange to me that this is a for loop and not a while !converged
  for(n in 1:niter){    
    # Update Beta
    if(!is.intercept){
      
      if(is.null(is.expressed)){
        for(i in 1:dim(Z)[2]){
          Beta[,i] = glm(cbind(Z[,i],1-Z[,i]) ~ 0 + X,quasibinomial)$coefficients
        }
        
      }else{ # Only over non-control genes
        for(i in (1:dim(Z)[2])[!is.expressed]){
          Beta[,i] = glm(cbind(Z[,i],1-Z[,i]) ~ 0 + X,quasibinomial)$coefficients
        }
      }
      
      
    }else{
      # Short-cut!
      Beta = t(matrix(log(colMeans(Z)/(1-colMeans(Z)))))
    }
    
    # Update Z
    Z = Post_Z(Y,W,Alpha,X,Beta)
  
    if(is.null(is.expressed)){
      # Update W
      for (i in 1:dim(Z)[1]){
        W[i,] = glm(cbind(Y[i,],(1-Y[i,])*Z[i,]) ~ 0 + t(Alpha),quasibinomial)$coefficients
      }
    
      # Update Z
      Z = Post_Z(Y,W,Alpha,X,Beta)
    }
  }
  
  # Stick to our initial assumptions
  if(!is.null(is.expressed)){
    Beta[,is.expressed] = NA
    Z[,is.expressed] = 1
  }
  EL = sum(na.omit(as.vector(Z*log(Lik_Fun(Y,W,Alpha))))) + sum(na.omit(as.vector(Z*log(Lik_Fun(1,X,Beta))))) + sum(na.omit(as.vector((1-Z)*log(Lik_Fun(0,X,Beta)))))
  
  return(list(Z = t(Z), EL = EL, P = t(Lik_Fun(1,W,Alpha)), W = W, X = X, Alpha = Alpha, Beta = Beta))
}