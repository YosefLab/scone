#' Parameter estimation of zero-inflated negative binomial model
#' 
#' This function implements an EM algorithm to estimate the parameters of a zero-inflated negative binomial model
#' 
#' This function implements an expectation-maximization algorithm for a zero-inflated
#' model, with a constant gene-specific mean parameter, gene-specific dispersion and
#' a mixing probability that depends solely on the mean.
#'  
#' @param Y matrix. The data matrix (genes in rows, cells in columns)
#' @param maxiter numeric. The maximum number of iterations.
#' @param verbose logical. Whether or not to print the value of the likelihood at each value.
#' 
#' @importFrom MASS glm.nb
#' @export
#' 
#' @return a list with the following elements:
#' \itemize{
#' \item{mu}{the mean of the negative binomial component}
#' \item{theta}{the dispersion parameter of the negative binomial component}
#' \item{pi}{the mixing probability}
#' \item{coefs}{the coefficients of the logistic regression}
#' \item{p_z}{the probability P(Z=0|Y), useful as weights in weighted principal components analysis}
#' \item{expected_value}{the expected value E[Y]}
#' \item{variance}{the variance Var(Y)}
#' \item{loglik}{the log-likelihood}
#' \item{convergence}{0 if the algorithm converged and 1 if maxiter was reached}
#' }
estimate_zinb <- function(Y, maxiter=10, verbose=FALSE) {
  
  n <- ncol(Y)
  J <- nrow(Y)
  
  # 0. Initial estimate of mu, pi and z
  mu0 <- rowMeans(Y)
  
  pi0 <- sapply(seq_len(n), function(i) {
    y <- as.numeric(Y[,i]<=0)
    fit <- glm(y~log(mu0), family = binomial)
    return(fitted.values(fit))
  })
    
  zhat <- pi0/(pi0 + (1 - pi0) * dnbinom(0, size = 1, mu = matrix(mu0, nrow=J, ncol=n)))
  zhat[Y>0] <- 0
  thetahat <- rep(1, J)
  
  coefs_mu <- log(mu0)
  coefs_pi <- sapply(seq_len(n), function(i) {
    fit <- suppressWarnings(glm(zhat[,i]~log(mu0), family = binomial(link = logit)))
    return(coefficients(fit))
  })
 
  X <- matrix(rep(1, n), ncol=1)
  W <- model.matrix(~log(mu0))
  linkobj <- binomial()
  
  ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, 0, 0, linkobj)
  ll_old <- 2 * ll_new
  
  ## EM iteration
  iter <- 0
  while (abs((ll_old - ll_new)/ll_old) > 1e-4 & iter<maxiter) {
    ll_old <- ll_new
    fit_mu <- bplapply(seq_len(J), function(i) {
      fit <- suppressWarnings(glm.nb(Y[i,] ~ 1, weights = (1 - zhat[i,]), init.theta = thetahat[i], start=coefs_mu[i]))
      return(list(coefs=coefficients(fit), theta=fit$theta))
    })
    coefs_mu <- sapply(fit_mu, function(x) x$coefs)
    muhat <- exp(coefs_mu)
    thetahat <- sapply(fit_mu, function(x) x$theta)
    
    fit_pi <- bplapply(seq_len(n), function(i) {
      fit <- suppressWarnings(glm(zhat[,i]~log(muhat), family = binomial(link = logit), start=coefs_pi[,i]))
      return(list(coefs=coefficients(fit), fitted=fitted.values(fit)))
    })
    coefs_pi <- sapply(fit_pi, function(x) x$coefs)
    pihat <- sapply(fit_pi, function(x) x$fitted)
    
    zhat <- pihat/(pihat + (1 - pihat) * dnbinom(0, size = matrix(thetahat, nrow=J, ncol=n), mu = matrix(muhat, nrow=J, ncol=n)))
    zhat[Y>0] <- 0
    
    W <- model.matrix(~log(muhat))
    ll_new <- loglik_small(c(coefs_mu, coefs_pi[1,], coefs_pi[2,], log(thetahat)), Y, Y>0, X, W, J, n*2, 0, 0, linkobj)
    
    if(verbose) {
      print(ll_new)
    }
    
    iter <- iter + 1
  }
  
  convergence <- 0
  if(iter == maxiter) {
    convergence <- 1
  }
  
  p_z <- 1 - (Y == 0) * pihat / (pihat + (1 - pihat) * (1 + muhat / thetahat)^(-thetahat))
  eval <- (1 - pihat) * muhat
  variance <- (1 - pihat) * muhat * (1 + muhat * (1/thetahat + pihat))
  
  return(list(mu=muhat, theta=thetahat, pi=pihat, coefs=coefs_pi,
              p_z=p_z, expected_value=eval, variance=variance, 
              loglik=ll_new, convergence=convergence))
}

#' Log-likelihood function of the zero-inflated negative binomial model
#' 
#' This function computes the log-likelihood of a standard regression model
#' 
#' This is a (hopefully) memory-efficient implementation of the log-likelihood of a 
#' zero-inflated negative binomial regression model.
#' In this attempt, the design matrices don't have n*J rows, but n and J, respectively.
#' The computation is a bit slower, but the memory usage should be much smaller for
#' large J and n.
#' 
#' @param parms a vector of parameters: should contain the values of beta, followed by those of alpha, followed by the log(1/phi)
#' @param Y the data matrix (genes in rows, cells in columns)
#' @param Y1 a logical indicator of Y>0
#' @param X the design matrix for the regression on mu (n x k_X)
#' @param W the design matrix for the regression on pi (J x k_W)
#' @param kx the number of beta parameters
#' @param kw the number of alpha parameters
#' @param offsetx the offset for the regression on X
#' @param offsetw the offset for the regression on W
#' @param linkobj the link function object for the regression on pi (typically the result of binomial())
loglik_small <- function(parms, Y, Y1, X, W, kx, kw, offsetx, offsetw, linkobj) {
  
  J <- nrow(Y)
  n <- ncol(Y)
  
  beta <- matrix(parms[1:kx], ncol=J, nrow=ncol(X))
  mu <- t(exp(X %*% beta + offsetx))
  
  alpha <- matrix(parms[(kx + 1):(kw+kx)], nrow=ncol(W), ncol=n, byrow=TRUE)
  eta <- W %*% alpha + offsetw
  pi <- logistic(eta)
  
  theta <- matrix(exp(parms[(kw + kx + 1):(kw + kx + J)]), nrow=J, ncol=n)
  
  loglik0 <- log(pi + exp(log(1 - pi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - pi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  
  return(sum(loglik0[which(Y==0)]) + sum(loglik1[which(Y>0)]))
}

#' Imputation of zero counts based on zero-inflated negative binomial model
#' 
#' This function is used to impute the data, by weighting the observed zeros by their
#' probability of coming from the zero-inflation part of the distribution P(Z=1|Y).
#'
#' @details The imputation is carried out with the following formula:
#' y_{ij}* =  y_{ij} * Pr(Z_{ij} = 0 | Y_{ij}=y_{ij}) + mu_{ij} * Pr(Z_{ij} = 1 | Y_{ij} = y_{ij}).
#' Note that for y_{ij} > 0, Pr(Z_{ij} = 0 | Y_{ij}=y_{ij}) = 1 and hence the data are not imputed.
#'  
#' @export
#' 
#' @param expression the data matrix (genes in rows, cells in columns)
#' 
#' @return the imputed expression matrix.
impute_zinb <- function(expression) {
  pars <- estimate_zinb(expression)
  w <- pars$p_z
  imputed <- expression * w + pars$mu * (1 - w)
  return(imputed)
}

## change the above to be a closure inside the main function so that we can
## specify the number of clusters (or use bplapply?)

logit <- binomial()$linkfun
logistic <- binomial()$linkinv