.fit_pw <- function(X, w=0.01){
  X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square

  ## USE PART OF DATA FOR FIT
  N <- ncol(X)
  n <- nrow(X) # Resolution

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_lag_cov(X)
  C_1 <- .estim_lag_cov(X, lag = 1)

  eigval_C_0 <- eigen(C_0)
  eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)
  eigval_C_0$values <- eigval_C_0$values / n

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)
  C_eps_dag <- matrix(0, nrow=nrow(X), ncol=nrow(X))
  for(i in 1:K){
    C_eps_dag <- C_eps_dag +
      ( 1 / eigval_C_0$values[i] ) *
      (eigval_C_0$vectors[,i] %o% eigval_C_0$vectors[,i] )
  }

  #   Estimate Alpha/Delta
  m2 <- 1 / N * rowSums(X)

  hat_alpha <- C_1 %*% C_eps_dag / n
  hat_delta <- pmax(0,m2 - (hat_alpha %*% m2) / n)


  list('alpha'=hat_alpha,
       'delta'=hat_delta)
       # 'estimate'=.predictARCH(X2=X, alpha=ests$hat_alpha,
       #                         delta=ests$hat_delta, est2 = FALSE),
       # 'fitting_params'=list('lambda'=lambda,
       #                       'C_0'=C_0, 'C_1'=C_1,
       #                       'eigval_C_0'=eigval_C_0,
       #                       'K'=K,
       #                       'projector'=projector))
}

.pred_pw <- function(X, alpha, delta, C_eps){
  X2 <- X^2

  estimate <- rep(NA,nrow(X2))

  for(j in 1:nrow(X2)){
    vector_op <- alpha[j,] * X2[,ncol(X2)] # X^2
    estimate[j] <- delta[j] + dot_integrate(as.matrix(vector_op))
  }
  estimate <- ifelse(estimate>0, estimate, 0)

  sqrt(estimate) %*% t(sqrt(estimate)) * C_eps
  # as.vector( t(sqrt(estimate) %*% C_eps) %*% sqrt(estimate) ) / length(estimate)
}
