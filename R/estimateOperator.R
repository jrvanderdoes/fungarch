#' Operator ARCH Estimation
#'
#' @param X
#' @param w
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
operatorARCH <- function(X, w=0.001, lambda=NA){
  X <- X - rowMeans(X) # De-mean
  X <- .tensor(X) # Pointwise square
  #X <- X - rowMeans(X) # De-mean
  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # Lag-h Covariance operator estimates
  C_0 <- .estim_lag_cov(X)
  C_1 <- .estim_lag_cov(X, lag = 1)

  eigval_C_0 <- eigen(C_0)

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

  # Standardize Vectors in integral
  #   apply( eigval_C_0$vectors^2,MARGIN = 2,.approx_int)
  eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)

  # Compute Projection
  projector <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K){
    projector <- projector +
      ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
  }

  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(.tensor(X))

  # Estimate alpha and delta
  error_plot <- NULL
  if(is.na(lambda)){
    optim_results <-
      optim(10000,function(lambda, X, N, n, C_0, C_1, projector, m2){

        ests <- .estimateARCH(N=N, n=n, C_0=C_0, C_1=C_1,
                              lambda=lambda, projector=projector, m2=m2)

        .computeARCHError(X = X, alpha = ests$hat_alpha, delta = ests$hat_delta)
      },
      method="Brent",
      lower=10^-7,upper=500000,
      X=X,N=N,n=n,C_0=C_0,C_1=C_1,projector=projector,m2=m2)

    lambda <- optim_results$par
  }

  ests <- .estimateARCH(N=N, n=n, C_0=C_0, C_1=C_1,
                        lambda=lambda, projector=projector, m2=m2)

  list('hat_alpha'=ests$hat_alpha,
       'hat_delta'=ests$hat_delta,
       'estimate'=.predictARCH(X2=X, alpha=ests$hat_alpha,
                               delta=ests$hat_delta, est_sqr=FALSE),
       'fitting_params'=list('lambda'=lambda,
                             'C_0'=C_0, 'C_1'=C_1,
                             'eigval_C_0'=eigval_C_0,
                             'K'=K,
                             'projector'=projector,
                             'error_plot'=error_plot) )
}


#' Title
#'
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
.tensor <- function(X,Y=X){
  ## TODO:: Discuss this
  if(is.vector(X) && is.vector(Y)){
    result <- X %*% t(Y)
  } else{
    result <- list()
    for(i in 1:ncol(X)){
      result[[i]] <- X[,i] %*% t(Y[,i])
    }
    # Create Tensor Matrix
    tensorVal <- do.call("cbind",result)
    for(i in 2:ncol(X)){
      # Build row
      row <- result[[i]]
      for(j in 2:ncol(X)){
        func <- ifelse(j<i,t,function(X){X})
        row <- cbind(row,func(result[[abs(i-j)+1]]))
      }
      tensorVal <- rbind(tensorVal, row)
    }
  }

  result
}
