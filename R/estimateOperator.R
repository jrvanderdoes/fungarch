#' Title
#'
#' @param X
#' @param w
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
vecARCH <- function(X, w=0.001,
                    lambda=c(10^-7,10^-6,10^-5,10^-4,
                             seq(0.001,0.01,0.002),
                             seq(0.02,0.09,0.01),
                             seq(0.1,5,0.1),seq(6,20,1),
                             seq(22,50,2),
                             seq(60,100,10),seq(125,250,25),
                             seq(300,500,50))){

  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution
  X <- .tensor(X) # Tensor product square

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_tensor_lag_cov(X)
  C_1 <- .estim_tensor_lag_cov(X, lag = 1)

  eigval_C_0 <- eigen(C_0)

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

  # Standardize Vectors in integral
  #   apply( eigval_C_0$vectors^2,MARGIN = 2,.approx_int)
  ## TODO:: Would I want the tensor square to be 1 now?
  eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)

  # Compute Projection
  projector <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K){
    projector <- projector +
      ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
  }

  # Prep
  I_H <- diag(n)

  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X^2)

  # Estimate alpha and delta
  error_plot <- NULL
  if(length(lambda)>1){
    # optim_results <-
    #   optim(0.1,function(lambda, X, N, C_0, C_1, projector, m2, I_H){
    #     # Formula (3.26)
    #     invMat <- solve( sqrt(N) * C_0 %*% C_0 + lambda * I_H )
    #     #   Formula (4.8)
    #     hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector
    #     #hat_alpha <- apply(hat_alpha,MARGIN = 2, function(x){pmax(0,x)})
    #
    #     # Estimate delta, (4.9)
    #     hat_delta <- m2 - hat_alpha %*% m2
    #     #hat_delta <- pmax(0,hat_delta)
    #
    #     # Compute Error
    #     .computeError(X, hat_alpha, hat_delta)
    #   },
    # method="L-BFGS-B",
    # lower=10^-7,upper=10,
    # control = list('ndeps'=10^-6),
    # X=X,N=N,C_0=C_0,C_1=C_1,projector=projector, m2=m2, I_H=I_H)

    lambda_result <- .findLambdaARCH(lambda,X=X,C_0=C_0,C_1=C_1,I_H=I_H,N=N,
                                     projector=projector,m2=m2)
    error_plot <- lambda_result$error_plot
    lambda <- lambda_result$lambda
  }


  # Formula (3.26)
  invMat <- solve( sqrt(N) * C_0 %*% C_0 + lambda * I_H )
  #   Formula (4.8)
  hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector
  #hat_alpha <- apply(hat_alpha,MARGIN = 2, function(x){pmax(0,x)})

  # Estimate delta, (4.9)
  hat_delta <- m2 - hat_alpha %*% m2
  #hat_delta <- pmax(0,hat_delta)

  list('hat_alpha'=hat_alpha,
       'hat_delta'=hat_delta,
       'fitting_params'=list('lambda'=lambda,
                             'C_0'=C_0, 'C_1'=C_1,
                             'invMat'=invMat,
                             'eigval_C_0'=eigval_C_0,
                             'K'=K,
                             'projector'=projector,
                             'error_plot'=error_plot))
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
  }

  result
}


#' Title
#'
#' @param X
#' @param lag
#'
#' @return
#' @export
#'
#' @examples
.estim_tensor_lag_cov <- function(X, lag=0){
  N <- ncol(X)

  # TODO:: Confirm if use definition from (3.3) or after (3.6)?
  m1 <- (1/(N-lag)) * rowSums(X[,1:(N-lag)])

  covMat <- matrix(0, nrow=nrow(X),ncol=nrow(X))

  for(k in 1:(N-lag)){
    covMat <- covMat + ( (X[,k]-m1) %*% t(X[,k+lag]-m1) )
  }

  covMat /((N-lag)-1)
}
