#' Fit a pointwise ARCH model
#'
#' Estimate pointwise autoregressive coefficients and innovation variances
#' after squaring and de-meaning the functional observations.
#'
#' @param X Finite numeric matrix of functional observations, with observations
#'   in columns.
#' @param w Approximation threshold used for covariance-basis selection.
#' @param inv_type Inversion method: `"M"` for the truncated covariance
#'   inverse or `"T"` for Tikhonov regularization.
#' @param training Positive fraction of observations used for fitting when
#'   `inv_type = "T"`. Values below 1 reserve the remaining observations for
#'   tuning the regularization parameter; values at least 1 use all data.
#'
#' @return A list containing the estimated `alpha` matrix and pointwise
#'   innovation variances `delta`.
#'
#' @noRd
.fit_pw <- function(X, w=0.01, inv_type=c('M','T'), training=0.8){
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix", call. = FALSE)
  }
  if (ncol(X) < 3L) {
    stop("X must contain at least three observations", call. = FALSE)
  }
  if (length(w) != 1L || !is.numeric(w) || !is.finite(w) || w < 0 || w > 1) {
    stop("w must be a single finite value between 0 and 1", call. = FALSE)
  }
  if (length(training) != 1L || !is.numeric(training) ||
      !is.finite(training) || training <= 0) {
    stop("training must be a single positive finite value", call. = FALSE)
  }
  ## Setup Vars
  # inv_type <- inv_type[1]
  inv_type <- match.arg(inv_type, c("M", "T"))

  # Setup Data
  X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square

  N <- ncol(X)
  n <- nrow(X) # Resolution


  # Select K Projections
  if(inv_type=='M'){
    # Lag-h Covariance operator estimates
    #   Formula (3.6)
    C_0 <- .estim_lag_cov(X)
    C_1 <- .estim_lag_cov(X, lag = 1)

    # Eigenvector to Eigenfunctions
    eigval_C_0 <- eigen(C_0)
    eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)
    eigval_C_0$values <- eigval_C_0$values / n

    # Important Number
    K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

    # Inverse
    C_eps_dag <- matrix(0, nrow=nrow(X), ncol=nrow(X))
    for(i in 1:K){
      C_eps_dag <- C_eps_dag +
        ( 1 / eigval_C_0$values[i] ) *
        (eigval_C_0$vectors[,i] %o% eigval_C_0$vectors[,i] )
    }
  }else if(inv_type=='T'){

    ## USE PART OF DATA FOR FIT
    if(training < 1){
      # OLD:
      # X_part <- X[,1:floor(training * ncol(X))]
      # X_pred <- X[,floor(training * ncol(X)+1):ncol(X)]
      n_train <- floor(training * ncol(X))
      if (n_train < 3L || n_train >= ncol(X)) {
        stop("training must leave at least three fitting observations and prediction data",
             call. = FALSE)
      }
      X_part <- X[, seq_len(n_train), drop = FALSE]
      X_pred <- X[, seq.int(n_train + 1L, ncol(X)), drop = FALSE]
    }else{
      X_part <- X
      X_pred <- X
    }

    # Lag-h Covariance operator estimates
    #   Formula (3.6)
    C_0 <- .estim_lag_cov(X_part)
    C_1 <- .estim_lag_cov(X_part, lag = 1)

    # Eigenvector to Eigenfunctions
    eigval_C_0 <- eigen(C_0)
    eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)
    eigval_C_0$values <- eigval_C_0$values / n

    # Important Number
    K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

    # Compute Projection
    projector <- matrix(0, nrow=n, ncol=n)
    for(i in 1:K){
      projector <- projector +
        ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
    }

    optim_results <-
      optim(10000,function(lambda, X1, N, n, C_0, C_1, projector, m2){
        ests <- .estimateARCH(N=N, n=n, C_0=C_0, C_1=C_1,
                              lambda=lambda, projector=projector, m2=m2)

        .computeARCHError(X = X1, alpha = ests$hat_alpha, delta = ests$hat_delta)
      },
      method="Brent",
      lower=10^-7, upper=500000,
      X1=X_pred, N=ncol(X_part), n=n,
      C_0=C_0, C_1=C_1, projector=projector, m2=rowMeans(X_part))

    lambda <- optim_results$par

    # Lag-h Covariance operator estimates
    #   Formula (3.6)
    C_0 <- .estim_lag_cov(X)
    C_1 <- .estim_lag_cov(X, lag = 1)

    eigval_C_0 <- eigen(C_0)
    eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)

    # Select K Projections
    #   Remark 3.2(c)
    K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

    # Compute Projection
    projector <- matrix(0, nrow=n, ncol=n)
    for(i in 1:K){
      projector <- projector +
        ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
    }

    # invMat <- solve( C_0 + diag(lambda, n, n) )
    # hat_alpha <- C_1 %*% invMat %*% projector
    C_eps_dag <- solve( C_0 + diag(lambda, n, n) )
  }

  #   Estimate Alpha/Delta
  hat_alpha <- C_1 %*% C_eps_dag / n
  hat_delta <- pmax(0, rowMeans(X) - (hat_alpha %*% rowMeans(X)) / n)


  list('alpha'=hat_alpha,
       'delta'=hat_delta)
}

#' Predict a pointwise conditional covariance
#'
#' @param X Numeric matrix of recent functional observations.
#' @param alpha Estimated pointwise autoregressive coefficient matrix.
#' @param delta Estimated pointwise innovation variances.
#' @param C_eps Innovation covariance matrix or correlation structure.
#'
#' @return A predicted covariance matrix.
#'
#' @noRd
.pred_pw <- function(X, alpha, delta, C_eps){
  # OLD: The estimate calculation is retained here for comparison.
  # X2 <- X^2
  # estimate <- rep(NA,nrow(X2))
  # for(j in 1:nrow(X2)){
  #   vector_op <- alpha[j,] * X2[,ncol(X2)] # X^2
  #   estimate[j] <- delta[j] + fChange:::dot_integrate(as.matrix(vector_op))
  # }
  # estimate <- ifelse(estimate>0, estimate, 0)

  estimate <- .pointwise_variance(X, alpha, delta)

  sqrt(estimate) %*% t(sqrt(estimate)) * C_eps
  # as.vector( t(sqrt(estimate) %*% C_eps) %*% sqrt(estimate) ) / length(estimate)
}

#' Predict pointwise conditional standard deviations
#'
#' @param X Numeric matrix of recent functional observations.
#' @param alpha Estimated pointwise autoregressive coefficient matrix.
#' @param delta Estimated pointwise innovation variances.
#' @param C_eps Unused compatibility argument retained for existing workflows.
#'
#' @return A vector of predicted pointwise standard deviations.
#'
#' @noRd
.pred_pw_simple <- function(X, alpha, delta, C_eps){
  # OLD: The estimate calculation is retained here for comparison.
  # X2 <- X^2
  # estimate <- rep(NA,nrow(X2))
  # for(j in 1:nrow(X2)){
  #   vector_op <- alpha[j,] * X2[,ncol(X2)] # X^2
  #   estimate[j] <- delta[j] + fChange:::dot_integrate(as.matrix(vector_op))
  # }
  # estimate <- ifelse(estimate>0, estimate, 0)

  estimate <- .pointwise_variance(X, alpha, delta)

  sqrt(estimate)
  # as.vector( t(sqrt(estimate) %*% C_eps) %*% sqrt(estimate) ) / length(estimate)
}

#' Compute pointwise conditional variances for the prediction helpers.
#'
#' @param X Numeric matrix of recent functional observations.
#' @param alpha Estimated pointwise autoregressive coefficient matrix.
#' @param delta Estimated pointwise innovation variances.
#'
#' @return A nonnegative vector of pointwise conditional variance estimates.
#'
#' @noRd
.pointwise_variance <- function(X, alpha, delta){
  X2 <- X^2
  estimate <- numeric(nrow(X2))

  for (j in seq_len(nrow(X2))) {
    vector_op <- alpha[j, ] * X2[, ncol(X2)]
    estimate[j] <- delta[j] +
      fChange:::dot_integrate(as.matrix(vector_op))
  }

  pmax(estimate, 0)
}

#' Estimate ARCH Parameters
#'
#' This (internal) function computes the alpha and delta parameters for a
#'  ARCH(1) process
#'
#' @param N Numeric for the total number of observations.
#' @param n Numeric for the resolution.
#' @param C_0 Matrix of numerics for the lag 0 covariance operator.
#' @param C_1 Matrix of numerics for the lag 1 covariance operator.
#' @param lambda Numeric for the smoothing parameter of Tikhonov
#' @param projector Matrix of numerics for projection
#' @param m2 Vector of numerics for average squared values
#'
#' @return List of parameter estimates for ARCH(1):
#'  \enumerate{
#'    \item hat_alpha: Estimate for alpha
#'    \item hat_delta: Estimate for delta
#'  }
#'
#' @noRd
.estimateARCH <- function(N, n, C_0, C_1, lambda, projector, m2){
  # Formula (3.26)
  invMat <- solve( sqrt(N) * C_0 %*% C_0 + diag(lambda, n, n) )
  #   Formula (4.8)
  hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector

  # Estimate delta, (4.9)
  hat_delta <- pmax(0,m2 - (hat_alpha %*% m2) / n)

  list('hat_alpha'=hat_alpha,
       'hat_delta'=hat_delta)
}

#' Predict ARCH Estimates
#'
#' This (internal) function estimates given X^2, alpha, and delta
#'  of an ARCH(1) process.
#'
#' @param X2 Matrix of numerics indicating the squared process
#' @param alpha Matrix of numerics estimating the alpha values
#' @param delta Matrix of numerics estimating the delta values
#' @param est2 If `TRUE`, return conditional variance estimates; otherwise
#'   return their pointwise square roots after truncating negative values.
#'
#' @return Estimate ARCH(1) based on data with given alpha and delta
#'
#' @noRd
.predictARCH <- function(X2, alpha, delta, est2=TRUE){
  estimate <- matrix(0, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- delta

  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j,] * X2[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + fChange:::dot_integrate(as.matrix(vector_op))
      #pracma::trapz(x=seq(0,1,length.out=length(vector_op)), y=vector_op)
    }
  }

  if(est2) return(estimate)

  estimate <- ifelse(estimate>0,estimate,0)

  sqrt(estimate)
}


#' Compute ARCH Error
#'
#' This (internal) function estimates values and then compares the estimates
#'  to the known values of an ARCH(1,1) process. The functional L2 error is
#'  estimated.
#'
#' @param X Matrix of numerics indicating the squared process
#' @inheritParams .predictARCH
#'
#' @return  Estimate of error based on data with given alpha and delta
#'
#' @noRd
.computeARCHError <- function(X, alpha, delta){
  error <- rep(0, ncol(X))
  estimate <- .predictARCH(X,alpha,delta,est2=TRUE)

  for(i in 2:ncol(X)){
    error[i] <- fChange:::dot_integrate(abs(X[,i] - estimate[,i]))
  }

  sum(error)
}
