#' Pointwise Yule Walker ARCH(1) Estimation
#'
#' Estimates the parameters for ARCH(1) model using pointwise Yule Walker
#'  estimation.
#'
#' @param X Data.frame of functional data--curves as observations and rows as
#'  discretization
#' @param w (Optional) Numeric. Weight for choice of K that indicates the number
#'  of projectors to use.
#' @param lambda (Optional) Numeric. Amount of tikhonov regularation. Default is
#'  NA which uses optimazation to find the value. Alternatively many values can
#'  be given and tested, or a single value just used as given value.
#' @param training
#'
#' @return List of fit and some fitting parameters:
#'  \enumerate{
#'    \item hat_alpha: Estimated alpha estimate
#'    \item hat_delta: Estimated delta estimate
#'    \item estimate: Estimates from the estimated parameters
#'    \item fitting_params: List of parameters used in fitting
#'      \enumerate{
#'        \item lambda: Parameter for YW estimation
#'        \item C_0: Lag-0 covariance operator estimate
#'        \item C_1: Lag-1 covariance operator estimate
#'        \item eigval_C_0: Eigenvalues of C_0
#'        \item K: Number of projectors (i.e. dimension) to use
#'        \item projector:
#'        \item error_plot
#'      }
#'  }
#' @export
#'
#' @references Functional ARCH and GARCH models: A Yule-Walker approach.
#'  Sebastian K\"uhnert (TODO:: Complete reference)
#'
#' @examples
#' set.seed(12345)
#' alpha <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
#' delta <- function(t){ rep(0.1,length(t)) }
#' dat <- generateFunctionalCondVol(51,200,'arch',
#'           alpha_par = alpha, delta_par = delta)
#' est <- pointwiseARCH(X = dat$garch_mat)
#'
#' # s <- seq (0,1,length.out=nrow(dat$garch_mat))
#' # alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
#' # for(s1 in 1:length(s)){
#' #   for(s2 in 1:length(s)){
#' #     alpha_true[s1,s2] <- alpha(s[s1],s[s2])
#' #   }
#' # }
#' ## Compare Alpha
#' # plot3D::persp3D(z=alpha_true,main=parse(text = paste0('~ z == a(s,t)')),
#' #                 theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
#' # plot3D::persp3D(z=est$hat_alpha,main=bquote(~z == hat(a)(s,t)),
#' #                 theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
#' ## Compare Delta
#' # d <- seq(0,1,length.out=nrow(dat$garch_mat))
#' # delta_true <- delta(d)
#' # plot(x=d,delta_true,
#' #      ylim=c(min(delta_true, est$hat_delta),
#' #             max(delta_true,est$hat_delta)),
#' #      type='l')
#' # lines(x=d,y=est$hat_delta, col='red')
pointwiseARCH <- function(X, w=0.001, lambda=NA, training=0.8){
  X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square
  n <- nrow(X) # Resolution

  ## USE PART OF DATA FOR FIT
  if(training < 1){
    X_part <- X[,1:(training * ncol(X))]
    X_pred <- X[,(training * ncol(X)+1):ncol(X)]
  }else{
    X_part <- X
    X_pred <- X
  }
  N <- ncol(X_part)

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_lag_cov(X_part)
  C_1 <- .estim_lag_cov(X_part, lag = 1)

  eigval_C_0 <- eigen(C_0)

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

  # Standardize Vectors in integral
  #   apply( eigval_C_0$vectors^2,MARGIN = 2,dot_integrate)
  eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)

  # Compute Projection
  projector <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K){
    projector <- projector +
      ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
  }

  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X_part)

  optim_results <-
    optim(10000,function(lambda, X1, N, n, C_0, C_1, projector, m2){
      ests <- .estimateARCH(N=N, n=n, C_0=C_0, C_1=C_1,
                            lambda=lambda, projector=projector, m2=m2)

      .computeARCHError(X = X1, alpha = ests$hat_alpha, delta = ests$hat_delta)
    },
    method="Brent",
    lower=10^-7,upper=500000,
    X1=X_pred,N=N,n=n,C_0=C_0,C_1=C_1,projector=projector,m2=m2)

  lambda <- optim_results$par

  ## USE ALL DATA FOR FINAL SETUP
  N <- ncol(X) # Number of observations

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_lag_cov(X)
  C_1 <- .estim_lag_cov(X, lag = 1)

  eigval_C_0 <- eigen(C_0)

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)

  # Standardize Vectors in integral
  #   apply( eigval_C_0$vectors^2,MARGIN = 2,dot_integrate)
  eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)

  # Compute Projection
  projector <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K){
    projector <- projector +
      ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
  }
  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X)

  ests <- .estimateARCH(N=N, n=n, C_0=C_0, C_1=C_1,
                        lambda=lambda, projector=projector, m2=m2)

  list('hat_alpha'=ests$hat_alpha,
       'hat_delta'=ests$hat_delta,
       'estimate'=.predictARCH(X2=X, alpha=ests$hat_alpha,
                               delta=ests$hat_delta, est2 = FALSE),
       'fitting_params'=list('lambda'=lambda,
                             'C_0'=C_0, 'C_1'=C_1,
                             'eigval_C_0'=eigval_C_0,
                             'K'=K,
                             'projector'=projector))
}

#' Estimate ARCH Parameters
#'
#' This (internal) function computes the alpha and delta parameters for a
#'  ARCH(1) process
#'
#' @param N Numeric for the total number of observations.
#' @param n Numeric for the resolution.
#' @param c_0 Matrix of numerics for the lag 0 covariance operator
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
      estimate[j,i] <- delta[j] + dot_integrate(as.matrix(vector_op))
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
    error[i] <- dot_integrate(abs(X[,i] - estimate[,i]))
  }

  sum(error)
}
