#' Pointwise Yule Walker ARCH Estimation
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
pointwiseARCH <- function(X, w=0.001, lambda=NA){
  #X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square
  #X <- X - rowMeans(X) # De-mean
  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
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
  m2 <- 1 / N * rowSums(X^2)

  # Estimate alpha and delta
  error_plot <- NULL
  if(length(lambda)>1){
    lambda_result <- .findLambdaARCH(lambda,X=X,C_0=C_0,C_1=C_1,#I_H=I_H, N=N,
                                     projector=projector,m2=m2)
    error_plot <- lambda_result$error_plot
    lambda <- lambda_result$lambda
  } else if(is.na(lambda)){
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
                             'error_plot'=error_plot))
}


#' Find Lambda Values from Various Options
#'
#' XXXX
#'
#' @param X
#' @param C_0
#' @param C_1
#' @param I_H
#' @param N
#' @param projector
#' @param m2
#'
#' @return
#' @export
#'
#' @examples
.findLambdaARCH <- function(lambdas,X,C_0,C_1,I_H,projector,m2){
  errors <- rep(NA,length(lambdas))

  for(i in 1:length(lambdas)){
    ests <- .estimateARCH(N=ncol(X), n=nrow(X), C_0=C_0, C_1=C_1,
                  lambda=lambdas[i], projector=projector, m2=m2)

    errors[i] <- .computeARCHError(X = X,
                               alpha = ests$hat_alpha,
                               delta = ests$hat_delta)
  }

  tmp <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(lambdas,errors)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw()

  list('lambda'=lambdas[which.min(errors)],
       'error_plot'= tmp)
}


#' Title
#'
#' @param N
#' @param n
#' @param c_0
#' @param lambda
#' @param projector
#' @param m2
#'
#' @return
#' @export
#'
#' @examples
.estimateARCH <- function(N, n, C_0, C_1, lambda, projector, m2){
  # Formula (3.26)
  invMat <- solve( sqrt(N) * C_0 %*% C_0 + diag(lambda, n, n) )
  #   Formula (4.8)
  hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector

  # Estimate delta, (4.9)
  hat_delta <- m2 - (hat_alpha %*% m2) / n

  list('hat_alpha'=hat_alpha,
       'hat_delta'=hat_delta)
}



################################################################################



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
pointwiseGARCH <- function(X, w=0.001, lambdas=NA, L=3){
  # Setup Data
  X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square
  # X <- X - rowMeans(X) # De-mean

  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X^2)

  ## Estimate G_hat (3.9)
  C <- list()
  for(l in 0:L){
    C[[l+1]] <- .estim_lag_cov(X,lag=l)
    # t(matrix(c(l-1,l,l+1,
    #            l,l+1,l+2,
    #            l+1,l+2,l+3),3,3))
  }
  GL_hat <- do.call("cbind",C[-(L+1)])
  GL1_hat <- do.call("cbind",C[-1])
  for(l in 2:L){
    # Build row
    row <- t(C[[l-1+1]])
    row1 <- t(C[[l-1+2]])
    for(i in 2:L){
      func <- ifelse(i<l,t,function(X){X})
      row <- cbind(row,func(C[[abs(l-i)+1]]))
      row1 <- cbind(row1,func(C[[abs(l-i)+2]]))
    }
    GL_hat <- rbind(GL_hat,row)
    GL1_hat <- rbind(GL1_hat,row1)
  }

  # Projection
  eigval_GL <- eigen(GL_hat)

  # Select K Projections
  #   Remark 3.2(c)
  K <- which.min(cumsum(eigval_GL$values)/sum(eigval_GL$values) < 1-w)

  eigval_GL$vectors <- eigval_GL$vectors * sqrt(n) #TODO:: check

  # Compute Projection
  projector_pi <- matrix(0, nrow=n*L, ncol=n*L)
  for(i in 1:K){
    projector_pi <- projector_pi +
      ( eigval_GL$vectors[,i] %*% t(eigval_GL$vectors[,i]) )
  }

  # Two param optimization, lambda_pi and lambda_beta
  if(is.na(lambdas)){
    optim_results <- optim(c(3,10),
                           fn = function(lambdas,GL_hat, GL1_hat,
                                         projector_pi,
                                         N, n, m2, L, w){
                             ests <- .estimateGARCH(lambdas,
                                                    GL_hat, GL1_hat,
                                                    projector_pi,
                                                    N, n, m2, L, w)
                             .computeGARCHError(X = X, alpha = ests$hat_alpha,
                                                beta = ests$hat_beta,
                                                delta = ests$hat_delta)
                             },
                           GL_hat=GL_hat, GL1_hat=GL1_hat,
                           projector_pi=projector_pi,
                           N=N, n=n, m2=m2, L=L, w=w)
    lambdas <- optim_results$par
  }else if(length(lambdas)!= 2){
    stop('Errors: lambdas must be NA or length 2')
  }


  ests <- .estimateGARCH(lambdas,
                         GL_hat, GL1_hat,
                         projector_pi,
                         N, n, m2, L, w)

  list('hat_alpha'=ests$hat_alpha,
       'hat_beta'=ests$hat_beta,
       'hat_delta'=ests$hat_delta,
       'fitting_params'=list('lambdas'=lambdas,
                             'GL_hat'=GL_hat, 'GL1_hat'=GL1_hat,
                             'eigval_GL'=eigval_GL,
                             'projector_pi'=projector_pi,
                             'K'=K,
                             'm2'=m2))
}


#' Title
#'
#' @param lambda_pi
#' @param lambda_beta
#' @param GL_hat
#' @param GL1_hat
#' @param projector_pi
#' @param projector_beta
#' @param N
#' @param n
#'
#' @return
#' @export
#'
#' @examples
.estimateGARCH <- function(lambdas,
                           GL_hat, GL1_hat,
                           projector_pi,
                           N, n, m2, L, w){
  lambda_pi <- lambdas[1]
  lambda_beta <- lambdas[2]

  ## Estimate Pis (3.42)
  invMat_pi <- solve( sqrt(N) * GL_hat %*% GL_hat + diag(lambda_pi,L*n,L*n) )
  pi_hat <- sqrt(N) * GL1_hat %*% GL_hat %*% invMat_pi %*% projector_pi

  pis <- list()
  for(i in 0:(L-1)){
    pis[[i+1]] <- pi_hat[1:n,i*n+1:n]
  }
  # Formula (3.41)
  bigPI <- pis[[1]]
  boldPI <- pis[[2]]


  ## Alpha (3.43)
  hat_alpha <- pis[[1]]

  ## Find projector_beta
  eigval_pi <- eigen(bigPI)

  # Select K Projections, Remark 3.2(c)
  K <- which.min(Re(cumsum(eigval_pi$values)/sum(eigval_pi$values)) < 1-w)
  eigval_pi$vectors <- eigval_pi$vectors * sqrt(n)

  projector_beta <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K){
    projector_beta <- projector_beta +
      ( eigval_pi$vectors[,i] %*% t(eigval_pi$vectors[,i]) )
  }


  ## Estimate Beta (3.39)/(3.43)
  invMat_beta <- solve( sqrt(N) * bigPI %*% t(bigPI) +
                          diag(lambda_beta, n, n) )
  hat_beta <- sqrt(N) * boldPI %*% t(bigPI) %*%
    invMat_beta %*% projector_beta

  # Estimate delta, (3.30)
  hat_delta <- m2 - ((hat_alpha + hat_beta) %*% m2) / n

  list('hat_alpha'=(hat_alpha),
       'hat_beta'=Re(hat_beta),
       'hat_delta'=Re(hat_delta))
}
