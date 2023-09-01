#' Pointwise Yule Walker ARCH Estimation
#'
#' TODO:: Add description
#'
#' @param X Functional data, curves are the columns
#' @param c (Optional) Numeric. Preliminary factor (at least 10^6)
#' @param w (Optional) Numeric. Weight for choice of K.
#'
#' @return
#' @export
#'
#' @references Functional ARCH and GARCH models: A Yule-Walker approach.
#'  Sebastian K\"uhnert (TODO:: Complete reference)
#'
#' @examples
#' set.seed(12345)
#' # alpha <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
#' # delta <- function(t){ 0.1* t * (1-t) }
#' alpha <- function(s,t){(1/8)*((s-1/2)^2 + (t-1/2)^2)}
#' delta <- function(t){ rep(0.1,length(t)) }
#' dat <- generateFunctionalCondVol(51,500,'arch',
#'   alpha_par = alpha, delta_par = delta)
#' est <- estimateYuleWalkerARCH(X = dat$garch_mat[-nrow(dat$garch_mat),])
#'
#' s <- seq (0,1,length.out=nrow(dat$garch_mat))
#' alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
#' for(s1 in 1:length(s)){
#'   for(s2 in 1:length(s)){
#'     alpha_true[s1,s2] <- alpha(s[s1],s[s2])
#'   }
#' }
#' # Compare Alpha
#' plot3D::persp3D(z=alpha_true,main=parse(text = paste0('~ z == a(s,t)')),
#'                 theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
#' plot3D::persp3D(z=est$hat_alpha,main=bquote(~z == hat(a)(s,t)),
#'                 theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
#' # Compare Delta
#' d <- seq(0,1,length.out=nrow(dat$garch_mat))
#' delta_true <- delta(d)
#' plot(x=d,delta_true,
#'      ylim=c(min(delta_true, est$hat_delta),
#'             max(delta_true,est$hat_delta)),
#'      type='l')
#' lines(x=d,y=est$hat_delta, col='red')
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


#' Title
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
pointwiseGARCH <- function(X, w=0.001, lambda=NA, L=3){
  #X <- X - rowMeans(X) # De-mean
  X <- X^2 # Pointwise square
  #X <- X - rowMeans(X) # De-mean
  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # # Lag-h Covariance operator estimates
  # #   Formula (3.6)
  # C_0 <- .estim_lag_cov(X)
  # C_1 <- .estim_lag_cov(X, lag = 1)
  #
  # eigval_C_0 <- eigen(C_0)
  #
  # # Select K Projections, Remark 3.2(c)
  # K <- which.min(cumsum(eigval_C_0$values)/sum(eigval_C_0$values) < 1-w)
  #
  # # Standardize Vectors in integral
  # #   apply( eigval_C_0$vectors^2,MARGIN = 2,.approx_int)
  # eigval_C_0$vectors <- eigval_C_0$vectors * sqrt(n)
  #
  # # Compute Projection
  # projector <- matrix(0, nrow=n, ncol=n)
  # for(i in 1:K){
  #   projector <- projector +
  #     ( eigval_C_0$vectors[,i] %*% t(eigval_C_0$vectors[,i]) )
  # }
  #
  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X^2)

  ## Estimate G_hat (3.9)
  C <- list()
  for(l in 0:L){
    C[[l+1]] <- .estim_lag_cov(X,lag=l)
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
  optim_results <- optim(c(1000,1000),
                         fn = function(lambdas){
                           ests <- .estimateGARCH()
                           .computeGARCHError(X = X, alpha = ests$hat_alpha,
                                         beta = ests$hat_beta,
                                         delta = ests$hat_delta)})

  ests <- .estimateGARCH(optim_results$par[1], optim_results$par[2],
                         GL_hat, GL1_hat,
                         projector_pi, projector_beta,
                         N, n, m2)

  list('hat_alpha'=ests$hat_alpha,
       'hat_beta'=ests$hat_beta,
       'hat_delta'=ests$hat_delta,
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
.estimateGARCH <- function(lambda_pi, lambda_beta,
                           GL_hat, GL1_hat,
                           projector_pi, projector_beta,
                           N, n, m2){
  ## Estimate Pis (3.42)
  invMat_pi <- solve( sqrt(N) * GL_hat %*% GL_hat +
                        diag(lambda_pi,L*n,L*n) )
  pi_hat <- sqrt(N) * GL1_hat %*% GL_hat %*% invMat_pi %*% projector_pi
  for(i in 0:L){
    pis[[i]] <- pi_hat[1:n,i*n+1:n]
  }
  # Formula (3.41)
  #   Discuss
  bigPI <- pis[[1]]
  boldPI <- pis[[1]]


  ## Alpha (3.43)
  hat_alpha <- pis[[1]]

  ## Find projector_beta
  eigval_pi <- eigen(pi_hat)

  projector_beta <- matrix(0, nrow=n*L, ncol=n*L)
  for(i in 1:K){
    projector_beta <- projector_beta +
      ( eigval_pi$vectors[,i] %*% t(eigval_pi$vectors[,i]) )
  }


  ## Estimate Beta (3.39)/(3.43)
  #   Ask pi^*?
  invMat_beta <- solve( sqrt(N) * bigPI %*% bigPI +
                          diag(lambda_beta, n, n) )
  hat_beta <- sqrt(N) * boldPI %*% bigPI %*%
    invMat_beta %*% projector_beta

  # Estimate delta, (3.30)
  hat_delta <- m2 - ((hat_alpha + hat_beta) %*% m2) / n

  list('hat_alpha'=hat_alpha,
       'hat_beta'=hat_beta,
       'hat_delta'=hat_delta)
}
