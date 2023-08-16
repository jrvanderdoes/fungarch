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
pointwiseARCH <- function(X, w=0.001,
                          lambda=c(10^-7,10^-6,10^-5,10^-4,
                                   seq(0.001,0.01,0.002),
                                   seq(0.02,0.09,0.01),
                                   seq(0.1,5,0.1),seq(6,20,1),
                                   seq(22,50,2),
                                   seq(60,100,10),seq(125,250,25),
                                   seq(300,500,50))){
  X <- X^2 # Pointwise square
  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_pointwise_lag_cov(X)
  C_1 <- .estim_pointwise_lag_cov(X, lag = 1)

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

  # Estimate delta, (4.9) / (3.30)
  hat_delta <- m2 - (hat_alpha %*% m2) / n
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
#' @param w
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
pointwiseGARCH <- function(X, w=0.001,
                           lambda=c(10^-7,10^-6,10^-5,10^-4,
                                    seq(0.001,0.01,0.002),
                                    seq(0.02,0.09,0.01),
                                    seq(0.1,5,0.1),seq(6,20,1),
                                    seq(22,50,2),
                                    seq(60,100,10),seq(125,250,25),
                                    seq(300,500,50))){
  X <- X^2 # Pointwise square
  N <- ncol(X) # Number of observations
  n <- nrow(X) # Resolution

  # Lag-h Covariance operator estimates
  #   Formula (3.6)
  C_0 <- .estim_pointwise_lag_cov(X)
  C_1 <- .estim_pointwise_lag_cov(X, lag = 1)

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

  # Prep
  I_H <- diag(n)

  # Define m2, (3.3)
  m2 <- 1 / N * rowSums(X^2)

  # Estimate alpha and delta
  error_plot <- NULL
  if(length(lambda)>1){
    lambda_result <- .findLambdaGARCH(lambda,X=X,C_0=C_0,C_1=C_1,I_H=I_H,N=N,
                                      projector=projector,m2=m2)
    error_plot <- lambda_result$error_plot
    lambda <- lambda_result$lambda
  }

  ## Estimate Pis (3.42)
  GL_hat
  GL1_hat

  invMat <- solve( sqrt(N) * C_0 %*% C_0 + lambda * I_H )
  #   Formulas (3.39), (3.42), and (3.43)
  hat_pi[[i]] <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector


  hat_beta <- 0
  hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector

  # Estimate delta, (4.9) / (3.30)
  hat_delta <- (m2 - (hat_alpha + hat_beta) %*% m2) / n

  list('hat_alpha'=hat_alpha,
       'hat_beta'=hat_beta,
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
#' @param alpha
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
.computeError <- function(X,alpha,delta){
  error <- rep(0, ncol(X))
  estimate <- .estimateARCH(X,alpha,delta)

  # for(i in 2:ncol(X)){
  #   for(j in 1:nrow(X)){
  #     vector_op <- alpha[j, ] * X[,(i-1)] # X^2
  #     estimate[j] <- delta[j] + #.approx_int(as.matrix(vector_op))
  #       pracma::trapz(x=seq(0,1,length.out=length(vector_op)), y=vector_op)
  #   }
  #
  #   #error[i] <- sum((X[,i] - estimate)^2)
  #   #error[i] <- sqrt(.approx_int((X[,i] - estimate)^2))
  #   error[i] <- pracma::trapz((X[,i] - estimate)^2)
  # }
  for(i in 2:ncol(X)){
    error[i] <- pracma::trapz((X[,i] - estimate^2)^2)
  }


  sum(error)
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
.estim_pointwise_lag_cov <- function(X, lag=0){
  N <- ncol(X)

  # TODO:: Confirm if use definition from (3.3) or after (3.6)?
  m1 <- (1/(N-lag)) * rowSums(X[,1:(N-lag)])

  covMat <- matrix(0, nrow=nrow(X),ncol=nrow(X))

  for(k in 1:(N-lag)){
    covMat <- covMat + ( (X[,k]-m1) %*% t(X[,k+lag]-m1) )
  }

  covMat /((N-lag)-1)
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
.findLambdaARCH <- function(lambdas,X,C_0,C_1,I_H,N,projector,m2){
  errors <- rep(NA,length(lambdas))

  for(i in 1:length(lambdas)){
    # Formula (3.26)
    invMat <- solve( sqrt(N) * C_0 %*% C_0 + lambdas[i] * I_H )
    #   Formula (4.8)
    hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector

    # Estimate delta, (4.9)
    hat_delta <- m2 - (hat_alpha %*% m2) / nrow(X)

    errors[i] <- .computeError(X = X, alpha = hat_alpha, delta = hat_delta)
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
#' @param lambdas
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
.findLambdaGARCH <- function(lambdas,X,C_0,C_1,I_H,N,projector,m2){
  errors <- rep(NA,length(lambdas))

  for(i in 1:length(lambdas)){
    # Formula (3.26)
    invMat <- solve( sqrt(N) * C_0 %*% C_0 + lambdas[i] * I_H )
    #   Formula (4.8)
    hat_alpha <- sqrt(N) * C_1 %*% C_0 %*% invMat %*% projector

    # Estimate delta, (4.9)
    hat_delta <- (m2 - (hat_alpha + hat_beta) %*% m2) / nrow(X)

    errors[i] <- .computeError(X = X, alpha = hat_alpha, delta = hat_delta)
  }

  tmp <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(lambdas,errors)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw()
  list('lambda'=lambdas[which.min(errors)],
       'error_plot'= tmp)
}
