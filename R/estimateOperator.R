#' Operator ARCH Estimation
#'
#' Estimates the parameters for ARCH(1) model using operator Yule Walker
#'  estimation.
#'
#' @inheritParams pointwiseARCH
#'
#' @return
#' @export
#'
#' @examples
operatorARCH <- function(X, w=0.001, theta=NA){
  tmp <- X - rowMeans(X) # De-mean
  X2 <- array(apply(tmp, 2, tcrossprod),
             dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))
  #X2 <- apply(X2,3, function(XX){array(XX)})

  N <- ncol(X)
  n <- nrow(X)

  eps <- .generateFunctionalIID(N,100)

  # C_eps
  resolution <- nrow(X)
  C_eps <- matrix(NA,resolution,resolution)
  times <- 1:resolution/resolution
  for (i in 1:resolution){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }
  #C_eps_long <- stack(data.frame(C_eps))$values

  C_X <- rowMeans(X2,dims=2)
  # Define Y
  Y <- Y_eps <- X2
  for(i in 1:N){
    Y[,,i] <- Y[,,i] - C_X
    Y_eps[,,i] <- Y[,,i] %*% C_eps
  }
  # C_eps %*% Y
  # Y_eps <-  apply(Y, MARGIN=2,#apply(Y, MARGIN=3,
  #                     FUN = function(y,ce){y %*% ce},
  #                     ce=C_eps, simplify = T)
  # Y_eps <- array(NA,dim=c(nrow(X),nrow(X),ncol(X)))
  # for(k in 1:dim(Y_eps)[3]){
  #   Y_eps[,,k] <- tmp_Y_eps[[k]]
  # }

  # Switch Y to long
  #Y <- apply(Y,3, function(XX){array(XX)})
  #Y_eps <- apply(Y_eps,3, function(XX){array(XX)})

  # Define Z (45 / 46)
  #   TODO:: Deal with larger p
  Z <- Y # stack(data.frame(Y))$values
  Z_eps <- Y_eps

  # Lag-h Covariance operator estimates
  C_ZZe0 <- .estim_lag_cov_tensor(Z, Z_eps, lag=0)
  C_ZY1 <- .estim_lag_cov_tensor(Z, Y, lag=1)

  eigval_C_ZZe0 <- eigen(C_ZZe0)
  eigval_C_eps <- eigen(C_eps)

  # Select K Projections
  #   Remark 3.2(c)
  #   TODO:: Check Imaginary Size
  K <- which.min(cumsum(eigval_C_ZZe0$values)/sum(eigval_C_ZZe0$values) < 1-w)
  K1 <- which.min(cumsum(eigval_C_eps$values)/sum(eigval_C_eps$values) < 1-w)

  # Standardize Vectors in integral
  #   apply( eigval_C_ZZe0$vectors^2,MARGIN = 2,.approx_int)
  eigval_C_ZZe0$vectors <- Re(eigval_C_ZZe0$vectors) * sqrt(n)
  eigval_C_eps$vectors <- Re(eigval_C_eps$vectors) * sqrt(sqrt(n))

  # Compute Projection
  projector <- projector_eps <-matrix(0, nrow=n, ncol=n)
  # matrix(0, nrow=sqrt(n), ncol=sqrt(n))
  for(i in 1:K){
    projector <- projector +
      ( eigval_C_ZZe0$vectors[,i] %*% t(eigval_C_ZZe0$vectors[,i]) )
  }
  for(i in 1:K1){
    projector_eps <- projector_eps +
      ( eigval_C_eps$vectors[,i] %*% t(eigval_C_eps$vectors[,i]) )
  }

  # Define m2, (3.3)
  #     TODO:: Check this out (before 50)
  #C_X <- 1 / N * rowSums(Y)

  # Estimate alpha and delta
  error_plot <- NULL
  if(is.na(theta)){
    optim_results <-
      optim(10000,function(theta, X, C_ZZe0, C_ZY1, projector, C_X,
                           C_eps, projector_eps, n){
        ests <- .estimateARCH_tensor(
          C_ZZe0=C_ZZe0, C_ZY1=C_ZY1, theta=theta,
          projector=projector, C_X=C_X, C_eps=C_eps,
          projector_eps=projector_eps, n=n)

        .computeARCHError_tensor_new(
          X = X2, alpha = ests$hat_alpha, delta = ests$hat_delta,
          C_eps=C_eps)
      },
      method="Brent",
      lower=10^-7,upper=500000,
      X=X,C_ZZe0=C_ZZe0,C_ZY1=C_ZY1,projector=projector,
      C_X=C_X, C_eps=C_eps, projector_eps=projector_eps, n=n)

    theta <- optim_results$par
  }

  ests <- .estimateARCH_tensor(
    C_ZZe0=C_ZZe0, C_ZY1=C_ZY1, theta=theta,
    projector=projector, C_X=C_X, C_eps=C_eps,
    projector_eps=projector_eps, n=n)

  list('hat_alpha'=ests$hat_alpha,
       'hat_delta'=ests$hat_delta,
       # 'estimate'=.predictARCH_tensor(X2=X, alpha=ests$hat_alpha,
       #                         delta=ests$hat_delta, est_sqr=FALSE),
       'error'=.computeARCHError_tensor_new(
            X = Z, alpha = ests$hat_alpha, delta = ests$hat_delta, C_eps = C_eps),
       'fitting_params'=list('theta'=theta,
                             'C_ZZe0'=C_ZZe0, 'C_ZY1'=C_ZY1,
                             'eigval_C_ZZe0'=eigval_C_ZZe0,
                             'eigval_C_eps'=eigval_C_eps,
                             'K'=K, 'K1'=K1,
                             'projector'=projector,
                             'projector_eps'=projector_eps,
                             'C_X'=C_X, 'C_eps'=C_eps) )
}


#' Title
#'
#' @param C_ZZe0
#' @param C_ZY1
#' @param theta
#' @param projector
#' @param C_X
#' @param C_eps
#' @param projector_eps
#' @param n
#'
#' @return
#' @export
#'
#' @examples
.estimateARCH_tensor <- function(C_ZZe0, C_ZY1, theta, projector,
                                 C_X, C_eps, projector_eps,
                                 n=nrow(C_ZZe0)){
  # Formula (49)
  invMat <- solve( C_ZZe0 %*% t(C_ZZe0) + diag(theta, n, n) )
  hat_alpha <- C_ZY1 %*% t(C_ZZe0) %*% invMat %*% projector

  # Estimate delta, (51)
  C_eps_d <- C_eps %*% solve(C_eps %*% t(C_eps) + diag(theta, n, n))
  # t(replicate(length(C_X),C_X))
  hat_delta <- ( C_X - hat_alpha %*% t(C_X) %*% C_eps ) %*%
    C_eps_d %*% projector_eps

  list('hat_alpha'=hat_alpha,
       'hat_delta'=hat_delta)
}


#'
#' #' Tensor
#' #'
#' #' This (interal) function computes tensor product.
#' #'
#' #' CURRENTLY UNUSED!
#' #'
#' #' @param X
#' #' @param Y
#' #'
#' #' @return
#' #'
#' #' @examples
#' .tensor <- function(X,Y=X){
#'   ## TODO:: Discuss this
#'   if(is.vector(X) && is.vector(Y)){
#'     result <- X %*% t(Y)
#'   } else{
#'     result <- list()
#'     for(i in 1:ncol(X)){
#'       result[[i]] <- X[,i] %*% t(Y[,i])
#'     }
#'     # Create Tensor Matrix
#'     tensorVal <- do.call("cbind",result)
#'     for(i in 2:ncol(X)){
#'       # Build row
#'       row <- result[[i]]
#'       for(j in 2:ncol(X)){
#'         func <- ifelse(j<i,t,function(X){X})
#'         row <- cbind(row,func(result[[abs(i-j)+1]]))
#'       }
#'       tensorVal <- rbind(tensorVal, row)
#'     }
#'   }
#'
#'   tensorVal
#' }
