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

.computeARCHError_tensor_new <- function(X, alpha, delta, C_eps){
  error <- rep(0, dim(X)[3])
  estimate <- .predictARCH_tensor_new(X,alpha,delta, C_eps)

  for(i in 2:dim(X)[3]){
    error[i] <- .approx_int(abs(X[,,i] - estimate[,,i]))
  }

  sum(error)
}

.predictARCH_tensor_new <- function(X2, alpha, delta, C_eps){
  dims <- dim(X2)
  Sigma <- array(NA, dims)
  #delta_l <- stack(as.data.frame(delta))$values

  #hd <- matrix(ests$hat_delta,nrow=4,ncol=3)
  #tmp = matrix(delta + alpha %*% X2[,1],nrow = sqrt(n),ncol = sqrt(n))
  #matrix(ests$hat_delta + ests$hat_alpha %*% X2[,1],nrow = sqrt(n),ncol = sqrt(n))
  Sigma[,,1] <- delta #+ alpha %*% X2[,1]
  # c_e is cov OU

  ## TODO:: Review this
  tmp <- 0
  for(i in 2:dims[3]){
    # TODO:: Change to matrix mult
    for(j in 1:dims[1]){
      ## TODO:: Discuss this...
      vector_op <- alpha[j, ] * X2[,,(i-1)] # X^2
      Sigma[j,,i] <- delta[j,] + dot_integrate_col(vector_op)
    }
  }

  #Sigma <- ifelse(Sigma>0,Sigma,0)

  Estimate <- array(NA, dims)
  for(i in 1:dims[3]){
    Sigma_svd <- svd(Sigma[,,i])
    Sigma_sqrt <- Sigma_svd$u %*% diag(sqrt(Sigma_svd$d)) %*% t(Sigma_svd$v)
    # Make sure to add sum of Est[,i-1]

    Estimate[,,i] <- Sigma_sqrt %*% (C_eps) %*% Sigma_sqrt
  }

  Estimate
}

###############################################

#' Predict GARCH Estimates
#'
#' This (internal) function estimates given X^2, alpha, beta, and delta
#'  of an GARCH(1,1) process.
#'
#' @inheritParams .predictARCH
#'
#' @return Estimate of GARCH(1,1) based on data with given alpha, beta, and delta
#'
#' @noRd
.predictGARCH <- function(X2, alpha, beta, delta, est_sqr=FALSE){
  estimate <- matrix(0, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- delta

  for(i in 2:ncol(X2)){
    for(j in 1:nrow(X2)){
      vector_for_alpha_op <- alpha[j, ] * X2[,(i-1)]  # X^2
      tmp <- ifelse(X2[,(i-1)]>0,X2[,(i-1)],0)
      vector_for_beta_op <- beta[j, ] * tmp  # X

      estimate[j,i] <- delta[j] +
        dot_integrate(as.matrix(vector_for_alpha_op)) +
        dot_integrate(as.matrix(vector_for_beta_op))
    }
  }

  if(est_sqr) return(estimate)

  estimate <- ifelse(estimate>0,estimate,0)

  sqrt(estimate)
}


.predictGARCH <- function(X2, alpha, beta, delta, est_sqr=FALSE){
  estimate <- matrix(0, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- delta

  for(i in 2:ncol(X2)){
    for(j in 1:nrow(X2)){
      vector_for_alpha_op <- alpha[j, ] * X2[,(i-1)]  # X^2
      tmp <- ifelse(X2[,(i-1)]>0,X2[,(i-1)],0)
      vector_for_beta_op <- beta[j, ] * tmp  # X

      estimate[j,i] <- delta[j] +
        dot_integrate(as.matrix(vector_for_alpha_op)) +
        dot_integrate(as.matrix(vector_for_beta_op))
    }
  }

  if(est_sqr) return(estimate)

  estimate <- ifelse(estimate>0,estimate,0)

  sqrt(estimate)
}

#' Compute GARCH Error
#'
#' This (internal) function estimates values and then compares the estimates
#'  to the known values of a GARCH(1,1) process. The functional L2 error is
#'  estimated.
#'
#' @inheritParams .computeARCHError
#'
#' @return  Estimate error based on data with given alpha, beta, and delta
#'
#' @noRd
.computeGARCHError <- function(X,alpha,beta,delta){
  error <- rep(0, ncol(X))
  estimate <- .predictGARCH(X,alpha,beta,delta,est_sqr=TRUE)

  for(i in 2:ncol(X)){
    error[i] <- sqrt(dot_integrate((sqrt(X[,i]) - estimate[,i])^2))
  }

  sum(error)
}
