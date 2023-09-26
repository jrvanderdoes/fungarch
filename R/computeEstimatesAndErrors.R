#' Predict ARCH Estimates
#'
#' This (internal) function estimates given X^2, alpha, and delta
#'  of an ARCH(1) process.
#'
#' @param X2 Matrix of numerics indicating the squared process
#' @param alpha Matrix of numerics estimating the alpha values
#' @param delta Matrix of numerics estimating the delta values
#'
#' @return
.predictARCH <- function(X2, alpha, delta, est_sqr=FALSE){
  estimate <- matrix(0, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- delta

  for(i in 2:ncol(X2)){
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * X2[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + .approx_int(as.matrix(vector_op))
      #pracma::trapz(x=seq(0,1,length.out=length(vector_op)), y=vector_op)
    }
  }

  if(est_sqr) return(estimate)

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
#' @return
.computeARCHError <- function(X, alpha, delta){
  error <- rep(0, ncol(X))
  estimate <- .predictARCH(X,alpha,delta,est_sqr=TRUE)

  for(i in 2:ncol(X)){
    # TODO:: Ask about this
    error[i] <- sqrt(.approx_int((sqrt(X[,i]) - estimate[,i])^2))
    #error[i] <- sqrt(.approx_int((X[,i] - estimate[,i])^2))
  }

  sum(error)
}


###############################################

#' Predict GARCH Estimates
#'
#' This (internal) function estimates given X^2, alpha, beta, and delta
#'  of an GARCH(1,1) process.
#'
#' @inheritParams .predictARCH
#'
#' @return
.predictGARCH <- function(X2, alpha, beta, delta, est_sqr=FALSE){
  estimate <- matrix(0, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- delta

  for(i in 2:ncol(X2)){
    for(j in 1:nrow(X2)){
      vector_for_alpha_op <- alpha[j, ] * X2[,(i-1)]  # X^2
      tmp <- ifelse(X2[,(i-1)]>0,X2[,(i-1)],0)
      vector_for_beta_op <- beta[j, ] * tmp  # X

      estimate[j,i] <- delta[j] +
        .approx_int(as.matrix(vector_for_alpha_op)) +
        .approx_int(as.matrix(vector_for_beta_op))
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
#' @return
.computeGARCHError <- function(X,alpha,beta,delta){
  error <- rep(0, ncol(X))
  estimate <- .predictGARCH(X,alpha,beta,delta,est_sqr=TRUE)

  for(i in 2:ncol(X)){
    error[i] <- sqrt(.approx_int((sqrt(X[,i]) - estimate[,i])^2))
    #error[i] <- sqrt(.approx_int(((X[,i]) - estimate[,i])^2))
  }

  sum(error)
}
