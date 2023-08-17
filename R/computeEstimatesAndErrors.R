#' Title
#'
#' @param X2
#' @param alpha
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
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
  estimate <- .predictARCH(X,alpha,delta,est_sqr=TRUE)

  for(i in 2:ncol(X)){
    error[i] <- sqrt(.approx_int((sqrt(X[,i]) - estimate[,i])^2))
  }

  sum(error)
}
