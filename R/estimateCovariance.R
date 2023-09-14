
#' Title
#'
#' @param X
#' @param lag
#'
#' @return
#' @export
#'
#' @examples
.estim_lag_cov <- function(X, lag=0){
  # Build from (3.9)/(3.11)
  N <- ncol(X)

  m1 <- (1/(N-lag)) * rowSums(X[,1:(N-lag)])

  covMat <- matrix(0, nrow=nrow(X),ncol=nrow(X))

  for(k in 1:(N-lag)){
    covMat <- covMat + ( (X[,k]-m1) %*% t(X[,k+lag]-m1) )
  }

  covMat /((N-lag)-1)
}
