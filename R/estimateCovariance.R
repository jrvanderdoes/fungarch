
#' Title
#'
#' @param X
#' @param lag
#'
#' @return
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


#' Title
#'
#' @param X
#' @param lag
#'
#' @return
#'
#' @examples
.estim_lag_cov_tensor <- function(X, Y, lag=0){
  if(prod(dim(X) != dim(Y)))
    stop('Error: X and Y have different dimensions')
  # Build from (3.9)/(3.11)
  N <- dim(X)[3]

  covMat <- matrix(0, nrow=nrow(X),ncol=nrow(X))

  for(k in 1:(N-lag)){
    covMat <- covMat + ( X[,,k] %*% t(Y[,,k+lag]) )
  }

  # var.tensor(X[[b=1]],X[[a=~"a'",b=2]],along="i")
  covMat / (N-lag) # TODO:: Check
}
