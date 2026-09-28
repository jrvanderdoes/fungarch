
#' Estimate a lagged covariance operator
#'
#' Estimate the lag-`lag` covariance matrix for functional observations
#' stored in the columns of `X`. At lag 0 this is the usual sample
#' covariance; for a positive lag it estimates the cross-covariance between
#' observations separated by that lag.
#'
#' @param X Numeric matrix. Rows are evaluation points and columns are
#'   observations.
#' @param lag Nonnegative integer lag. The lag must leave at least two paired
#'   observations.
#'
#' @return A numeric square matrix with `nrow(X)` rows and columns.
#'
#' @details The centering vector is estimated from the first `N - lag`
#'   observations and is used for both members of each lagged pair. The
#'   positive-lag estimate need not be symmetric. Missing values are not
#'   handled.
#'
#' @noRd
.estim_lag_cov <- function(X, lag=0){
  # Build from (3.9)/(3.11)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix", call. = FALSE)
  }
  if (length(lag) != 1L || !is.numeric(lag) || is.na(lag) ||
      lag < 0 || lag != floor(lag)) {
    stop("lag must be a single nonnegative integer", call. = FALSE)
  }

  N <- ncol(X)
  lag <- as.integer(lag)
  n_pairs <- N - lag
  if (n_pairs < 2L) {
    stop("X must contain at least two observations after applying lag", call. = FALSE)
  }
  
  m1 <- rowMeans(X[, seq_len(n_pairs), drop = FALSE])
  
  covMat <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  for (k in seq_len(n_pairs)) {
    covMat <- covMat + ((X[, k] - m1) %*% t(X[, k + lag] - m1))
  }
  
  covMat / (n_pairs - 1)
}
