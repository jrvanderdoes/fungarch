
#' OU Covariance
#'
#' Generate matrix with OU covariance
#'
#' @param n Positive integer for the dimension of the covariance matrix.
#' @param denom Positive finite decay scale used in the exponent.
#'
#' @return Matrix with OU covariance.
#'
#' @example
#' fgarch::.generate_ceps(10)
#'
#' @noRd
.generate_ceps <- function(n, denom=2){
  if (length(n) != 1L || !is.numeric(n) || is.na(n) ||
      n < 1 || n != floor(n)) {
    stop("n must be a single positive integer", call. = FALSE)
  }
  if (length(denom) != 1L || !is.numeric(denom) || is.na(denom) ||
      !is.finite(denom) || denom <= 0) {
    stop("denom must be a single positive finite number", call. = FALSE)
  }

  n <- as.integer(n)
  C_eps <- matrix(NA, n, n)
  times <- seq_len(n) / n

  for (i in 1:n){
    C_eps[i,] <- exp(-abs(times[i]-times)/denom) # min(times[i], times)
  }

  C_eps
}
