#' Apply a discretized integral operator
#'
#' Multiply a matrix operator by the vectorized input and divide by the
#' number of grid points, preserving the package's equal-grid integration
#' scaling.
#'
#' @param X Numeric operator matrix.
#' @param Y Numeric vector or matrix to vectorize before applying `X`.
#'
#' @return A numeric vector containing the integrated operator result.
#'
#' @noRd
operator <- function(X, Y){
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix", call. = FALSE)
  }
  if (!is.numeric(Y) || any(!is.finite(Y))) {
    stop("Y must contain finite numeric values", call. = FALSE)
  }
  if (nrow(X) < 1L || ncol(X) != length(Y)) {
    stop("X columns must match the vectorized length of Y", call. = FALSE)
  }

  X %*% as.vector(Y) / nrow(X)
}

#' Vectorize a tensor product
#'
#' Construct the outer product of two numeric inputs and return it as a row
#' matrix. The vector order matches the previous implementation after
#' `as.vector()` is applied.
#'
#' @param X Numeric vector or matrix.
#' @param Y Numeric vector or matrix.
#'
#' @return A `1 x (length(X) * length(Y))` numeric matrix.
#'
#' @noRd
tensor_vec <- function(X,Y){
  if (!is.numeric(X) || any(!is.finite(X)) ||
      !is.numeric(Y) || any(!is.finite(Y))) {
    stop("X and Y must contain finite numeric values", call. = FALSE)
  }

  # OLD: t(as.vector(X) %o% as.vector(Y))
  matrix(
    as.vector(t(as.vector(X) %o% as.vector(Y))),
    nrow = 1L
  )
  # kronecker(as.vector(X),t(as.vector(Y)))
}
